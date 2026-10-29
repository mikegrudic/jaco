"""A model as an explicit ODE right-hand side for a generic integrator (scipy.integrate.solve_ivp or any other).

    sys = OdeSystem(make_model())
    rates = sys.rates(T, x, params)             # heat [erg cm^-3 s^-1], dx/dt [per H per s], Td, u, P, rho, ...
    dydt = sys.rhs(y, params)                   # y = (u, x_1, ..., x_m) in sys.integrated order; du/dt, dx/dt
    T, x = sys.temperature(u, x, params)        # the gas temperature and every abundance behind a state vector

The system is the model's reduced network as the code generator sees it: the steady-state closures, the conservation
reductions, the derived parameters and the intermediates are applied; the backward-Euler terms are removed by setting
every start-of-step value equal to the current one. The state a generic integrator carries is the specific internal
energy and the time-dependent abundances (model.time_dependent), as the C solver's is; the other solve species are
steady-state variables and are solved algebraically at every evaluation, together with the temperature behind u,
by Newton's method with the network's analytic Jacobian, warm-started from the previous evaluation. A balance
variable (the dust temperature) is solved at every evaluation as the stable root of its row nearest the previous
value, as the C solver and GIZMO's dust temperature search do. Expressions are evaluated with the generated C code's
floating-point semantics (jaco.model_check), tables included, so a model that checks also integrates.

The parameters a host supplies (n_Htot, pdv_work, columns, velocity gradients, ...) are the caller's: see
jaco.onezone for a one-zone parcel that closes them.
"""

import numpy as np
import sympy as sp
from scipy.optimize import brentq

from .declarations import implied_parameters, Species
from .interpolation import tables_in
from .model_check import _c_lambdify, _parameter_values
from .symbols import sanitize_name, sanitize_symbols

__all__ = ["OdeSystem"]

BALANCE_STEP = 1.1  # bracketing walk of a balance variable: steps of 10%, growing by 10%
BALANCE_MAXITER = 200
NEWTON_MAXITER = 100
NEWTON_TOL = 1e-10  # on the relative change of every unknown
NEWTON_MAX_STEP = {"T": np.log(2.0), "x": np.log(10.0)}  # per iteration, in the logs of the unknowns
BALANCE_PASSES = 20  # Newton passes with the balance variables re-solved in between
PSEUDO_MAXSTEPS = 400  # pseudo-transient fallback: implicit steps, and Newton iterations within each
PSEUDO_INNER = 12


class OdeSystem:
    """The explicit rate equations of a model.

    Parameters
    ----------
    model: jaco.model.Model
    params: dict, optional
        Default parameter values by name (as in the expressions, e.g. "G_0", "∇v"), overriding the declarations'
        defaults; the values passed to rates() override these in turn. Parameters left without a value take the
        declared default, solar abundances for element totals, or 1 (reported in `missing`).

    Attributes
    ----------
    species: list of str
        Every solve species (e.g. ["H+", "He+", "He++", "H_2"]).
    integrated: list of str
        The time-dependent ones, in the state-vector order after u (e.g. ["H_2"]).
    steady: list of str
        The steady-state ones, solved at every evaluation (e.g. ["H+", "He+", "He++"]).
    balances: list of str
        The balance variables solved at every evaluation (e.g. ["Td"]).
    parameters: list of str
        Every parameter the right-hand side reads, by its expression name.
    """

    def __init__(self, model, params=None):
        self.model = model
        net = model.network
        solve_vars, td = list(model.solve_vars), list(model.time_dependent)
        rhs, jac, idx = net.solver_functions(solve_vars, td, return_jac=True)
        variables = list(idx)
        prelude = list(net._intermediate_prelude)  # the intermediates and their derivatives
        self._variables = variables
        names = [str(v) for v in variables]
        self.species = [n[2:] for n in names if n.startswith("x_")]
        self.integrated = [s for s in self.species if s in td]
        self.steady = [s for s in self.species if s not in td]
        self.balances = [n for n in names if n not in ("u", "T") and not n.startswith("x_")]
        decl_vars = getattr(net, "variables", {})
        self._balance_bounds = {b: (decl_vars[b].floor, decl_vars[b].ceiling) for b in self.balances}
        self._balance_row = {b: names.index(b) for b in self.balances}
        meta = net._solver_metadata(idx, td)
        self._floor = {}
        self._ceiling = {}
        for v, m in zip(variables, meta["variables"]):
            n = str(v)
            if n.startswith("x_"):
                self._floor[n[2:]], self._ceiling[n[2:]] = m["floor"], m["ceiling"]
        # the conservation budgets of the eliminated abundances, total - sum w x >= 0 (an H budget x_H+ + 2 x_H2 <= 1):
        # the state is projected into them before every evaluation, as the C solver keeps it
        by_c_name = {sanitize_name(f"x_{s}"): s for s in self.species}
        self._budgets = [(b["total"], b["total_param"], [(by_c_name[t], w) for t, w in b["terms"] if t in by_c_name])
                         for b in meta["budgets"]]

        # the rows with the backward-Euler terms removed: start-of-step values equal to the current ones
        u_sym, u_expr = sp.Symbol("u"), net._reduce_expression(net.eos.internal_energy)
        same = {sp.Symbol("u_initial"): u_sym}
        for s in self.species:
            same[sp.Symbol(f"x_{s}_initial")] = sp.Symbol(f"x_{s}")

        def reduce(e):
            return sp.sympify(e).xreplace(same).xreplace({u_sym: u_expr})

        rows = [reduce(e) for e in rhs]
        self._i_heat = names.index("T")
        self._i_species = [names.index(f"x_{s}") for s in self.species]
        P_expr = net._reduce_expression(net.eos.pressure)
        rho_expr = net._reduce_expression(net.eos.density)
        du_dT = sp.diff(u_expr, sp.Symbol("T"))
        du_dx = [sp.diff(u_expr, sp.Symbol(f"x_{s}")) for s in self.species]
        extras = [u_expr, P_expr, rho_expr, du_dT] + du_dx
        self._i_extras = len(rows)
        # every process's heat row, reduced like the equations, for heat_breakdown()
        self._process_names = []
        process_heat = []
        for proc in model.subprocesses:
            e = proc.network.get("heat")
            if e is None or sp.sympify(e.rhs) == 0:
                continue
            self._process_names.append(proc.name)
            process_heat.append(reduce(net._reduce_expression(sp.sympify(e.rhs))))
        self._i_process = len(rows) + len(extras)
        # the Jacobian of the algebraic rows (steady-state species, balance variables) with respect to T, the
        # steady-state abundances and the balance variables, with u's dependence on them carried through (the network
        # differentiates at fixed u, one of its variables)
        self._i_jac = self._i_process + len(process_heat)
        self._alg = [f"x_{s}" for s in self.steady] + list(self.balances)
        self._jac_cols = ["T"] + self._alg
        jac_block = []
        for rname in self._alg:
            i = names.index(rname)
            d_du = sp.diff(sp.sympify(rhs[i]), u_sym)
            for c in self._jac_cols:
                j = names.index(c)
                e = sp.sympify(jac[i][j])
                if d_du != 0 and c != "T" and c.startswith("x_"):
                    e = e + d_du * du_dx[self.species.index(c[2:])]
                elif d_du != 0 and c == "T":
                    e = e + d_du * du_dT
                jac_block.append(reduce(e))
        exprs = rows + extras + process_heat + jac_block

        syms = set().union(*[e.free_symbols for e in exprs + [W for _, W in prelude]])
        local = {w for w, _ in prelude}
        self.parameters = sorted({str(s) for s in syms} - set(names) - {str(w) for w in local} - {"Δt"})
        species_decl = {s: Species(s, k) for s, k in net.species_kinds.items()}
        declared = set(net.parameters) | {p.name for p in implied_parameters(species_decl, td)} | \
            {"n_Htot", "T", "pdv_work", "y", "Y", "Z", "C_2", "C_3"}
        self.undeclared = [p for p in self.parameters if p not in declared]
        self._defaults, self.missing = _parameter_values(net, self.parameters, td, dict(params or {}))

        tables = tables_in(exprs + [W for _, W in prelude])
        var_args = [sanitize_symbols(v) for v in variables]
        par_args = [sanitize_symbols(sp.Symbol(p)) for p in self.parameters]
        pre_args = [sanitize_symbols(w) for w, _ in prelude]
        self._pre = [_c_lambdify(var_args + par_args + pre_args[:k], [sanitize_symbols(W)], tables)
                     for k, (_, W) in enumerate(prelude)]
        self._f = _c_lambdify(var_args + par_args + pre_args, [sanitize_symbols(e) for e in exprs], tables)
        self._names = names
        self.reset()

    # ------------------------------------------------------------------ evaluation

    def _evaluate(self, T, x, balance, params):
        """every expression at one state: rows, then u, P, rho, du/dT, du/dx_s, the process heats, the Jacobian block"""
        vals = []
        for n in self._names:
            if n == "u":
                vals.append(0.0)  # u is substituted by u(T, x) in every row; the slot is only the argument order
            elif n == "T":
                vals.append(float(T))
            elif n.startswith("x_"):
                vals.append(float(x[n[2:]]))
            else:
                vals.append(float(balance[n]))
        pv = [float(params[p]) if p != "Δt" else 1e300 for p in self.parameters]
        pre = []
        for f in self._pre:
            pre.append(f(*(vals + pv + pre))[0])
        out = self._f(*(vals + pv + pre))
        return [float(o) for o in out]

    def _clip(self, x, params):
        """the state within its floors, ceilings and budgets (a full budget scaled to 1 - 1e-10 of its total)"""
        out = {}
        for s in self.species:
            out[s] = min(max(float(x[s]), self._floor[s]), self._ceiling[s])
        for total, total_param, terms in self._budgets:
            cap = (float(params[total_param]) if total_param else total) * (1 - 1e-10)
            used = sum(w * out[s] for s, w in terms)
            if used > cap > 0:
                for s, _ in terms:
                    out[s] = max(self._floor[s], out[s] * cap / used)
        return out

    def _solve_balance(self, name, T, x, balance, params):
        """the stable root of the balance row nearest the previous value: walk from it in the direction the row
        points, by steps of 10% growing by 10%, to the first sign change, then Brent"""
        lo, hi = self._balance_bounds[name]
        row = self._balance_row[name]

        def F(v):
            b = dict(balance)
            b[name] = v
            return self._evaluate(T, x, b, params)[row]

        a = min(max(balance[name], lo), hi)
        fa = F(a)
        if fa == 0:
            return a
        up = fa > 0
        fac = BALANCE_STEP
        for _ in range(BALANCE_MAXITER):
            if (up and a >= hi) or (not up and a <= lo):
                return a  # the row keeps its sign to the bound: the bound is the answer
            b = min(hi, a * fac) if up else max(lo, a / fac)
            fb = F(b)
            if fb == 0:
                return b
            if (fb > 0) != (fa > 0):
                return brentq(F, min(a, b), max(a, b), xtol=1e-12, rtol=1e-10)
            a, fa = b, fb
            fac *= BALANCE_STEP
        return a

    # The algebraic unknowns, in the logs: [ln T if u is prescribed] + ln x of the steady species + ln of the balance
    # variables. The residuals: (u - u_target)/u_target, the steady rows, the balance rows.

    def _system(self, T, u_target, x, balance, params):
        """the residuals and the Jacobian (d/d ln of the unknowns) of the algebraic system at a state"""
        solve_T = u_target is not None
        n_a = len(self._alg)
        k0 = 1 if solve_T else 0
        out = self._evaluate(T, x, balance, params)
        u = out[self._i_extras]
        du_dT, du_dx = out[self._i_extras + 3], out[self._i_extras + 4:self._i_extras + 4 + len(self.species)]
        vals = [x[s] for s in self.steady] + [balance[b] for b in self.balances]
        F, J = np.empty(k0 + n_a), np.zeros((k0 + n_a, k0 + n_a))
        if solve_T:
            F[0] = (u - u_target) / u_target
            J[0, 0] = du_dT * T / u_target
            for c, s in enumerate(self.steady):
                J[0, 1 + c] = du_dx[self.species.index(s)] * x[s] / u_target
        for r, name in enumerate(self._alg):
            F[k0 + r] = out[self._names.index(name)]
            jrow = out[self._i_jac + r * (1 + n_a):self._i_jac + (r + 1) * (1 + n_a)]
            if solve_T:
                J[k0 + r, 0] = jrow[0] * T
            for c in range(n_a):
                J[k0 + r, k0 + c] = jrow[1 + c] * vals[c]
        return F, J, u

    @staticmethod
    def _solve_linear(J, F):
        """J dz = -F with the rows equilibrated first: the rows span many decades (a rate row at 1e-38 next to the
        energy row at 1), and a least-squares fallback would otherwise drop the small ones; an identically zero row
        (a balance variable whose physics is switched off) gets a zero step"""
        scale = np.max(np.abs(J), axis=1)
        keep = scale > 0
        Js, Fs = J[keep] / scale[keep, None], F[keep] / scale[keep]
        dz = np.zeros(len(F))
        try:
            dz[keep] = np.linalg.solve(Js[:, keep], -Fs)
        except np.linalg.LinAlgError:
            dz[keep] = np.linalg.lstsq(Js[:, keep], -Fs, rcond=1e-14)[0]
        return dz

    def _apply(self, T, x, balance, dz, solve_T, params):
        """the unknowns moved by dz in their logs (damped) and clipped; returns T, x, balance, the changes made"""
        k0 = 0
        if solve_T:
            dz[0] = np.clip(dz[0], -NEWTON_MAX_STEP["T"], NEWTON_MAX_STEP["T"])
            T = T * np.exp(dz[0])
            k0 = 1
        n_s = len(self.steady)
        dz[k0:k0 + n_s] = np.clip(dz[k0:k0 + n_s], -NEWTON_MAX_STEP["x"], NEWTON_MAX_STEP["x"])
        dz[k0 + n_s:] = np.clip(dz[k0 + n_s:], -NEWTON_MAX_STEP["T"], NEWTON_MAX_STEP["T"])
        new = dict(x)
        for c, s in enumerate(self.steady):
            new[s] = x[s] * np.exp(dz[k0 + c])
        new = self._clip(new, params)
        nb = dict(balance)
        for c, b in enumerate(self.balances):
            lo, hi = self._balance_bounds[b]
            nb[b] = min(max(balance[b] * np.exp(dz[k0 + n_s + c]), lo), hi)
        change = ([abs(dz[0])] if solve_T else []) + \
            [abs(np.log(new[s] / x[s])) if x[s] > 0 and new[s] > 0 else abs(new[s] - x[s]) for s in self.steady] + \
            [abs(np.log(nb[b] / balance[b])) for b in self.balances]
        return T, new, nb, change

    def _at_a_bound_wrongly(self, x, F, J, dz, k0):
        """True if a steady species sits at a bound that is not a root: at its floor with its row still producing it,
        unless the coupled step wants it lower and its own rate is destruction-dominated there (a species fed only by
        another floored species, falling with it); or at its ceiling with its row destroying it. A balance variable
        at its bound is the answer, as in the bracketing walk."""
        for r, s in enumerate(self.steady):
            k = k0 + r
            if x[s] <= self._floor[s] * (1 + 1e-9) and F[k] > 0 and not (dz[k] < 0 and J[k, k] < 0):
                return True
            if x[s] >= self._ceiling[s] * (1 - 1e-9) and F[k] < 0:
                return True
        return False

    def _newton(self, T, u_target, x, balance, params):
        """all the algebraic unknowns at once: damped Newton in their logs with the analytic Jacobian; returns
        (T, x, balance) or None if it does not converge"""
        solve_T = u_target is not None
        if not self._alg and not solve_T:
            return T, x, balance
        k0 = 1 if solve_T else 0
        x, balance = dict(x), dict(balance)
        for _ in range(NEWTON_MAXITER):
            F, J, _ = self._system(T, u_target, x, balance, params)
            if not np.all(np.isfinite(F)) or not np.all(np.isfinite(J)):
                return None
            dz = self._solve_linear(J, F)
            wanted = dz.copy()
            T, x, balance, change = self._apply(T, x, balance, dz, solve_T, params)
            if max(change) < NEWTON_TOL:
                return None if self._at_a_bound_wrongly(x, F, J, wanted, k0) else (T, x, balance)
        return None

    def _pseudo_transient(self, T, u_target, x, balance, params):
        """the steady state reached by evolving the steady species in time from x (implicit Euler in their logs with
        a step growing tenfold on success), T and the balance variables kept on their algebraic rows: finds the stable
        branch time evolution would, where Newton from a cold start does not. Returns (T, x, balance) or None."""
        solve_T = u_target is not None
        k0 = 1 if solve_T else 0
        n_s = len(self.steady)
        if n_s == 0:
            return self._newton(T, u_target, x, balance, params)
        x, balance = dict(x), dict(balance)
        F, J, _ = self._system(T, u_target, x, balance, params)
        rates = np.array([abs(J[k0 + r, k0 + r]) / max(x[s], self._floor[s]) for r, s in enumerate(self.steady)])
        dt = 0.1 / max(rates.max(), 1e-300)
        for _ in range(PSEUDO_MAXSTEPS):
            z_k = np.array([np.log(max(x[s], 1e-300)) for s in self.steady])
            Tn, xn, bn = T, dict(x), dict(balance)
            ok = False
            for _ in range(PSEUDO_INNER):
                F, J, _ = self._system(Tn, u_target, xn, bn, params)
                if not np.all(np.isfinite(F)) or not np.all(np.isfinite(J)):
                    break
                xv = np.array([max(xn[s], 1e-300) for s in self.steady])
                sl = slice(k0, k0 + n_s)
                # G = z - z_k - dt F/x for the species rows, dG/dz = I - dt (J - diag F)/x; the other rows stay Newton's
                G = F.copy()
                G[sl] = np.log(xv) - z_k - dt * F[sl] / xv
                M = J.copy()
                M[sl, :] = -dt * J[sl, :] / xv[:, None]
                M[sl, sl] += np.eye(n_s) + dt * np.diag(F[sl] / xv)
                dz = self._solve_linear(M, G)
                wanted = dz.copy()
                Tn, xn, bn, change = self._apply(Tn, xn, bn, dz, solve_T, params)
                if max(change) < 1e-8:
                    ok = True
                    break
            if not ok:
                dt /= 10.0
                if dt * rates.max() < 1e-12:
                    return None
                continue
            moved = max(abs(np.log(max(xn[s], 1e-300)) - z_k[r]) for r, s in enumerate(self.steady))
            T, x, balance = Tn, xn, bn
            if moved < NEWTON_TOL:
                F, J, _ = self._system(T, u_target, x, balance, params)
                return None if self._at_a_bound_wrongly(x, F, J, wanted, k0) else (T, x, balance)
            dt *= 10.0
        return None

    def _solve_state(self, T, u_target, x, balance, params):
        """the algebraic unknowns at a state: the balance variables walked to a root when none is known, then Newton
        on everything from the warm start; if that fails, the pseudo-transient fallback with the balances re-walked
        between passes"""
        if self.balances and any(balance[b] is None for b in self.balances):
            for b in self.balances:
                balance[b] = min(max(T, self._balance_bounds[b][0]), self._balance_bounds[b][1])
            for b in self.balances:
                balance[b] = self._solve_balance(b, T, x, balance, params)
        res = self._newton(T, u_target, x, balance, params)
        if res is not None:
            return res
        for _ in range(BALANCE_PASSES):
            T, x, balance = self._newton_from_scratch(T, u_target, x, balance, params)
            moved = False
            for b in self.balances:
                old = balance[b]
                balance[b] = self._solve_balance(b, T, x, balance, params)
                moved |= abs(balance[b] - old) > 1e-8 * abs(old)
            if not moved:
                break
        return T, x, balance

    def _newton_from_scratch(self, T, u_target, x, balance, params):
        """when Newton from the warm start fails: the pseudo-transient from it, then from cold starts"""
        starts = [x, {**x, **{s: self._floor[s] for s in self.steady}}, {**x, **{s: 1e-3 for s in self.steady}}]
        T_starts = [T] + ([1e3, 1e4, 1e2, 1e6] if u_target is not None else [])
        for T0 in T_starts:
            for st in starts:
                res = self._pseudo_transient(T0, u_target, dict(st), dict(balance), params)
                if res is not None:
                    return res
        raise RuntimeError(f"OdeSystem: the steady-state solve did not converge at T = {T:.4g}, "
                           f"n_Htot = {params.get('n_Htot')}, x = {x}")

    def _prepare(self, T, x, params, solved):
        """the parameters, abundances and balance variables of one evaluation, the steady-state species solved unless
        the caller says they already are"""
        p = dict(self._defaults)
        p.update(params)
        xx = {s: x.get(s, None) for s in self.species}
        for s in self.species:
            if xx[s] is None:
                xx[s] = self._last_x.get(s) if s in self.steady and self._last_x.get(s) is not None else self._floor[s]
        xx = self._clip(xx, p)
        balance = dict(self._last_balance)
        if not solved:
            T, xx, balance = self._solve_state(T, None, xx, balance, p)
        elif self.balances and any(balance[b] is None for b in self.balances):  # solved elsewhere: the balances
            for b in self.balances:                                              # of the last solve stand
                balance[b] = min(max(T, self._balance_bounds[b][0]), self._balance_bounds[b][1])
            for b in self.balances:
                balance[b] = self._solve_balance(b, T, xx, balance, p)
        self._remember(T, xx, balance)
        return p, xx, balance

    def _remember(self, T, x, balance):
        self._last_T = T
        self._last_x = dict(x)
        self._last_balance = dict(balance)

    def rates(self, T, x, params, solved=False):
        """Everything at one state.

        Parameters
        ----------
        T: float
            Gas temperature [K].
        x: dict
            Abundance per H nucleus of each solve species; a missing time-dependent one is at its floor, a missing
            steady-state one starts from the previous evaluation's value.
        params: dict
            Parameter values by expression name for this evaluation (n_Htot, pdv_work, ...); the rest from the
            system's defaults.
        solved: bool
            True if x already holds the steady-state abundances at this state (from temperature()); otherwise they
            are solved here, the values in x the starting guess.

        Returns
        -------
        dict with heat [erg cm^-3 s^-1] (the gas's net heating including pdv_work), x (every abundance, the
        steady-state ones solved), dxdt {species: per H per s} (zero for a solved steady-state species), dudt
        [erg g^-1 s^-1], dTdt [K s^-1] at fixed steady-state abundances, u [erg g^-1], P [erg cm^-3], rho [g cm^-3],
        cv (du/dT) and the balance variables' values.
        """
        p, xx, balance = self._prepare(T, x, params, solved)
        out = self._evaluate(T, xx, balance, p)
        k = self._i_extras
        u, P, rho, cv = out[k], out[k + 1], out[k + 2], out[k + 3]
        du_dx = out[k + 4:k + 4 + len(self.species)]
        n = float(p["n_Htot"])
        dxdt = {s: out[i] / n for s, i in zip(self.species, self._i_species)}
        heat = out[self._i_heat]
        dTdt = (heat / rho - sum(d * dxdt[s] for d, s in zip(du_dx, self.species) if s in self.integrated)) / cv
        return dict(heat=heat, x=dict(xx), dxdt=dxdt, dudt=heat / rho, dTdt=dTdt, u=u, P=P, rho=rho, cv=cv, **balance)

    def thermo(self, T, x, params):
        """u, P, rho at a state as given (nothing solved; missing species as in rates())"""
        p = dict(self._defaults)
        p.update(params)
        xx = {s: x.get(s, None) for s in self.species}
        for s in self.species:
            if xx[s] is None:
                xx[s] = self._last_x.get(s) if s in self.steady and self._last_x.get(s) is not None else self._floor[s]
        xx = self._clip(xx, p)
        balance = {b: (self._last_balance[b] if self._last_balance[b] is not None else T) for b in self.balances}
        out = self._evaluate(T, xx, balance, p)
        k = self._i_extras
        return dict(u=out[k], P=out[k + 1], rho=out[k + 2], x=xx)

    def temperature(self, u, x, params, T_guess=None):
        """The gas temperature and every abundance behind a state (u, time-dependent abundances): T and the
        steady-state species solved together from u = u(T, x) and their rows. Returns (T, x)."""
        p = dict(self._defaults)
        p.update(params)
        xx = {s: x.get(s, None) for s in self.species}
        for s in self.species:
            if xx[s] is None:
                xx[s] = self._last_x.get(s) if s in self.steady and self._last_x.get(s) is not None else self._floor[s]
        xx = self._clip(xx, p)
        T = T_guess if T_guess is not None else (self._last_T if self._last_T is not None else 1e3)
        T, xx, balance = self._solve_state(T, float(u), xx, dict(self._last_balance), p)
        self._remember(T, xx, balance)
        return T, xx

    def heat_breakdown(self, T, x, params, solved=False):
        """process name -> its heat row [erg cm^-3 s^-1] at the state (the steady-state and balance variables solved as
        in rates()); positive heats the gas. The sum is rates()['heat']."""
        p, xx, balance = self._prepare(T, x, params, solved)
        out = self._evaluate(T, xx, balance, p)
        return dict(zip(self._process_names, out[self._i_process:self._i_jac]))

    def rhs(self, y, params):
        """dy/dt for y = (u, x_1, ..., x_m) in self.integrated order"""
        T, x = self.temperature(y[0], {s: v for s, v in zip(self.integrated, y[1:])}, params)
        r = self.rates(T, x, params, solved=True)
        return np.array([r["dudt"]] + [r["dxdt"][s] for s in self.integrated])

    def state(self, T, x, params):
        """the state vector (u, time-dependent abundances) for T and a dict of abundances (missing time-dependent
        species at their floor, steady-state ones solved)"""
        r = self.rates(T, x, params)
        return np.array([r["u"]] + [r["x"][s] for s in self.integrated], dtype=float)

    def reset(self):
        """forget the warm starts (balance variables, steady-state abundances, temperature)"""
        self._last_balance = {b: None for b in self.balances}
        self._last_x = {}
        self._last_T = None

    @property
    def floors(self):
        return dict(self._floor)

    @property
    def defaults(self):
        """the parameter values used when a call gives none"""
        return dict(self._defaults)
