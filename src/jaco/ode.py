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
from .model_check import _parameter_values
from .ode_codegen import emit_source, compile_source
from .symbols import sanitize_name, sanitize_symbols

__all__ = ["OdeSystem"]

# the parameters a host closes from the cell's state (jaco.onezone.Closures, GizmoHostClosures): the derivative
# block carries d/d(parameter) for these, so a Jacobian can chain through a closure
HOST_PARAMETERS = ("pdv_work", "rho", "∇v", "∇v_tf", "N_H", "N_H_cool", "Δx", "G_0", "A_dust", "kappa_thick", "lambda_J")
BALANCE_STEP = 1.1  # bracketing walk of a balance variable: steps of 10%, growing by 10%
BALANCE_MAXITER = 200
NEWTON_MAXITER = 100
NEWTON_TOL = 1e-10  # on the relative change of every unknown
NEWTON_MAX_STEP = {"T": np.log(2.0), "x": np.log(10.0)}  # per iteration, in the logs of the unknowns
BALANCE_PASSES = 20  # Newton passes with the balance variables re-solved in between
PSEUDO_MAXSTEPS = 400  # pseudo-transient fallback: implicit steps, and Newton iterations within each
PSEUDO_INNER = 12
STEP_TOL = 1e-6  # solve_step: relative change of every unknown at convergence (GIZMO's JACO_TOL_DEFAULT)
STEP_STALL = 10.0  # a full step within this many tol that the line search cannot improve on is round-off: converged
STEP_MAXITER = 100
STEP_BACKTRACK = 12  # halvings of the Newton step in the line search
STEP_T_BOUNDS = (1.0, 1e10)  # the gas temperature's bounds in a step [K], as GIZMO's JACO solver's
STEP_LINEAR_SCALE = 1e-30  # the smallest value a linear unknown's relative change is measured against


class OdeSystem:
    """The explicit rate equations of a model.

    Parameters
    ----------
    model: jaco.model.Model
    params: dict, optional
        Default parameter values by name (as in the expressions, e.g. "G_0", "∇v"), overriding the declarations'
        defaults; the values passed to rates() override these in turn. Parameters left without a value take the
        declared default, solar abundances for element totals, or 1 (reported in `missing`).
    step: bool
        Keep the backward-Euler terms: the rows are those of one implicit step over Δt from the start-of-step values
        (u_initial, x_<species>_initial), which become parameters with Δt, as the generated code's are; solve_step()
        solves them. For a model written as a step (starforge_legacy_RT: the kick's absorption factors exp(-a Δt), a
        rate on x_H+ - x_H+_initial, the IR band's temperature in and out), which has no explicit rate equations. The
        model's declared outputs are evaluated at the solution (outputs).

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

    def __init__(self, model, params=None, backend="numba", step=False):
        self.model = model
        self.backend = backend
        self.step_mode = step
        net = model.network
        solve_vars, td = list(model.solve_vars), list(model.time_dependent)
        rhs, _, idx = net.solver_functions(solve_vars, td, return_jac=True)  # return_jac builds the intermediates' derivatives
        variables = list(idx)
        prelude = list(net._intermediate_prelude)  # the intermediates and their derivatives d{w}_d{j}
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

        # the rows with the backward-Euler terms removed: start-of-step values equal to the current ones (kept in step
        # mode, where the start-of-step values are parameters)
        u_sym, u_expr = sp.Symbol("u"), net._reduce_expression(net.eos.internal_energy)
        same = {sp.Symbol("u_initial"): u_sym}
        for s in self.species:
            same[sp.Symbol(f"x_{s}_initial")] = sp.Symbol(f"x_{s}")
        if step:
            same = {}

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
        # the derivatives the solves need: of du/dt = heat/rho, the time-dependent rates row/n, the steady-state
        # rows and the balance rows, with respect to T, every abundance and the balance variables; the intermediates'
        # derivatives come from the prelude's symbols d{w}_d{j} (j the variable's index), zero where the network found
        # none
        for _, W in prelude:
            if u_sym in W.free_symbols:
                raise NotImplementedError("an intermediate depends on u; the explicit right-hand side substitutes u(T, x)")
        prelude_names = {str(w) for w, _ in prelude}
        intermediates = [w for w, _ in net._intermediates]
        # the host parameters a one-zone closure may tie to the state also get derivative columns, with the
        # intermediates' derivatives with respect to them added to the prelude
        used = set().union(*[e.free_symbols for e in rows + extras + process_heat + [W for _, W in prelude]])
        host = [sp.Symbol(q) for q in HOST_PARAMETERS if sp.Symbol(q) in used]
        Dp = {}
        for k, (w, W) in enumerate(net._intermediates):
            for q in host:
                d = sp.diff(W, q)
                for wj, _ in net._intermediates[:k]:
                    if (wj, q) in Dp and wj in W.free_symbols:
                        d = d + sp.diff(W, wj) * Dp[(wj, q)]
                if d != 0:
                    Dp[(w, q)] = sp.Symbol(f"d{w}_dp_{sanitize_name(str(q))}")
                    prelude.append((Dp[(w, q)], d))

        def total_derivative(expr, var):
            d = sp.diff(expr, var)
            is_param = str(var) not in names
            for w in intermediates:
                if w in expr.free_symbols:
                    D = Dp.get((w, var)) if is_param else sp.Symbol(f"d{w}_d{names.index(str(var))}")
                    if D is not None and str(D) in prelude_names | {str(v) for v in Dp.values()}:
                        d = d + sp.diff(expr, w) * D
            return d

        n_sym = sp.Symbol("n_Htot")
        self._deriv_rows = (["dudt"] + [f"x_{s}" for s in self.integrated] + [f"x_{s}" for s in self.steady]
                            + list(self.balances))
        self._deriv_params = [str(q) for q in host]
        self._deriv_cols = ["T"] + [f"x_{s}" for s in self.species] + list(self.balances) + self._deriv_params
        self._alg = [f"x_{s}" for s in self.steady] + list(self.balances)  # the algebraic unknowns besides T
        row_of = {n: rows[i] for i, n in enumerate(names)}
        F_exprs = ([rows[self._i_heat] / rho_expr] + [row_of[f"x_{s}"] / n_sym for s in self.integrated]
                   + [row_of[f"x_{s}"] for s in self.steady] + [row_of[b] for b in self.balances])
        self._i_deriv = self._i_process + len(process_heat)
        deriv_block = [total_derivative(F, sp.Symbol(c)) for F in F_exprs for c in self._deriv_cols]
        self._deriv = (len(F_exprs), len(self._deriv_cols))
        exprs = rows + extras + process_heat + deriv_block
        # step mode: the model's outputs (e.g. the IR band's radiation temperature after the step), at the solution
        self.output_names = [o[0] for o in net._outputs] if step else []
        self._i_outputs = len(exprs)
        if step:
            exprs = exprs + [reduce(o[1]) for o in net._outputs]

        syms = set().union(*[e.free_symbols for e in exprs + [W for _, W in prelude]])
        local = {w for w, _ in prelude}
        self.parameters = sorted({str(s) for s in syms} - set(names) - {str(w) for w in local}
                                 - (set() if step else {"Δt"}))
        species_decl = {s: Species(s, k) for s, k in net.species_kinds.items()}
        declared = set(net.parameters) | {p.name for p in implied_parameters(species_decl, td)} | \
            {"n_Htot", "T", "pdv_work", "y", "Y", "Z", "C_2", "C_3"}
        self.undeclared = [p for p in self.parameters if p not in declared]
        self._defaults, self.missing = _parameter_values(net, self.parameters, td, dict(params or {}))

        tables = tables_in(exprs + [W for _, W in prelude])
        var_args = [sanitize_name(n) for n in names]
        par_args = [sanitize_name(p) for p in self.parameters]
        self.source = emit_source([sanitize_symbols(e) for e in exprs], var_args, par_args,
                                  [(sanitize_symbols(w), sanitize_symbols(W)) for w, W in prelude], tables,
                                  jax=(backend == "jax"), numpy_math=(backend == "numpy"))
        self._f = compile_source(self.source, backend)
        self._names = names
        self._i_param = {p: i for i, p in enumerate(self.parameters)}
        self.reset()

    # ------------------------------------------------------------------ evaluation

    def _evaluate(self, T, x, balance, params):
        """every expression at one state: rows, then u, P, rho, du/dT, du/dx_s, the process heats, the derivatives"""
        X = np.empty(len(self._names))
        for i, n in enumerate(self._names):
            if n == "u":
                X[i] = 0.0  # u is substituted by u(T, x) in every row; the slot is only the argument order
            elif n == "T":
                X[i] = T
            elif n.startswith("x_"):
                X[i] = x[n[2:]]
            else:
                X[i] = balance[n]
        P = np.empty(len(self.parameters))
        for p, i in self._i_param.items():
            P[i] = params[p] if (p != "Δt" or self.step_mode) else 1e300
        return np.asarray(self._f(X, P), dtype=float)

    def _derivatives(self, out):
        """the derivative block of an evaluation as a (rows, cols) array, rows self._deriv_rows, cols self._deriv_cols"""
        return out[self._i_deriv:self._i_deriv + self._deriv[0] * self._deriv[1]].reshape(self._deriv)

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
        D = self._derivatives(out)
        rows = [self._deriv_rows.index(a) for a in self._alg]
        cols = [self._deriv_cols.index(a) for a in self._alg]
        F, J = np.empty(k0 + n_a), np.zeros((k0 + n_a, k0 + n_a))
        if solve_T:
            F[0] = (u - u_target) / u_target
            J[0, 0] = du_dT * T / u_target
            for c, s in enumerate(self.steady):
                J[0, 1 + c] = du_dx[self.species.index(s)] * x[s] / u_target
        for r, name in enumerate(self._alg):
            F[k0 + r] = out[self._names.index(name)]
            if solve_T:
                J[k0 + r, 0] = D[rows[r], 0] * T
            for c in range(n_a):
                J[k0 + r, k0 + c] = D[rows[r], cols[c]] * vals[c]
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

    def _at_a_bound_wrongly(self, x, F, J, dz, k0, blind=()):
        """True if a steady species sits at a bound that is not a root: at its floor with its row still producing it,
        unless the coupled step wants it lower and its own rate is destruction-dominated there (a species fed only by
        another floored species, falling with it); or at its ceiling with its row destroying it. A balance variable
        at its bound is the answer, as in the bracketing walk. blind: rows whose derivatives were non-finite at this
        state; they give no direction, so the clamped value stands."""
        for r, s in enumerate(self.steady):
            k = k0 + r
            if k in blind:
                continue
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
            if not np.all(np.isfinite(F)):
                return None
            blind = {k for k in range(len(F)) if not np.all(np.isfinite(J[k]))}
            J = np.where(np.isfinite(J), J, 0.0)  # a derivative that is 0 * inf (a dead dust row at Z_d = 0) carries no information
            dz = self._solve_linear(J, F)
            wanted = dz.copy()
            T, x, balance, change = self._apply(T, x, balance, dz, solve_T, params)
            if max(change) < NEWTON_TOL:
                return None if self._at_a_bound_wrongly(x, F, J, wanted, k0, blind) else (T, x, balance)
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
        J = np.where(np.isfinite(J), J, 0.0)
        rates = np.array([abs(J[k0 + r, k0 + r]) / max(x[s], self._floor[s]) for r, s in enumerate(self.steady)])
        dt = 0.1 / max(rates.max(), 1e-300)
        for _ in range(PSEUDO_MAXSTEPS):
            z_k = np.array([np.log(max(x[s], 1e-300)) for s in self.steady])
            Tn, xn, bn = T, dict(x), dict(balance)
            ok = False
            for _ in range(PSEUDO_INNER):
                F, J, _ = self._system(Tn, u_target, xn, bn, params)
                if not np.all(np.isfinite(F)):
                    break
                J = np.where(np.isfinite(J), J, 0.0)
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
                blind = {k for k in range(len(F)) if not np.all(np.isfinite(J[k]))}
                return None if self._at_a_bound_wrongly(x, F, np.where(np.isfinite(J), J, 0.0), wanted, k0, blind) \
                    else (T, x, balance)
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
        return dict(zip(self._process_names, out[self._i_process:self._i_deriv]))

    def jacobian(self, y, params, dparams=None):
        """d(dy/dt)/dy for y = (u, x_1, ..., x_m): the algebraic unknowns behind the state (T, the steady-state species,
        the balance variables) differentiated implicitly through their rows; one sitting at a floor, ceiling or bound
        is clamped there and does not move. dparams chains in host parameters a closure ties to the state (those in
        HOST_PARAMETERS): {parameter: {"T": d/dT, "x_<species>": d/dx, ...}}, the closure's own partial derivatives
        at fixed density. Without it the parameters are held fixed. Non-finite where a parameter derivative is
        singular at the state (a fractional power at zero) and the parameter moves."""
        T, x = self.temperature(y[0], {s: v for s, v in zip(self.integrated, y[1:])}, params)
        p = dict(self._defaults)
        p.update(params)
        balance = self._last_balance
        out = self._evaluate(T, x, balance, p)
        D = self._derivatives(out)
        du_dT, du_dx = out[self._i_extras + 3], out[self._i_extras + 4:self._i_extras + 4 + len(self.species)]
        td = [f"x_{s}" for s in self.integrated]
        free_steady = [s for s in self.steady
                       if self._floor[s] * (1 + 1e-9) < x[s] < self._ceiling[s] * (1 - 1e-9)]
        free_bal = [b for b in self.balances
                    if self._balance_bounds[b][0] * (1 + 1e-9) < balance[b] < self._balance_bounds[b][1] * (1 - 1e-9)]
        z = ["T"] + [f"x_{s}" for s in free_steady] + free_bal
        cz = [self._deriv_cols.index(c) for c in z]
        ct = [self._deriv_cols.index(c) for c in td]
        rF = [self._deriv_rows.index(r) for r in ["dudt"] + td]
        rG = [self._deriv_rows.index(r) for r in z[1:]]
        n_y, n_z = 1 + len(td), len(z)
        # G = (u(T, x) - u, the free steady rows, the free balance rows) = 0 defines z(y)
        G_z = np.zeros((n_z, n_z))
        G_y = np.zeros((n_z, n_y))
        G_z[0, 0] = du_dT
        for k, sname in enumerate(free_steady):
            G_z[0, 1 + k] = du_dx[self.species.index(sname)]
        G_y[0, 0] = -1.0
        for k, sname in enumerate(self.integrated):
            G_y[0, 1 + k] = du_dx[self.species.index(sname)]
        if rG:
            G_z[1:, :] = D[np.ix_(rG, cz)]
            G_y[1:, 1:] = D[np.ix_(rG, ct)]
        F_z = D[np.ix_(rF, cz)]
        F_y = np.zeros((n_y, n_y))
        F_y[:, 1:] = D[np.ix_(rF, ct)]
        if dparams:
            # the parameters' own derivatives along z and along y (y's u column is 0: u reaches them through T)
            for q, dq in dparams.items():
                if q not in self._deriv_params:
                    continue
                pz = np.array([dq.get(c, 0.0) for c in z], dtype=float)
                py = np.array([0.0] + [dq.get(c, 0.0) for c in td], dtype=float)
                if not (np.any(pz) or np.any(py)):
                    continue  # a parameter that does not move: its (possibly singular) derivative is irrelevant
                c = self._deriv_cols.index(q)
                if rG:
                    G_z[1:, :] += np.outer(D[rG, c], pz)
                    G_y[1:, :] += np.outer(D[rG, c], py)
                F_z += np.outer(D[rF, c], pz)
                F_y += np.outer(D[rF, c], py)
        # an unknown whose row is identically zero (dust physics switched off at Z_d = 0) cannot move either
        keep = np.max(np.abs(G_z), axis=1) > 0
        G_z, G_y, F_z = G_z[np.ix_(keep, keep)], G_y[keep], F_z[:, keep]
        scale = np.max(np.abs(G_z), axis=1)
        dz_dy = -np.linalg.solve(G_z / scale[:, None], G_y / scale[:, None])
        return F_y + F_z @ dz_dy

    # ------------------------------------------------------------------ step mode

    def _step_residual(self, T, x, balance, p, rows, cols, lin):
        """the step's residuals (the derivative block's rows: du/dt, the time-dependent rows / n, the steady rows, the
        balance rows) and their Jacobian in the step's unknowns (ln T, ln x or x for the linear ones, ln balances)"""
        out = self._evaluate(T, x, balance, p)
        n = float(p["n_Htot"])
        F = [out[self._i_heat] / out[self._i_extras + 2]]
        for s in self.species:
            v = out[self._names.index(f"x_{s}")]
            F.append(v / n if s in self.integrated else v)
        F += [out[self._balance_row[b]] for b in self.balances]
        F = np.array(F)
        J = self._derivatives(out)[np.ix_(rows, cols)]
        scale = np.array([T] + [1.0 if s in lin else x[s] for s in self.species] + [balance[b] for b in self.balances])
        return F, J * scale[None, :], out

    def _step_apply(self, T, x, balance, dz, lam, p, lin):
        """the unknowns moved by lam dz (the logs damped as in the state solves, the linear ones floored), clipped to
        their bounds and budgets"""
        n_s = len(self.species)
        dT = np.clip(lam * dz[0], -NEWTON_MAX_STEP["T"], NEWTON_MAX_STEP["T"])
        Tn = min(max(T * np.exp(dT), STEP_T_BOUNDS[0]), STEP_T_BOUNDS[1])
        xn = {}
        for c, s in enumerate(self.species):
            d = lam * dz[1 + c]
            xn[s] = x[s] + d if s in lin else x[s] * np.exp(np.clip(d, -NEWTON_MAX_STEP["x"], NEWTON_MAX_STEP["x"]))
        xn = self._clip(xn, p)
        bn = {}
        for c, b in enumerate(self.balances):
            lo, hi = self._balance_bounds[b]
            d = np.clip(lam * dz[1 + n_s + c], -NEWTON_MAX_STEP["T"], NEWTON_MAX_STEP["T"])
            bn[b] = min(max(balance[b] * np.exp(d), lo), hi)
        return Tn, xn, bn

    def _step_at_bounds(self, T, x, balance, lin):
        """the indices (in the step's unknowns: T, the species, the balance variables) at their lower and upper bounds"""
        lower, upper = set(), set()
        if T <= STEP_T_BOUNDS[0] * (1 + 1e-12):
            lower.add(0)
        if T >= STEP_T_BOUNDS[1] * (1 - 1e-12):
            upper.add(0)
        for c, s in enumerate(self.species):
            lo, hi = self._floor[s], self._ceiling[s]
            if (x[s] <= 0.0) if s in lin else (x[s] <= lo * (1 + 1e-9)):
                lower.add(1 + c)
            if x[s] >= hi * (1 - 1e-9):
                upper.add(1 + c)
        for c, b in enumerate(self.balances):
            lo, hi = self._balance_bounds[b]
            k = 1 + len(self.species) + c
            if balance[b] <= lo * (1 + 1e-12):
                lower.add(k)
            if balance[b] >= hi * (1 - 1e-12):
                upper.add(k)
        return lower, upper

    def _step_change(self, T, x, balance, Tn, xn, bn, lin):
        """the largest relative change of an unknown between two iterates"""
        change = [abs(np.log(Tn / T))]
        for s in self.species:
            if s in lin:
                change.append(abs(xn[s] - x[s]) / max(abs(xn[s]), abs(x[s]), STEP_LINEAR_SCALE))
            else:
                change.append(abs(np.log(xn[s] / x[s])))
        change += [abs(np.log(bn[b] / balance[b])) for b in self.balances]
        return max(change), int(np.argmax(change))

    def solve_step(self, T, x, balance, params, tol=STEP_TOL, maxiter=STEP_MAXITER, verbose=False):
        """One backward-Euler step (step mode): the gas temperature, every abundance and the balance variables at the
        end of a step of params["Δt"] from the start-of-step values in params (u_initial, x_<species>_initial, and the
        model's own, e.g. Td_initial), every row solved at once by damped Newton with the analytic Jacobian, a
        backtracking line search on the row-equilibrated residual, from the guess (T, x, balance).

        The unknowns move in their logs (T, the abundances, the balance variables), but for the species with a zero
        floor (the radiation bands), which can be exactly 0 and whose rows are linear in them, which move linearly.
        Converged when no unknown changes by more than tol (relative). A balance variable whose row is identically zero
        (the dust temperature with no dust) keeps its guess.

        Returns (T, x, balance, outputs, iterations) or None if Newton does not converge; outputs are the model's
        declared outputs at the solution (output_names)."""
        if not self.step_mode:
            raise RuntimeError("solve_step needs OdeSystem(..., step=True)")
        p = dict(self._defaults)
        p.update(params)
        lin = {s for s in self.species if self._floor[s] <= 0}
        x = self._clip({s: max(float(x[s]), 0.0) if s in lin else max(float(x[s]), self._floor[s])
                        for s in self.species}, p)
        balance = {b: min(max(float(balance[b]), self._balance_bounds[b][0]), self._balance_bounds[b][1])
                   for b in self.balances}
        T = min(max(float(T), STEP_T_BOUNDS[0]), STEP_T_BOUNDS[1])
        rows = [self._deriv_rows.index(r) for r in ["dudt"] + [f"x_{s}" for s in self.species] + self.balances]
        cols = [self._deriv_cols.index(c) for c in ["T"] + [f"x_{s}" for s in self.species] + self.balances]
        names = ["T"] + list(self.species) + list(self.balances)
        stalled = False
        F, J, out = self._step_residual(T, x, balance, p, rows, cols, lin)
        for it in range(1, maxiter + 1):
            if not np.all(np.isfinite(F)):
                if verbose:
                    print(f"  step it {it}: non-finite residual {[n for n, f in zip(names, F) if not np.isfinite(f)]}")
                return None
            J = np.where(np.isfinite(J), J, 0.0)
            row_scale = np.max(np.abs(J), axis=1)
            row_scale[row_scale == 0] = 1.0
            # an unknown at a bound whose Newton step points out of it is pinned there: its row and column leave the
            # linear system and the merit (GIZMO's JACO solver, tier 1), so that a residual nothing can remove (a
            # species at its floor still being destroyed) neither steers the step nor stalls the line search
            lower, upper = self._step_at_bounds(T, x, balance, lin)
            # a balance variable nothing depends on (its column is zero: the dust temperature with no dust) leaves
            # Newton's system, which it would make singular; it is walked along its own row at the end
            free = {1 + len(self.species) + c for c in range(len(self.balances))
                    if not np.any(J[:, 1 + len(self.species) + c])}
            pinned = set(free)
            for _ in range(len(names)):
                active = [k for k in range(len(names)) if k not in pinned]
                dz = np.zeros(len(names))
                dz[active] = self._solve_linear(J[np.ix_(active, active)], F[active])
                # (a species row is its net production: one at its floor that its own row would grow is not pinned)
                out_of_bounds = {k for k in active
                                 if (k in lower and dz[k] < 0 and not (0 < k <= len(self.species) and F[k] > 0))
                                 or (k in upper and dz[k] > 0 and not (0 < k <= len(self.species) and F[k] < 0))}
                if not out_of_bounds:
                    break
                pinned |= out_of_bounds
            keep = np.array([k not in pinned for k in range(len(names))])
            merit = np.linalg.norm((F / row_scale)[keep])
            # converged when the full Newton step (clipped to the bounds) moves no unknown by more than tol, or by
            # little more where the last line search found nothing better than round-off (stalled)
            Tn, xn, bn = self._step_apply(T, x, balance, dz, 1.0, p, lin)
            full_change, worst = self._step_change(T, x, balance, Tn, xn, bn, lin)
            if full_change < tol or (stalled and full_change < STEP_STALL * tol):
                for c, b in enumerate(self.balances):  # the walk of OdeSystem's state solves: to the row's root nearest
                    if 1 + len(self.species) + c in free:  # the guess, or to the bound it keeps its sign to
                        bn[b] = self._solve_balance(b, Tn, xn, bn, p)
                Fn, Jn, out = self._step_residual(Tn, xn, bn, p, rows, cols, lin)
                if np.all(np.isfinite(Fn)):
                    outputs = dict(zip(self.output_names,
                                       out[self._i_outputs:self._i_outputs + len(self.output_names)]))
                    self._remember(Tn, xn, bn)  # the answer is the state later evaluations (rates, heats) start from
                    return Tn, xn, bn, outputs, it
            lam = 1.0
            stalled = False
            for _ in range(STEP_BACKTRACK):
                Fn, Jn, outn = self._step_residual(Tn, xn, bn, p, rows, cols, lin)
                if np.all(np.isfinite(Fn)) and np.linalg.norm((Fn / row_scale)[keep]) <= (1 - 1e-4 * lam) * merit:
                    break
                lam *= 0.5
                Tn, xn, bn = self._step_apply(T, x, balance, dz, lam, p, lin)
            else:  # no decrease down to the smallest step: take that one
                Fn, Jn, outn = self._step_residual(Tn, xn, bn, p, rows, cols, lin)
                stalled = True
            if verbose:
                print(f"  step it {it}: merit {merit:.3e} -> {np.linalg.norm((Fn / row_scale)[keep]):.3e}, "
                      f"lam {lam:.3g}, full-step change {full_change:.2e} (largest: {names[worst]}), "
                      f"T {T:.6g} -> {Tn:.6g}, pinned {[names[k] for k in sorted(pinned)]}")
            if not np.all(np.isfinite(Fn)):
                return None
            T, x, balance, F, J, out = Tn, xn, bn, Fn, Jn, outn
        return None

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
