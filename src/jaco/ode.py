"""A model as an explicit ODE right-hand side for a generic integrator (scipy.integrate.solve_ivp or any other).

    sys = OdeSystem(make_model())
    rates = sys.rates(T, x, params)             # heat [erg cm^-3 s^-1], dx/dt [per H per s], Td, u, P, rho, ...
    dydt = sys.rhs(y, params)                   # y = (T, x_1, ..., x_n) in sys.species order; dT/dt, dx/dt

The system is the model's reduced network as the code generator sees it: the steady-state closures, the conservation
reductions, the derived parameters and the intermediates are applied; the backward-Euler terms are removed by setting
every start-of-step value equal to the current one. The gas temperature follows from the energy equation through the
EOS, dT/dt = (heat / rho - sum_s du/dx_s dx_s/dt) / (du/dT), and a balance variable (the dust temperature) is solved
at every evaluation as the stable root of its row nearest the previous value, as the C solver and GIZMO's dust
temperature search do. Expressions are evaluated with the generated C code's floating-point semantics
(jaco.model_check), tables included, so a model that checks also integrates.

The parameters a host supplies (n_Htot, pdv_work, columns, velocity gradients, ...) are the caller's: see
jaco.freefall for a one-zone collapse that closes them.
"""

import numpy as np
import sympy as sp
from scipy.optimize import brentq

from .declarations import implied_parameters, Species
from .interpolation import tables_in
from .model_check import _c_lambdify, _parameter_values
from .symbols import sanitize_symbols

__all__ = ["OdeSystem"]

BALANCE_STEP = 1.1  # bracketing walk of a balance variable: steps of 10%, growing by 10%
BALANCE_MAXITER = 200


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
        The solve species, in the state-vector order after T (e.g. ["H+", "He+", "He++", "H_2"]).
    balances: list of str
        The balance variables solved at every evaluation (e.g. ["Td"]).
    parameters: list of str
        Every parameter the right-hand side reads, by its expression name.
    """

    def __init__(self, model, params=None):
        self.model = model
        net = model.network
        solve_vars, td = list(model.solve_vars), list(model.time_dependent)
        rhs, idx = net.solver_functions(solve_vars, td, return_jac=False)
        variables = list(idx)
        prelude = list(net._intermediates)
        self._variables = variables
        names = [str(v) for v in variables]
        self.species = [n[2:] for n in names if n.startswith("x_")]
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

        # the rows with the backward-Euler terms removed: start-of-step values equal to the current ones
        u_sym, u_expr = sp.Symbol("u"), net._reduce_expression(net.eos.internal_energy)
        same = {sp.Symbol("u_initial"): u_sym}
        for s in self.species:
            same[sp.Symbol(f"x_{s}_initial")] = sp.Symbol(f"x_{s}")
        rows = [sp.sympify(e).xreplace(same).xreplace({u_sym: u_expr}) for e in rhs]
        self._i_heat = names.index("T")
        self._i_species = [names.index(f"x_{s}") for s in self.species]
        P_expr = net._reduce_expression(net.eos.pressure)
        rho_expr = net._reduce_expression(net.eos.density)
        du_dT = sp.diff(u_expr, sp.Symbol("T"))
        du_dx = [sp.diff(u_expr, sp.Symbol(f"x_{s}")) for s in self.species]
        extras = [u_expr, P_expr, rho_expr, du_dT] + du_dx
        self._i_extras = len(rows)
        exprs = rows + extras

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
        self._last_balance = {b: None for b in self.balances}

    # ------------------------------------------------------------------ evaluation

    def _evaluate(self, T, x, balance, params):
        """every expression at one state: rows, then u, P, rho, du/dT, du/dx_s"""
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

    def _clip(self, x):
        out = {}
        for s in self.species:
            out[s] = min(max(float(x[s]), self._floor[s]), self._ceiling[s])
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

    def rates(self, T, x, params):
        """Everything at one state.

        Parameters
        ----------
        T: float
            Gas temperature [K].
        x: dict
            Abundance per H nucleus of each solve species (missing ones at their floor).
        params: dict
            Parameter values by expression name for this evaluation (n_Htot, pdv_work, ...); the rest from the
            system's defaults.

        Returns
        -------
        dict with heat [erg cm^-3 s^-1] (the gas's net heating including pdv_work), dxdt {species: per H per s},
        dTdt [K s^-1], u [erg g^-1], P [erg cm^-3], rho [g cm^-3], cv (du/dT) and the balance variables' values.
        """
        p = dict(self._defaults)
        p.update(params)
        xx = self._clip({s: x.get(s, self._floor[s]) for s in self.species})
        balance = {}
        for b in self.balances:
            guess = self._last_balance[b]
            balance[b] = guess if guess is not None else min(max(T, self._balance_bounds[b][0]), self._balance_bounds[b][1])
        for b in self.balances:
            balance[b] = self._solve_balance(b, T, xx, balance, p)
            self._last_balance[b] = balance[b]
        out = self._evaluate(T, xx, balance, p)
        k = self._i_extras
        u, P, rho, cv = out[k], out[k + 1], out[k + 2], out[k + 3]
        du_dx = out[k + 4:k + 4 + len(self.species)]
        n = float(p["n_Htot"])
        dxdt = {s: out[i] / n for s, i in zip(self.species, self._i_species)}
        heat = out[self._i_heat]
        dTdt = (heat / rho - sum(d * dxdt[s] for d, s in zip(du_dx, self.species))) / cv
        return dict(heat=heat, dxdt=dxdt, dTdt=dTdt, u=u, P=P, rho=rho, cv=cv, **balance)

    def rhs(self, y, params):
        """dy/dt for y = (T, x_1, ..., x_n) in self.species order"""
        T, x = y[0], {s: v for s, v in zip(self.species, y[1:])}
        r = self.rates(T, x, params)
        return np.array([r["dTdt"]] + [r["dxdt"][s] for s in self.species])

    def state(self, T, x):
        """the state vector for T and a dict of abundances (missing species at their floor)"""
        return np.array([T] + [x.get(s, self._floor[s]) for s in self.species], dtype=float)

    def reset(self):
        """forget the balance variables' warm starts"""
        self._last_balance = {b: None for b in self.balances}

    @property
    def floors(self):
        return dict(self._floor)
