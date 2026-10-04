"""check(model): evaluate every rate of a model, its partials in the solve variables and the model's outputs over a
(T, n_Htot, composition) grid with the floating-point semantics of the generated C code, compare the symbolic Jacobian
with finite differences, and flag undeclared symbols and non-finite values, so that a 0 * inf is found before the host
code meets it.

    report = jaco.check(make_model())        # the standard grid; jaco.check(model, "quick") for a coarse one
    print(report.summary()); assert report.ok

"C semantics": exp arguments clamped to [-500, 500], Heaviside(x) -> (x > 0), Max/Min -> fmax/fmin (which ignore a
NaN argument), IEEE double arithmetic (underflow to 0, pow(0, negative) = inf, 0 * inf = nan), the 1D, 2D and 3D table
interpolations of the generated helpers, and floats rounded as the C printer prints them.
"""

from dataclasses import dataclass, field

import numpy as np
import sympy as sp
from sympy.printing.numpy import NumPyPrinter

from .codegen.printers import _InterpPrinterMixin
from .declarations import CORE_PARAMETERS, implied_parameters, merge_declarations, Species
from .interpolation import tables_in
from .symbols import sanitize_symbols, sanitize_name

__all__ = ["check", "CheckGrid", "CheckReport", "Finding"]

FD_STEP = 1e-6  # finite-difference step, relative to the variable
FD_ATOL = 1e-7  # absolute tolerance of a partial: this fraction of its row's gross magnitude per unit of the variable


# --------------------------------------------------------------------------- #
# C semantics in numpy
# --------------------------------------------------------------------------- #

class _CSemanticsPrinter(NumPyPrinter):
    """numpy code with the semantics of JacoCCodePrinter's output"""

    _print_Float = _InterpPrinterMixin._print_Float

    def _print_exp(self, expr):
        return f"numpy.exp(numpy.fmax(-500.0, numpy.fmin(500.0, {self._print(expr.args[0])})))"

    def _print_Heaviside(self, expr):
        return f"(1.0*numpy.greater({self._print(expr.args[0])}, 0))"

    def _nested(self, func, args):
        out = self._print(args[0])
        for a in args[1:]:
            out = f"numpy.{func}({out}, {self._print(a)})"
        return out

    def _print_Max(self, expr):
        return self._nested("fmax", expr.args)

    def _print_Min(self, expr):
        return self._nested("fmin", expr.args)

    def _print_PiecewiseLinearInterp(self, expr):
        X, Y = tuple(float(v) for v in expr.args[1]), tuple(float(v) for v in expr.args[2])
        return f"_jaco_interp1d({self._print(expr.args[0])}, {X}, {Y}, {int(expr.args[3])})"

    def _print_PiecewiseConstantInterp(self, expr):
        X, V = tuple(float(v) for v in expr.args[1]), tuple(float(v) for v in expr.args[2])
        return f"_jaco_interp1d_const({self._print(expr.args[0])}, {X}, {V})"

    def _table(self, expr, nd, grad):
        args = ", ".join(self._print(a) for a in expr.args[:nd])
        return f"_jaco_table(({args},), '{expr.args[nd]}', {grad})"

    def _print_TableInterp2D(self, expr):
        return self._table(expr, 2, -1)

    def _print_TableInterp2D_dx(self, expr):
        return self._table(expr, 2, 0)

    def _print_TableInterp2D_dy(self, expr):
        return self._table(expr, 2, 1)

    def _print_TableInterp3D(self, expr):
        return self._table(expr, 3, -1)

    def _print_TableInterp3D_dx(self, expr):
        return self._table(expr, 3, 0)

    def _print_TableInterp3D_dy(self, expr):
        return self._table(expr, 3, 1)

    def _print_TableInterp3D_dz(self, expr):
        return self._table(expr, 3, 2)


def _interp1d(x, X, Y, extrap):
    """jaco_interp1d"""
    X, Y, x = np.asarray(X), np.asarray(Y), np.asarray(x, dtype=float)
    out = np.interp(x, X, Y)
    if extrap:
        lo, hi = x <= X[0], x >= X[-1]
        out = np.where(lo, Y[0] + (Y[1] - Y[0]) / (X[1] - X[0]) * (x - X[0]), out)
        out = np.where(hi, Y[-1] + (Y[-1] - Y[-2]) / (X[-1] - X[-2]) * (x - X[-1]), out)
    return out


def _interp1d_const(x, X, V):
    """jaco_interp1d_const: the value on the interval holding x, the end values beyond the ends"""
    X, V, x = np.asarray(X), np.asarray(V), np.asarray(x, dtype=float)
    i = np.clip(np.searchsorted(X, x, side="right") - 1, 0, len(V) - 1)
    return V[i]


def _table_function(tables):
    """jaco_table2d_eval / jaco_table3d_eval (jaco_table.h), vectorized; grad = -1 for the value, else the axis of the
    derivative. A query outside the table is clamped to the edge, with a zero derivative along that axis."""
    def table(xs, name, grad):
        t = tables[name]
        nd = t["ndim"]
        shape = np.broadcast(*[np.asarray(x, dtype=float) for x in xs]).shape
        f, cl, d, idx = [], [], [], []
        for a in range(nd):
            x = np.broadcast_to(np.asarray(xs[a], dtype=float), shape)
            log = bool(t["log_axes"][a])
            v = np.log(x) if log else x
            lo = np.log(t[f"axis{a}_min"]) if log else t[f"axis{a}_min"]
            hi = np.log(t[f"axis{a}_max"]) if log else t[f"axis{a}_max"]
            n = t["shape"][a]
            da = (hi - lo) / (n - 1)
            fa = (v - lo) / da
            clamped = (fa < 0) | (fa > n - 1)
            fa = np.clip(fa, 0, n - 1)
            ia = np.minimum(fa.astype(int), n - 2)
            f.append(fa - ia)
            idx.append(ia)
            cl.append(clamped)
            d.append(da * (x if log else 1.0))
        data = t["data"]
        val = np.zeros(shape)
        for corner in np.ndindex(*(2,) * nd):
            w = np.ones(shape)
            for a in range(nd):
                if a == grad:
                    w = w * (1.0 if corner[a] else -1.0)
                else:
                    w = w * (f[a] if corner[a] else 1 - f[a])
            val = val + w * data[tuple(idx[a] + corner[a] for a in range(nd))]
        if grad < 0:
            return val
        return np.where(cl[grad], 0.0, val / d[grad])
    return table


def _c_lambdify(args, exprs, tables):
    """Vectorized numpy function of exprs (a list) with C semantics"""
    modules = [{"_jaco_interp1d": _interp1d, "_jaco_interp1d_const": _interp1d_const,
                "_jaco_table": _table_function(tables)}, "numpy"]
    f = sp.lambdify(args, exprs, modules=modules, printer=_CSemanticsPrinter, cse=True)

    def wrapped(*vals):
        with np.errstate(all="ignore"):
            out = f(*vals)
        shape = np.broadcast_shapes(*[np.shape(v) for v in vals]) if vals else ()
        return [np.broadcast_to(np.asarray(o, dtype=float), shape) for o in out]
    return wrapped


# --------------------------------------------------------------------------- #
# grid
# --------------------------------------------------------------------------- #

class CheckGrid:
    """Where check() evaluates a model: every combination of a temperature, a density and a composition.

    The compositions are built from the model's species bounds and budgets, with every budget left positive by a
    relative margin of 1e-10 as the solver keeps it: all species at their floors; all at 1e-6 of their caps; each
    species alone at its cap (fully ionized, fully molecular, ...); all at 30% of their caps, scaled into the budgets;
    and n_random log-uniform draws scaled into the budgets.

    Parameters
    ----------
    T, n_Htot: array_like
        Temperatures [K] and H densities [cm^-3].
    n_random: int
        Random compositions besides the systematic ones.
    dt: float
        Time step [s] of the backward-Euler terms.
    seed: int
        Of the random compositions.
    """

    def __init__(self, T, n_Htot, n_random=0, dt=1e10, seed=0):
        self.T = np.asarray(T, dtype=float)
        self.n_Htot = np.asarray(n_Htot, dtype=float)
        self.n_random, self.dt, self.seed = int(n_random), float(dt), seed

    @classmethod
    def standard(cls):
        return cls(np.logspace(np.log10(3.0), 9, 25), np.logspace(-3, 9, 13), n_random=20)

    @classmethod
    def quick(cls):
        return cls(np.logspace(np.log10(3.0), 9, 7), np.logspace(-3, 9, 4), n_random=0)


# --------------------------------------------------------------------------- #
# report
# --------------------------------------------------------------------------- #

@dataclass
class Finding:
    """One problem: what (an expression, e.g. "process 'X' row H+", "system row x_H_2" or "output L"), which quantity
    ("value", "d/dT", ...), in how many of the evaluated points, for a derivative the worst error in units of the
    tolerance, and one example point."""

    what: str
    quantity: str
    count: int
    total: int
    worst: float = float("nan")
    example: dict = field(default_factory=dict)

    def __str__(self):
        ex = ", ".join(f"{k}={v:.4g}" for k, v in self.example.items())
        worst = "" if np.isnan(self.worst) else f", worst {self.worst:.3g} x tolerance"
        return f"{self.what}: {self.quantity} at {self.count}/{self.total} points{worst} (e.g. {ex})"


@dataclass
class CheckReport:
    """The result of check(): ok if no undeclared symbols, no non-finite values or partials, and the symbolic
    Jacobian agrees with finite differences everywhere."""

    name: str
    n_points: int
    n_expressions: int
    undeclared: list = field(default_factory=list)
    missing_values: list = field(default_factory=list)
    nonfinite: list = field(default_factory=list)
    jacobian: list = field(default_factory=list)

    @property
    def ok(self):
        return not (self.undeclared or self.missing_values or self.nonfinite or self.jacobian)

    def summary(self, limit=20):
        lines = [f"jaco.check {self.name}: {'OK' if self.ok else 'PROBLEMS'} ({self.n_expressions} expressions at "
                 f"{self.n_points} points)"]
        if self.undeclared:
            lines.append(f"  undeclared symbols: {', '.join(self.undeclared)}")
        if self.missing_values:
            lines.append(f"  parameters without a default (evaluated at 1): {', '.join(self.missing_values)}")
        for title, items in (("non-finite", self.nonfinite), ("Jacobian vs finite differences", self.jacobian)):
            if items:
                lines.append(f"  {title}: {len(items)}")
                lines += [f"    {f}" for f in items[:limit]]
                if len(items) > limit:
                    lines.append(f"    ... {len(items) - limit} more")
        return "\n".join(lines)

    __str__ = summary


# --------------------------------------------------------------------------- #
# check
# --------------------------------------------------------------------------- #

def _parameter_values(net, names, time_dependent, given):
    """name -> value of every parameter in names, from given, the declarations' defaults and solar abundances for the
    element totals; the names left without one"""
    from .data import SolarAbundances
    decls = merge_declarations("parameter", CORE_PARAMETERS, net.parameters.values())
    species = {s: Species(s, k) for s, k in net.species_kinds.items()}
    implied = {p.name: p for p in implied_parameters(species, time_dependent)}
    values, missing = {}, []
    for name in names:
        if name in given:
            values[name] = float(given[name])
            continue
        d = decls.get(name) or implied.get(name)
        if d is not None and d.default is not None:
            values[name] = float(d.default)
            continue
        el = name[2:].removesuffix(",tot") if name.startswith("x_") else None
        value = None
        if el and net.species_kinds.get(el) == "radiation":
            value = 1e-6  # photons per H nucleus: small, but not zero, so that the rates they drive are exercised
        elif el == "H":
            value = 1.0
        elif el:
            try:
                value = float(SolarAbundances.x(el))
            except (KeyError, ValueError, TypeError, AttributeError, NotImplementedError):
                value = None
        if value is None or not np.isfinite(value):
            value = 1.0
            missing.append(name)
        values[name] = value
    return values, missing


def _compositions(meta, species_vars, params, grid):
    """[(label, {var: abundance})]: see CheckGrid"""
    floor = {v: m["floor"] for v, m in species_vars.items()}
    ceil = {v: m["ceiling"] for v, m in species_vars.items()}
    budgets = [(b["total"] if b["total_param"] is None else params[b["total_param"]],
                [(t, w) for t, w in b["terms"]]) for b in meta["budgets"]]
    margin = 1 - 1e-10

    def cap(v):
        c = min(ceil[v], 1.0)
        for total, terms in budgets:
            for t, w in terms:
                if t == v:
                    c = min(c, total * margin / w)
        return c

    def fit(x):
        for total, terms in budgets:
            s = sum(w * x[t] for t, w in terms)
            if s > total * margin:
                for t, _ in terms:
                    x[t] = max(floor[t], x[t] * total * margin / s)
        return x

    comps = [("floor", dict(floor)), ("trace", fit({v: max(floor[v], 1e-6 * cap(v)) for v in floor}))]
    for v in species_vars:
        comps.append((f"{v} at cap", fit({u: (cap(u) if u == v else floor[u]) for u in floor})))
    comps.append(("30% of caps", fit({v: max(floor[v], 0.3 * cap(v)) for v in floor})))
    rng = np.random.default_rng(grid.seed)
    for i in range(grid.n_random):
        x = {}
        for v in floor:
            lo, hi = np.log(max(floor[v], 1e-30)), np.log(cap(v))
            x[v] = float(np.exp(rng.uniform(lo, hi)))
        comps.append((f"random {i}", fit(x)))
    return comps


def check(model, grid="standard", params=None, rtol=1e-3, name=None):
    """Evaluate a model's rates, their partials, its outputs and its reduced system with the generated C code's
    floating-point semantics over grid, and compare the symbolic Jacobian with finite differences.

    Parameters
    ----------
    model: jaco.model.Model
    grid: "standard", "quick" or CheckGrid
    params: dict, optional
        Parameter values (by their names in the expressions, e.g. "G_0"), overriding the declared defaults. The
        start-of-step values follow the state, u_initial = u(T, x), so that the backward-Euler terms vanish.
    rtol: float
        Relative tolerance of a partial against central differences, plus FD_ATOL of the gross magnitude of its row (the
        sum of the magnitudes of the processes' contributions and the backward-Euler term) per unit of the variable,
        times the round-off amplification of the eliminated abundances near a budget's end (total / remainder); a
        partial between the one-sided differences, as at a kink, passes, and so does one that agrees at either of
        two steps (1e-6 and 1e-5 of the variable and of what its budgets leave).
    name: str, optional
        For the report.

    Returns
    -------
    CheckReport
    """
    if isinstance(grid, str):
        grid = {"standard": CheckGrid.standard, "quick": CheckGrid.quick}[grid]()
    params = dict(params or {})
    name = name or repr(model)
    report = CheckReport(name, 0, 0)
    if not model.species:
        report.undeclared.append("(the model declares no species)")
        return report

    net = model.network
    solve_vars, td = list(model.solve_vars), list(model.time_dependent)
    rhs, jac, idx = net.solver_functions(solve_vars, td, return_jac=True)
    variables = list(idx)
    meta = net._solver_metadata(idx, td)
    prelude = list(net._intermediate_prelude)
    local = {w for w, _ in prelude}
    j_of = {v: j for j, v in enumerate(variables)}
    u_sym = sp.Symbol("u")

    # every rate: the rows of the processes as assembled, reduced like the equations
    rows = []
    for p in model.subprocesses:
        for key, eq in p.network.items():
            e = sp.sympify(eq.rhs)
            if e != 0:
                rows.append((f"process '{p.name}' row {key}", net._reduce_expression(e)))
    outputs = [(f"output {n}", e) for n, e, _, _ in net._outputs]
    system = [(f"system row {sanitize_name(str(v))}", r) for v, r in zip(variables, rhs)]

    def chain(e, j):
        """d e / d variable j, through the intermediates"""
        v = variables[j]
        d = sp.diff(e, v)
        for w in local:
            if w in e.free_symbols:
                D = sp.Symbol(f"d{w}_d{j}")
                if D in local:
                    d += sp.diff(e, w) * D
        return d

    free = [j for j, v in enumerate(variables) if v != u_sym]
    partials = [(what, j, chain(e, j)) for what, e in rows for j in free]

    exprs = [e for _, e in rows + outputs + system] + [d for _, _, d in partials] + [x for r in jac for x in r]
    syms = set().union(*[sp.sympify(e).free_symbols for e in exprs + [W for _, W in prelude]])
    param_names = sorted({str(s) for s in syms} - {str(v) for v in variables} - {str(w) for w in local})
    declared = set(merge_declarations("parameter", CORE_PARAMETERS, net.parameters.values())) | \
        {p.name for p in implied_parameters({s: Species(s, k) for s, k in net.species_kinds.items()}, td)}
    report.undeclared = [p for p in param_names if p not in declared]
    species_vars = {v: m for v, m in zip(variables, meta["variables"]) if str(v).startswith("x_")}
    declared_vars = {sp.Symbol(name): d for name, d in getattr(net, "variables", {}).items()}
    by_grid = {"n_Htot", "Δt", "u_initial"} | {f"{v}_initial" for v in species_vars}  # set per point below
    values, missing = _parameter_values(net, param_names, td, {**{p: 0.0 for p in by_grid}, **params})
    report.missing_values = [p for p in missing if p not in report.undeclared]

    # the grid
    sname = {v: sanitize_name(str(v)) for v in species_vars}
    svalues = {sanitize_name(k): v for k, v in values.items()}
    comps = _compositions(meta, {sname[v]: m for v, m in species_vars.items()}, svalues, grid)
    T, n, c = np.meshgrid(grid.T, grid.n_Htot, np.arange(len(comps)), indexing="ij")
    T, n, c = T.ravel(), n.ravel(), c.ravel()
    npts = len(T)
    state = {sp.Symbol("T"): T}
    for v, d in declared_vars.items():  # e.g. a dust temperature: at the gas temperature, within its bounds
        state[v] = np.clip(T, d.floor * (1 + 1e-6), d.ceiling)
    for v in species_vars:
        state[v] = np.array([comps[i][1][sname[v]] for i in range(len(comps))])[c]
    pvals = {p: np.full(npts, values[p]) for p in param_names}
    pvals["n_Htot"] = n
    if "Δt" in pvals:
        pvals["Δt"] = np.full(npts, grid.dt)
    for v in species_vars:
        if f"{v}_initial" in pvals:
            pvals[f"{v}_initial"] = state[v]
    report.n_points = npts
    report.n_expressions = len(exprs)

    # lambdify (sanitized names; argument order: variables, parameters, prelude symbols)
    tables = tables_in(exprs + [W for _, W in prelude])
    var_args = [sanitize_symbols(v) for v in variables]
    par_syms = [sp.Symbol(p) for p in param_names]
    par_args = [sanitize_symbols(s) for s in par_syms]
    pre_args = [sanitize_symbols(w) for w, _ in prelude]
    pre_funcs = [_c_lambdify(var_args + par_args + pre_args[:k], [sanitize_symbols(W)], tables)
                 for k, (_, W) in enumerate(prelude)]
    all_args = var_args + par_args + pre_args
    main = [e for _, e in rows + outputs + system]
    f_main = _c_lambdify(all_args, [sanitize_symbols(e) for e in main], tables)
    f_part = _c_lambdify(all_args, [sanitize_symbols(d) for _, _, d in partials], tables) if partials else None
    f_jac = _c_lambdify(all_args, [sanitize_symbols(x) for r in jac for x in r], tables)

    def evaluate(st, funcs):
        vals = [st.get(v, np.zeros(npts)) for v in variables]
        pv = [pvals[p] for p in param_names]
        pre = []
        for f in pre_funcs:
            pre.append(f(*(vals + pv + pre))[0])
        return [fn(*(vals + pv + pre)) if fn else [] for fn in funcs]

    if u_sym in variables:  # u = u(T, x), the EOS row at u = 0; u_initial = u: no backward-Euler energy term
        (m0,) = evaluate({**state, u_sym: np.zeros(npts)}, [f_main])
        state[u_sym] = m0[len(rows) + len(outputs) + j_of[u_sym]]
        if "u_initial" in pvals:
            pvals["u_initial"] = state[u_sym]
    base_main, base_part, base_jac = evaluate(state, [f_main, f_part, f_jac])

    def example(mask):
        i = int(np.flatnonzero(mask)[0])
        out = {"T": T[i], "n_Htot": n[i]}
        out.update({str(v): state[v][i] for v in species_vars})
        return out

    # non-finite values and partials
    labels = [w for w, _ in rows + outputs + system]
    for what, val in zip(labels, base_main):
        bad = ~np.isfinite(val)
        if bad.any():
            report.nonfinite.append(Finding(what, "value", int(bad.sum()), npts, example=example(bad)))
    row_vals = dict(zip(labels, base_main))
    for (what, j, _), val in zip(partials, base_part):
        bad = ~np.isfinite(val) & np.isfinite(row_vals[what])
        if bad.any():
            report.nonfinite.append(Finding(what, f"d/d{variables[j]}", int(bad.sum()), npts, example=example(bad)))
    nv = len(variables)
    for k, val in enumerate(base_jac):
        i, j = divmod(k, nv)
        bad = ~np.isfinite(val) & np.isfinite(base_main[len(rows) + len(outputs) + i])
        if bad.any():
            report.nonfinite.append(Finding(f"system row {sanitize_name(str(variables[i]))}",
                                            f"Jacobian d/d{variables[j]}", int(bad.sum()), npts, example=example(bad)))

    # Finite differences, for the rates' partials and the system's Jacobian, in T and the species. An error matters
    # against the gross magnitude of its row: the sum of the magnitudes of every process's contribution to it, plus
    # the backward-Euler term; a derivative of a clamped exp(-500) is wrong but irrelevant.
    by_name = {sname[v]: v for v in species_vars}
    row_keys = [key for p in model.subprocesses for key, eq in p.network.items() if sp.sympify(eq.rhs) != 0]
    gross = {}
    for key, val in zip(row_keys, base_main[:len(rows)]):
        gross[key] = gross.get(key, 0.0) + np.abs(val)
    for key, kind in net.species_kinds.items():  # an energy reservoir exchanges with the gas: its scale is the heat's
        if kind == "energy" and key in gross:
            gross[key] = gross[key] + gross.get("heat", 0.0)
    sys_gross = []
    for v in variables:
        if v in declared_vars:
            sys_gross.append(gross.get(declared_vars[v].row, 0.0))
        elif v == u_sym:
            sys_gross.append(np.abs(state[u_sym]))
        elif v == sp.Symbol("T"):
            g = gross.get("heat", 0.0)
            if "Δt" in pvals and u_sym in state:  # rho (u - u_initial)/dt, rho ~ 1.4 m_H n_Htot
                g = g + 2.34e-24 * n * np.abs(state[u_sym]) / pvals["Δt"]
            sys_gross.append(g)
        else:
            species = str(v)[2:]
            g = gross.get(species, 0.0)
            if f"{v}_initial" in pvals:
                g = g + n * np.abs(state[v]) / pvals["Δt"]
            sys_gross.append(g)

    def room(v):
        """(what v's budgets leave, in units of v; their largest total, in units of v)"""
        r, big = np.full(npts, np.inf), 0.0
        for b in meta["budgets"]:
            total = b["total"] if b["total_param"] is None else svalues[b["total_param"]]
            val = total - sum(w * state[by_name[t]] for t, w in b["terms"])
            for t, w in b["terms"]:
                if t == sname[v]:
                    r, big = np.fmin(r, val / w), max(big, total / w)
        return r, big

    n_sys0 = len(rows) + len(outputs)
    labels_sys = labels[n_sys0:n_sys0 + nv]
    # An eliminated abundance total - sum w x carries the round-off of its total: relative precision 1e-16 total/value.
    # Near a budget's end every rate reading it is that noisy, so the absolute tolerance grows by that factor.
    cond = np.ones(npts)
    for b in meta["budgets"]:
        total = b["total"] if b["total_param"] is None else svalues[b["total_param"]]
        val = total - sum(w * state[by_name[t]] for t, w in b["terms"])
        with np.errstate(all="ignore"):
            cond = np.fmax(cond, np.where(val > 0, total / val, np.inf))
    part_index = {(what, j): k for k, (what, j, _) in enumerate(partials)}
    jac_bad = {}
    for j in free:
        v = variables[j]
        x = state[v]
        comparable, r = np.ones(npts, bool), np.full(npts, np.inf)
        if v in species_vars:
            # A step small against the variable and against what its budgets leave (h_lin) is lost in a budget's
            # total (total - w x rounds) once below ~1e-9 of it (h_round). Take h_lin if it is resolved; else h_round
            # if that is still small (curvature errors below 1e-3); else h_lin where the budget path it misses is below
            # the tolerance (what remains >= 1e7 x the variable); else, next to a budget's end, do not compare.
            floor = species_vars[v]["floor"]  # a band's floor is 0: step on its scale there
            s_j = np.abs(x) + (floor if floor > 0 else 1e-12 * species_vars[v]["scale"]) + 1e-300
            r, big = room(v)
            small = np.fmin(s_j, r)
            h_lin, h_round = FD_STEP * small, 1e-3 * FD_STEP * big
            use_round = (h_lin < h_round) & (h_round <= 1e-3 * small)
            h = np.where(use_round, h_round, h_lin)
            comparable = (h_lin >= h_round) | use_round | (r >= 1e7 * s_j)
        else:
            s_j = np.abs(x)
            h = FD_STEP * s_j
        cases = [(labels[k], k, base_part[part_index[(labels[k], j)]], gross[row_keys[k]], f"d/d{v}")
                 for k in range(len(rows)) if (labels[k], j) in part_index]
        cases += [(labels_sys[i], n_sys0 + i, base_jac[i * nv + j], sys_gross[i], f"Jacobian d/d{v}")
                  for i in range(nv)]
        # at steps h and 10 h: round-off in a rate (e.g. a cancelling difference) grows as 1/h, truncation as h^2,
        # while a wrong partial disagrees with both
        verdicts = []
        for hk in (h, 10 * h):
            hk = (x + hk) - x  # the step as represented
            can_up, can_down = r - hk >= 0, x - hk >= 0
            (m_up,), (m_dn,) = evaluate({**state, v: x + hk}, [f_main]), evaluate({**state, v: x - hk}, [f_main])
            out = []
            for what, k_row, d, g, quantity in cases:
                f0, fu, fd = base_main[k_row], m_up[k_row], m_dn[k_row]
                with np.errstate(all="ignore"):
                    fwd, bwd = (fu - f0) / hk, (f0 - fd) / hk
                    up_ok, dn_ok = can_up & np.isfinite(fwd), can_down & np.isfinite(bwd)
                    fdv = np.where(up_ok & dn_ok, 0.5 * (fwd + bwd), np.where(up_ok, fwd, bwd))
                    tol = rtol * np.abs(fdv) + FD_ATOL * cond * g / s_j
                    lo, hi = np.fmin(fwd, bwd) - tol, np.fmax(fwd, bwd) + tol
                    err = np.abs(d - fdv)
                    bad = comparable & (hk > 0) & (up_ok | dn_ok) & np.isfinite(d) & (err > tol) & \
                        ~(up_ok & dn_ok & (d >= lo) & (d <= hi))
                    out.append((bad, err / tol))
            verdicts.append(out)
        for (what, _, _, _, quantity), (bad1, e1), (bad10, e10) in zip(cases, verdicts[0], verdicts[1]):
            bad = bad1 & bad10
            if bad.any():
                key = (what, quantity)
                if key not in jac_bad:
                    worst = float(np.max(np.fmin(e1[bad], e10[bad])))
                    jac_bad[key] = Finding(what, quantity, int(bad.sum()), npts, worst, example(bad))
    report.jacobian = list(jac_bad.values())
    return report
