"""The expressions of a model as one compiled Python function, for the explicit right-hand side (jaco.ode).

    src = emit_source(exprs, var_names, param_names, prelude, tables)   # def evaluate(X, params) -> out
    f = compile_source(src, backend)                                   # "numpy" | "numba" | "jax"

The source is the code generator's Python target (jaco.codegen.printers.JacoPythonCodePrinter): common subexpressions
eliminated across every output, the intermediates and their derivatives as a prelude, the 1D interpolations and the
multi-dimensional tables as module-level arrays with the same clamping and derivative conventions as jaco_table.h.
The numba backend compiles that source with nopython and numpy error semantics (a division by zero is an inf, as in
C), from a module written to $JACO_CACHE_DIR (default ~/.cache/jaco) so the machine code is cached across processes;
the jax backend prints a variant with jnp.where in place of conditionals and jits it; the numpy backend runs the
source as plain Python (slow: every scalar operation is a Python call).
"""

import numpy as np
import sympy as sp

from .codegen.printers import JacoPythonCodePrinter, _PYTHON_INTERP_HELPERS
from .equation_system import EquationSystem

__all__ = ["emit_source", "compile_source", "BACKENDS"]

BACKENDS = ("numba", "numpy", "jax")

_TABLE_HELPERS = '''
def jaco_table2d(x, y, data, n0, lo0, d0, log0, n1, lo1, d1, log1, grad):
    """bilinear table lookup (jaco_table.h conventions): axes uniform in linear or log space, a query outside the
    table clamped to the edge with a zero derivative along that axis; grad = -1 for the value, 0 or 1 for d/dx, d/dy"""
    v0 = math.log(x) if log0 else x
    v1 = math.log(y) if log1 else y
    f0 = (v0 - lo0) / d0
    f1 = (v1 - lo1) / d1
    c0 = f0 < 0 or f0 > n0 - 1
    c1 = f1 < 0 or f1 > n1 - 1
    f0 = min(max(f0, 0.0), n0 - 1.0)
    f1 = min(max(f1, 0.0), n1 - 1.0)
    i0 = min(int(f0), n0 - 2)
    i1 = min(int(f1), n1 - 2)
    f0 -= i0
    f1 -= i1
    if grad == 0:
        if c0:
            return 0.0
        val = (1 - f1) * (data[i0 + 1, i1] - data[i0, i1]) + f1 * (data[i0 + 1, i1 + 1] - data[i0, i1 + 1])
        return val / (d0 * (x if log0 else 1.0))
    if grad == 1:
        if c1:
            return 0.0
        val = (1 - f0) * (data[i0, i1 + 1] - data[i0, i1]) + f0 * (data[i0 + 1, i1 + 1] - data[i0 + 1, i1])
        return val / (d1 * (y if log1 else 1.0))
    return ((1 - f0) * (1 - f1) * data[i0, i1] + f0 * (1 - f1) * data[i0 + 1, i1]
            + (1 - f0) * f1 * data[i0, i1 + 1] + f0 * f1 * data[i0 + 1, i1 + 1])
'''

_JAX_HELPERS = '''
def jaco_interp1d(x, xp, yp, n, extrap):
    xp, yp = jnp.asarray(xp), jnp.asarray(yp)
    inside = jnp.interp(x, xp, yp)
    lo = yp[0] + (yp[1] - yp[0]) / (xp[1] - xp[0]) * (x - xp[0]) if extrap else yp[0]
    hi = yp[n - 1] + (yp[n - 1] - yp[n - 2]) / (xp[n - 1] - xp[n - 2]) * (x - xp[n - 1]) if extrap else yp[n - 1]
    return jnp.where(x <= xp[0], lo, jnp.where(x >= xp[n - 1], hi, inside))

def jaco_interp1d_const(x, xp, vp, n, extrap):
    xp, vp = jnp.asarray(xp), jnp.asarray(vp)
    i = jnp.clip(jnp.searchsorted(xp, x, side="right") - 1, 0, n - 2)
    return jnp.where(x <= xp[0], vp[0], jnp.where(x >= xp[n - 1], vp[-1], vp[i]))

def jaco_table2d(x, y, data, n0, lo0, d0, log0, n1, lo1, d1, log1, grad):
    data = jnp.asarray(data)
    v0 = jnp.log(x) if log0 else x
    v1 = jnp.log(y) if log1 else y
    f0 = (v0 - lo0) / d0
    f1 = (v1 - lo1) / d1
    c0 = (f0 < 0) | (f0 > n0 - 1)
    c1 = (f1 < 0) | (f1 > n1 - 1)
    f0 = jnp.clip(f0, 0.0, n0 - 1.0)
    f1 = jnp.clip(f1, 0.0, n1 - 1.0)
    i0 = jnp.minimum(f0.astype(jnp.int32), n0 - 2)
    i1 = jnp.minimum(f1.astype(jnp.int32), n1 - 2)
    f0 = f0 - i0
    f1 = f1 - i1
    a, b, c, d = data[i0, i1], data[i0 + 1, i1], data[i0, i1 + 1], data[i0 + 1, i1 + 1]
    if grad == 0:
        return jnp.where(c0, 0.0, ((1 - f1) * (b - a) + f1 * (d - c)) / (d0 * (x if log0 else 1.0)))
    if grad == 1:
        return jnp.where(c1, 0.0, ((1 - f0) * (c - a) + f0 * (d - b)) / (d1 * (y if log1 else 1.0)))
    return (1 - f0) * (1 - f1) * a + f0 * (1 - f1) * b + (1 - f0) * f1 * c + f0 * f1 * d
'''


class OdePythonPrinter(JacoPythonCodePrinter):
    """The Python target with conditionals and the multi-dimensional tables printed as Python (the shared mixin prints
    them as C), and the table lookups as calls on module-level arrays."""

    def __init__(self, tables, jax=False, **kw):
        super().__init__(**kw)
        self._tables = tables
        self._jax = jax

    def _print_exp(self, expr):
        """the C target's exp: the argument clamped to [-500, 500], so a switched-off fitting branch cannot overflow
        into inf * 0"""
        arg = self._print(expr.args[0])
        if self._jax:
            return f"jnp.exp(jnp.clip({arg}, -500.0, 500.0))"
        return f"math.exp(min(max({arg}, -500.0), 500.0))"

    def _print_Heaviside(self, expr):
        arg = self._print(expr.args[0])
        return f"jnp.where({arg} > 0, 1.0, 0.0)" if self._jax else f"(1.0 if {arg} > 0 else 0.0)"

    def _print_Piecewise(self, expr):
        if expr.args[-1].cond is not sp.true:
            raise ValueError("Piecewise without a default branch")
        result = self._print(expr.args[-1].expr)
        for e, c in reversed(expr.args[:-1]):
            cond, val = self._print(c), self._print(e)
            result = f"jnp.where({cond}, {val}, {result})" if self._jax else f"({val} if {cond} else {result})"
        return result

    def _print_Min(self, expr):
        args = [self._print(a) for a in expr.args]
        if self._jax:
            out = args[0]
            for a in args[1:]:
                out = f"jnp.minimum({out}, {a})"
            return out
        return f"min({', '.join(args)})"

    def _print_Max(self, expr):
        args = [self._print(a) for a in expr.args]
        if self._jax:
            out = args[0]
            for a in args[1:]:
                out = f"jnp.maximum({out}, {a})"
            return out
        return f"max({', '.join(args)})"

    def _table_call(self, expr, grad):
        name = str(expr.args[-1])
        t = self._tables[name]
        if t["ndim"] != 2:
            raise NotImplementedError(f"{t['ndim']}D table {name}: only 2D tables are printed for the Python target")
        args = [self._print(expr.args[0]), self._print(expr.args[1]), f"{name}_data"]
        for a in range(2):
            log = bool(t["log_axes"][a])
            lo = np.log(t[f"axis{a}_min"]) if log else t[f"axis{a}_min"]
            hi = np.log(t[f"axis{a}_max"]) if log else t[f"axis{a}_max"]
            n = t["shape"][a]
            args += [str(n), repr(float(lo)), repr(float((hi - lo) / (n - 1))), str(log)]
        return f"jaco_table2d({', '.join(args)}, {grad})"

    def _print_TableInterp2D(self, expr):
        return self._table_call(expr, -1)

    def _print_TableInterp2D_dx(self, expr):
        return self._table_call(expr, 0)

    def _print_TableInterp2D_dy(self, expr):
        return self._table_call(expr, 1)

    def _print_TableInterp3D(self, expr):
        raise NotImplementedError("3D tables are not printed for the Python target")


def emit_source(exprs, var_names, param_names, prelude, tables, func_name="evaluate", jax=False, numpy_math=False):
    """Python source defining func_name(X, params) -> out (a float array of len(exprs)), X in var_names order and
    params in param_names order; prelude is the ordered (symbol, expression) list of the intermediates and their
    derivatives, tables the dict of tables_in(). Symbols must already be sanitized. numpy_math prints the scalar
    functions as numpy's (for the plain Python backend; math raises on overflow)."""
    printer = OdePythonPrinter(tables, jax=jax)
    ind = "    "
    lines = [f"def {func_name}(X, params):"]
    for i, v in enumerate(var_names):
        lines.append(f"{ind}{v} = X[{i}]")
    for i, p in enumerate(param_names):
        lines.append(f"{ind}{p} = params[{i}]")
    if prelude:
        lines += EquationSystem._prelude_lines(prelude, printer, "", "", ind, True)
    tmps, out = sp.cse(list(exprs), symbols=sp.numbered_symbols("c_"))
    for t, e in tmps:
        lines.append(f"{ind}{printer.doprint(t)} = {printer.doprint(e)}")
    if jax:
        lines.append(f"{ind}return jnp.array([" + ", ".join(printer.doprint(e) for e in out) + "])")
    else:
        lines.append(f"{ind}out = np.empty({len(out)})")
        for i, e in enumerate(out):
            lines.append(f"{ind}out[{i}] = {printer.doprint(e)}")
        lines.append(f"{ind}return out")
    body = "\n".join(lines)
    if jax:
        body = body.replace("math.", "jnp.")
    elif numpy_math:
        body = body.replace("math.", "np.")  # overflow and log(0) give inf with a warning instead of raising

    head = ["import math", "import numpy as np"]
    if jax:
        head += ["import jax.numpy as jnp", _JAX_HELPERS]
    else:
        head += [_PYTHON_INTERP_HELPERS.replace("import numpy as np\n", ""), _TABLE_HELPERS]
    for key, (name, values) in printer._interp_tables.items():
        head.append(f"{name} = np.array([{', '.join(repr(v) for v in values)}])")
    for name, t in tables.items():
        if t["ndim"] == 2:
            head.append(f"{name}_data = np.array({np.asarray(t['data']).tolist()!r})")
    return "\n".join(head) + "\n\n" + body + "\n"


def compile_source(src, backend="numba", func_name="evaluate"):
    """the function defined by the source, compiled for the backend (numba: nopython with numpy error semantics and
    the helpers compiled too; jax: jitted; numpy: as written)"""
    if backend not in BACKENDS:
        raise ValueError(f"backend {backend!r}: one of {BACKENDS}")
    ns = {}
    if backend == "numba":
        # numba caches machine code only for functions defined in a file: the source goes to a module in the cache
        # directory, named by its hash, so a model compiles once per machine rather than once per process
        import hashlib
        import importlib.util
        import os
        src = "import numba\n" + src
        src = src.replace("\ndef jaco_", "\n@numba.njit(cache=True, error_model='numpy')\ndef jaco_")
        src = src.replace(f"\ndef {func_name}(", f"\n@numba.njit(cache=True, error_model='numpy')\ndef {func_name}(")
        cache = os.environ.get("JACO_CACHE_DIR", os.path.join(os.path.expanduser("~"), ".cache", "jaco"))
        os.makedirs(cache, exist_ok=True)
        name = "jaco_ode_" + hashlib.sha256(src.encode()).hexdigest()[:16]
        path = os.path.join(cache, name + ".py")
        if not os.path.exists(path):
            tmp = path + f".{os.getpid()}.tmp"
            with open(tmp, "w") as f:
                f.write(src)
            os.replace(tmp, path)
        import sys
        if name in sys.modules:
            return getattr(sys.modules[name], func_name)
        spec = importlib.util.spec_from_file_location(name, path)
        module = importlib.util.module_from_spec(spec)
        sys.modules[name] = module  # numba's cache loader imports the module by name
        spec.loader.exec_module(module)
        return getattr(module, func_name)
    elif backend == "jax":
        import jax
        jax.config.update("jax_enable_x64", True)
        ns["jax"] = jax
    exec(compile(src, f"<jaco {func_name} {backend}>", "exec"), ns)
    f = ns[func_name]
    if backend == "jax":
        f = ns["jax"].jit(f)
    return f
