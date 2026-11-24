"""jaco's interpolants as array functions, in numpy or JAX, with the conventions of the generated code (jaco_interp1d,
jaco_interp1d_const, jaco_table.h): a 1D interpolant held at its end values or extrapolated linearly beyond them, a
2D/3D table query clamped to the table's edge with a zero derivative along that axis. InterpolantPrinter prints the
symbolic nodes as calls to them; jax_lambdify is sympy's lambdify for JAX with those calls available."""

import itertools

import jax.numpy as jnp
import numpy as np
import sympy as sp
from sympy.printing.numpy import JaxPrinter

from ..interpolation import (PiecewiseConstantInterp, PiecewiseLinearInterp, TableInterp2D, TableInterp2D_dx,
                             TableInterp2D_dy, TableInterp3D, TableInterp3D_dx, TableInterp3D_dy, TableInterp3D_dz,
                             tables_in)

__all__ = ["interp1d", "interp1d_const", "table_function", "InterpolantPrinter", "has_interpolants", "jax_lambdify"]


def _float(x, xp):
    return np.asarray(x, dtype=float) if xp is np else xp.asarray(x)  # JAX: its default float, 32 bits unless x64


def interp1d(x, X, Y, extrap, xp=np):
    """jaco_interp1d"""
    X, Y, x = _float(X, xp), _float(Y, xp), _float(x, xp)
    out = xp.interp(x, X, Y)
    if extrap:
        out = xp.where(x <= X[0], Y[0] + (Y[1] - Y[0]) / (X[1] - X[0]) * (x - X[0]), out)
        out = xp.where(x >= X[-1], Y[-1] + (Y[-1] - Y[-2]) / (X[-1] - X[-2]) * (x - X[-1]), out)
    return out


def interp1d_const(x, X, V, xp=np):
    """jaco_interp1d_const: the value on the interval holding x, the end values beyond the ends"""
    X, V, x = _float(X, xp), _float(V, xp), _float(x, xp)
    i = xp.clip(xp.searchsorted(X, x, side="right") - 1, 0, len(V) - 1)
    return V[i]


def table_function(tables, xp=np):
    """jaco_table2d_eval / jaco_table3d_eval (jaco_table.h), vectorized, for the tables of tables_in(): the returned
    table(xs, name, grad) gives the value for grad = -1, else the derivative along axis grad."""
    def table(xs, name, grad):
        t = tables[name]
        nd = t["ndim"]
        xs = xp.broadcast_arrays(*[_float(x, xp) for x in xs])
        f, cl, d, idx = [], [], [], []
        for a in range(nd):
            x = xs[a]
            log = bool(t["log_axes"][a])
            v = xp.log(x) if log else x
            lo = np.log(t[f"axis{a}_min"]) if log else t[f"axis{a}_min"]
            hi = np.log(t[f"axis{a}_max"]) if log else t[f"axis{a}_max"]
            n = t["shape"][a]
            da = (hi - lo) / (n - 1)
            fa = (v - lo) / da
            clamped = (fa < 0) | (fa > n - 1)
            fa = xp.clip(fa, 0, n - 1)
            ia = xp.minimum(fa.astype(int), n - 2)
            f.append(fa - ia)
            idx.append(ia)
            cl.append(clamped)
            d.append(da * (x if log else 1.0))
        data = _float(t["data"], xp)
        val = 0.0
        for corner in itertools.product((0, 1), repeat=nd):
            w = 1.0
            for a in range(nd):
                if a == grad:
                    w = w * (1.0 if corner[a] else -1.0)
                else:
                    w = w * (f[a] if corner[a] else 1 - f[a])
            val = val + w * data[tuple(idx[a] + corner[a] for a in range(nd))]
        if grad < 0:
            return val
        return xp.where(cl[grad], 0.0, val / d[grad])
    return table


class InterpolantPrinter:
    """Printer mixin: the interpolant nodes as calls to _jaco_interp1d, _jaco_interp1d_const and _jaco_table"""

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


INTERPOLANTS = (PiecewiseLinearInterp, PiecewiseConstantInterp, TableInterp2D, TableInterp2D_dx, TableInterp2D_dy,
                TableInterp3D, TableInterp3D_dx, TableInterp3D_dy, TableInterp3D_dz)


def has_interpolants(expr):
    return bool(expr.atoms(*INTERPOLANTS))


def interpolant_namespace(tables, xp=np):
    """the functions InterpolantPrinter's output calls, for the tables of tables_in()"""
    return {"_jaco_interp1d": lambda x, X, Y, extrap: interp1d(x, X, Y, extrap, xp),
            "_jaco_interp1d_const": lambda x, X, V: interp1d_const(x, X, V, xp),
            "_jaco_table": table_function(tables, xp)}


class _JaxPrinter(InterpolantPrinter, JaxPrinter):
    pass


def jax_lambdify(args, expr, **kwargs):
    """sympy.lambdify(args, expr, modules="jax") that also evaluates jaco's interpolants and tables. They are
    differentiable: jax.jacfwd of a table lookup is the derivative jaco_table.h returns."""
    flat = [e for e in sp.flatten([expr]) if isinstance(e, sp.Basic)]
    printer = _JaxPrinter({"fully_qualified_modules": False, "inline": True, "allow_unknown_functions": True,
                           "user_functions": {}})  # lambdify's own settings for its default printer
    return sp.lambdify(args, expr, modules=[interpolant_namespace(tables_in(flat), jnp), "jax"], printer=printer,
                       **kwargs)
