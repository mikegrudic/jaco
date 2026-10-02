"""Lambdify sympy expressions with the floating-point semantics of jaco's generated C code.

Matches JacoCCodePrinter: exp() arguments clamped to [-500, 500], Heaviside(x) -> (x > 0), Max/Min -> fmax/fmin,
and IEEE double arithmetic (underflow to 0, pow(0, negative) = inf, 0 * inf = nan). Not valid for table lookups.
"""

import numpy as np
import sympy as sp
from sympy.printing.numpy import NumPyPrinter


class _CLikePrinter(NumPyPrinter):
    def _print_exp(self, expr):
        return f"numpy.exp(numpy.clip({self._print(expr.args[0])}, -500.0, 500.0))"

    def _print_Heaviside(self, expr):
        return f"numpy.greater({self._print(expr.args[0])}, 0).astype(float)"

    def _nested(self, func, args):
        out = self._print(args[0])
        for a in args[1:]:
            out = f"numpy.{func}({out}, {self._print(a)})"
        return out

    def _print_Max(self, expr):
        return self._nested("fmax", expr.args)

    def _print_Min(self, expr):
        return self._nested("fmin", expr.args)


def c_lambdify(args, expr):
    f = sp.lambdify(args, expr, modules="numpy", printer=_CLikePrinter)

    def wrapped(*vals):
        with np.errstate(all="ignore"):
            return np.asarray(f(*[np.float64(v) for v in vals]), dtype=float)

    return wrapped
