"""Custom sympy Function classes for piecewise interpolation that generate clean code."""

import hashlib
import itertools

import sympy as sp
from sympy.core.symbol import Str
import numpy as np
from scipy.interpolate import RegularGridInterpolator


# --------------------------------------------------------------------------- #
# Tables for multi-dimensional interpolation
# --------------------------------------------------------------------------- #

class Table(Str):
    """A uniformly spaced N-d table as a sympy atom: its name is what the generated code calls it, and it carries its
    data, so code generation finds a model's tables in its expressions (``expr.atoms(Table)``).

    Two tables are equal only if their names and contents are; tables are ordered by creation (``serial``).

    Parameters
    ----------
    name : str
        Identifier in the generated code.
    data : np.ndarray
        N-dimensional array of table values (e.g. shape (n_axis0, n_axis1) for 2D).
    axes : list of np.ndarray
        1D arrays for each axis, in order, uniformly spaced in linear or log space.
    log_axes : list of bool, optional
        Whether each axis is log-spaced. Defaults to False for all axes.
    """

    _serials = itertools.count()

    def __new__(cls, name, data, axes, log_axes=None):
        obj = super().__new__(cls, name)
        data = np.ascontiguousarray(data, dtype=np.float64)
        ndim = data.ndim
        if len(axes) != ndim:
            raise ValueError(f"Expected {ndim} axes for {ndim}D data, got {len(axes)}")
        log_axes = [False] * ndim if log_axes is None else list(log_axes)
        info = {"data": data, "axes": [np.asarray(a, dtype=np.float64) for a in axes], "log_axes": log_axes,
                "ndim": ndim, "shape": data.shape}
        for i, (ax, is_log) in enumerate(zip(axes, log_axes)):
            vals = np.log(ax) if is_log else np.array(ax)
            diffs = np.diff(vals)
            if not np.allclose(diffs, diffs[0], rtol=1e-10):
                raise ValueError(f"Axis {i} is not uniformly spaced (in {'log' if is_log else 'linear'} space)")
            info[f"axis{i}_min"] = float(ax[0])
            info[f"axis{i}_max"] = float(ax[-1])
        h = hashlib.sha256(np.asarray(data.shape, dtype=np.int64).tobytes())
        for a, is_log in zip(info["axes"], log_axes):
            h.update(a.tobytes())
            h.update(bytes([int(is_log)]))
        h.update(data.tobytes())
        obj.info = info
        obj.digest = h.hexdigest()
        obj.serial = next(cls._serials)
        return obj

    def __getnewargs__(self):
        return (self.name, self.info["data"], self.info["axes"], self.info["log_axes"])

    def _hashable_content(self):
        return (self.name, self.digest)

    def save_hdf5(self, filename=None):
        """Save the table to an HDF5 file (default '{name}.h5')."""
        import h5py
        with h5py.File(filename or f"{self.name}.h5", "w") as f:
            f.create_dataset("data", data=self.info["data"])
            for i, ax in enumerate(self.info["axes"]):
                f.create_dataset(f"axis{i}", data=ax)
            f.attrs["ndim"] = self.info["ndim"]
            for i, is_log in enumerate(self.info["log_axes"]):
                f.attrs[f"axis{i}_log"] = int(is_log)


def tables_in(exprs):
    """name -> table info of every Table in exprs, in creation order; different tables of one name raise"""
    found = set()
    for e in exprs:
        found |= sp.sympify(e).atoms(Table)
    tables = {}
    for t in sorted(found, key=lambda t: t.serial):
        if t.name in tables:
            raise ValueError(f"two different tables are named {t.name}")
        tables[t.name] = t.info
    return tables


class PiecewiseLinearInterp(sp.Function):
    """Piecewise-linear interpolation of tabulated data.

    Arguments: (x, X_tuple, Y_tuple, extrapolate_flag[, name])
        x: sympy expression for the interpolation variable
        X_tuple: sympy.Tuple of Float values (breakpoints, strictly increasing)
        Y_tuple: sympy.Tuple of Float values (function values at breakpoints)
        extrapolate_flag: sympy.Integer(0) to clamp, sympy.Integer(1) to extrapolate
        name: sympy.core.symbol.Str, descriptive identifier for code generation (optional)
    """
    is_commutative = True  # scalar-valued; the Str/Tuple arguments otherwise leave this undetermined

    @classmethod
    def eval(cls, x, X_tuple, Y_tuple, extrapolate_flag, name=None):
        # If x is a pure number, evaluate numerically
        if x.is_Number:
            return cls._evaluate(float(x), X_tuple, Y_tuple, int(extrapolate_flag))
        return None

    @staticmethod
    def _evaluate(xval, X_tuple, Y_tuple, extrapolate):
        """Numerically evaluate the piecewise-linear interpolant at xval."""
        X = [float(v) for v in X_tuple]
        Y = [float(v) for v in Y_tuple]
        n = len(X)
        if xval <= X[0]:
            if extrapolate:
                slope = (Y[1] - Y[0]) / (X[1] - X[0])
                return sp.Float(Y[0] + slope * (xval - X[0]))
            return sp.Float(Y[0])
        if xval >= X[-1]:
            if extrapolate:
                slope = (Y[-1] - Y[-2]) / (X[-1] - X[-2])
                return sp.Float(Y[-1] + slope * (xval - X[-1]))
            return sp.Float(Y[-1])
        # Binary search
        lo, hi = 0, n - 1
        while hi - lo > 1:
            mid = (lo + hi) // 2
            if X[mid] <= xval:
                lo = mid
            else:
                hi = mid
        t = (xval - X[lo]) / (X[hi] - X[lo])
        return sp.Float(Y[lo] + t * (Y[hi] - Y[lo]))

    @property
    def table_name(self):
        """Return the descriptive name string, or None if not set."""
        if len(self.args) > 4:
            return str(self.args[4])
        return None

    def _eval_derivative(self, s):
        x = self.args[0]
        if not x.has(s):
            return sp.S.Zero
        X_tuple = self.args[1]
        Y_tuple = self.args[2]
        extrap = self.args[3]
        X = [float(v) for v in X_tuple]
        Y = [float(v) for v in Y_tuple]
        slopes = [(Y[i + 1] - Y[i]) / (X[i + 1] - X[i]) for i in range(len(X) - 1)]
        slopes_tuple = sp.Tuple(*[sp.Float(s_) for s_ in slopes])
        name = self.table_name
        deriv_args = [x, X_tuple, slopes_tuple, extrap]
        if name is not None:
            deriv_args.append(Str(name + "_slopes"))
        return PiecewiseConstantInterp(*deriv_args) * x.diff(s)

    def _eval_evalf(self, prec):
        x = self.args[0]
        if x.is_Number:
            return self._evaluate(
                float(x), self.args[1], self.args[2], int(self.args[3])
            )
        return self

    def doit(self, **hints):
        """Expand to a sympy Piecewise expression (fallback for printers that don't support this type)."""
        x = self.args[0]
        X = [float(v) for v in self.args[1]]
        Y = [float(v) for v in self.args[2]]
        extrapolate = bool(self.args[3])

        slopes = np.diff(Y) / np.diff(X)
        cases = []
        if extrapolate:
            cases.append((Y[0] + slopes[0] * (x - X[0]), x < X[0]))
            cases.append((Y[-1] + slopes[-1] * (x - X[-1]), x >= X[-1]))
        else:
            cases.append((sp.Float(Y[0]), x < X[0]))
            cases.append((sp.Float(Y[-1]), x >= X[-1]))

        for i in range(len(X) - 1):
            cases.append((Y[i] + slopes[i] * (x - X[i]), (x >= X[i]) & (x < X[i + 1])))
        return sp.Piecewise(*cases)


class PiecewiseConstantInterp(sp.Function):
    """Piecewise-constant interpolation (used for derivatives of piecewise-linear interpolants).

    Arguments: (x, X_tuple, values_tuple, extrapolate_flag[, name])
        x: sympy expression for the interpolation variable
        X_tuple: sympy.Tuple of Float values (breakpoints)
        values_tuple: sympy.Tuple of Float values (constant values on each interval, length = len(X)-1)
        extrapolate_flag: sympy.Integer(0) or sympy.Integer(1)
        name: sympy.core.symbol.Str, descriptive identifier for code generation (optional)
    """
    is_commutative = True  # scalar-valued; the Str/Tuple arguments otherwise leave this undetermined

    @property
    def table_name(self):
        """Return the descriptive name string, or None if not set."""
        if len(self.args) > 4:
            return str(self.args[4])
        return None

    @classmethod
    def eval(cls, x, X_tuple, values_tuple, extrapolate_flag, name=None):
        if x.is_Number:
            return cls._evaluate(float(x), X_tuple, values_tuple)
        return None

    @staticmethod
    def _evaluate(xval, X_tuple, values_tuple):
        """Numerically evaluate the piecewise-constant interpolant at xval."""
        X = [float(v) for v in X_tuple]
        vals = [float(v) for v in values_tuple]
        n = len(X)
        if xval <= X[0]:
            return sp.Float(vals[0])
        if xval >= X[-1]:
            return sp.Float(vals[-1])
        # Binary search
        lo, hi = 0, n - 1
        while hi - lo > 1:
            mid = (lo + hi) // 2
            if X[mid] <= xval:
                lo = mid
            else:
                hi = mid
        return sp.Float(vals[lo])

    def _eval_derivative(self, s):
        # Derivative of piecewise-constant is zero
        return sp.S.Zero

    def _eval_evalf(self, prec):
        x = self.args[0]
        if x.is_Number:
            return self._evaluate(float(x), self.args[1], self.args[2])
        return self

    def doit(self, **hints):
        """Expand to a sympy Piecewise expression."""
        x = self.args[0]
        X = [float(v) for v in self.args[1]]
        vals = [float(v) for v in self.args[2]]
        extrapolate = bool(self.args[3])

        cases = []
        # Below first breakpoint
        cases.append((sp.Float(vals[0]), x < X[0]))
        # Above last breakpoint
        cases.append((sp.Float(vals[-1]), x >= X[-1]))
        # Interior intervals
        for i in range(len(X) - 1):
            cases.append((sp.Float(vals[i]), (x >= X[i]) & (x < X[i + 1])))
        return sp.Piecewise(*cases)


# --------------------------------------------------------------------------- #
# Multi-dimensional table interpolation
# --------------------------------------------------------------------------- #

class TableInterp2D(sp.Function):
    """Bilinear interpolation on a 2D regular grid.

    Arguments: (x, y, table)
        x, y: sympy expressions for the interpolation variables
        table: the Table
    """
    is_commutative = True  # scalar-valued; the Str/Tuple arguments otherwise leave this undetermined

    @classmethod
    def eval(cls, x, y, table):
        if x.is_Number and y.is_Number:
            return cls._evaluate(float(x), float(y), table)
        return None

    @staticmethod
    def _evaluate(xval, yval, table):
        table = table.info
        interp = RegularGridInterpolator(
            tuple(table["axes"]), table["data"],
            method="linear", bounds_error=False, fill_value=None,
        )
        return sp.Float(float(interp([[xval, yval]])[0]))

    @property
    def table_name(self):
        return str(self.args[2])

    def _eval_derivative(self, s):
        x, y = self.args[0], self.args[1]
        table = self.args[2]
        result = sp.S.Zero
        if x.has(s):
            result += TableInterp2D_dx(x, y, table) * x.diff(s)
        if y.has(s):
            result += TableInterp2D_dy(x, y, table) * y.diff(s)
        return result

    def _eval_evalf(self, prec):
        x, y = self.args[0], self.args[1]
        if x.is_Number and y.is_Number:
            return self._evaluate(float(x), float(y), self.args[2])
        return self


class TableInterp2D_dx(sp.Function):
    """Partial derivative of TableInterp2D w.r.t. its first argument (x).

    Computed analytically from the bilinear formula — no extra table needed.
    """
    is_commutative = True  # scalar-valued; the Str/Tuple arguments otherwise leave this undetermined

    @classmethod
    def eval(cls, x, y, table):
        if x.is_Number and y.is_Number:
            return cls._evaluate(float(x), float(y), table)
        return None

    @staticmethod
    def _evaluate(xval, yval, table):
        """Evaluate x-partial derivative via finite difference of the bilinear interp."""
        info = table.info
        eps = (info["axis0_max"] - info["axis0_min"]) / (info["shape"][0] - 1) * 1e-6
        f_plus = float(TableInterp2D._evaluate(xval + eps, yval, table))
        f_minus = float(TableInterp2D._evaluate(xval - eps, yval, table))
        return sp.Float((f_plus - f_minus) / (2 * eps))

    @property
    def table_name(self):
        return str(self.args[2])

    def _eval_derivative(self, s):
        return sp.S.Zero  # second derivatives not supported


class TableInterp2D_dy(sp.Function):
    """Partial derivative of TableInterp2D w.r.t. its second argument (y)."""
    is_commutative = True  # scalar-valued; the Str/Tuple arguments otherwise leave this undetermined

    @classmethod
    def eval(cls, x, y, table):
        if x.is_Number and y.is_Number:
            return cls._evaluate(float(x), float(y), table)
        return None

    @staticmethod
    def _evaluate(xval, yval, table):
        info = table.info
        eps = (info["axis1_max"] - info["axis1_min"]) / (info["shape"][1] - 1) * 1e-6
        f_plus = float(TableInterp2D._evaluate(xval, yval + eps, table))
        f_minus = float(TableInterp2D._evaluate(xval, yval - eps, table))
        return sp.Float((f_plus - f_minus) / (2 * eps))

    @property
    def table_name(self):
        return str(self.args[2])

    def _eval_derivative(self, s):
        return sp.S.Zero


class TableInterp3D(sp.Function):
    """Trilinear interpolation on a 3D regular grid.

    Arguments: (x, y, z, table)
    """
    is_commutative = True  # scalar-valued; the Str/Tuple arguments otherwise leave this undetermined

    @classmethod
    def eval(cls, x, y, z, table):
        if x.is_Number and y.is_Number and z.is_Number:
            return cls._evaluate(float(x), float(y), float(z), table)
        return None

    @staticmethod
    def _evaluate(xval, yval, zval, table):
        table = table.info
        interp = RegularGridInterpolator(
            tuple(table["axes"]), table["data"],
            method="linear", bounds_error=False, fill_value=None,
        )
        return sp.Float(float(interp([[xval, yval, zval]])[0]))

    @property
    def table_name(self):
        return str(self.args[3])

    def _eval_derivative(self, s):
        x, y, z = self.args[0], self.args[1], self.args[2]
        table = self.args[3]
        result = sp.S.Zero
        if x.has(s):
            result += TableInterp3D_dx(x, y, z, table) * x.diff(s)
        if y.has(s):
            result += TableInterp3D_dy(x, y, z, table) * y.diff(s)
        if z.has(s):
            result += TableInterp3D_dz(x, y, z, table) * z.diff(s)
        return result

    def _eval_evalf(self, prec):
        x, y, z = self.args[0], self.args[1], self.args[2]
        if x.is_Number and y.is_Number and z.is_Number:
            return self._evaluate(float(x), float(y), float(z), self.args[3])
        return self


class TableInterp3D_dx(sp.Function):
    """Partial derivative of TableInterp3D w.r.t. x."""
    is_commutative = True  # scalar-valued; the Str/Tuple arguments otherwise leave this undetermined

    @classmethod
    def eval(cls, x, y, z, table):
        return None
    @property
    def table_name(self):
        return str(self.args[3])
    def _eval_derivative(self, s):
        return sp.S.Zero


class TableInterp3D_dy(sp.Function):
    """Partial derivative of TableInterp3D w.r.t. y."""
    is_commutative = True  # scalar-valued; the Str/Tuple arguments otherwise leave this undetermined

    @classmethod
    def eval(cls, x, y, z, table):
        return None
    @property
    def table_name(self):
        return str(self.args[3])
    def _eval_derivative(self, s):
        return sp.S.Zero


class TableInterp3D_dz(sp.Function):
    """Partial derivative of TableInterp3D w.r.t. z."""
    is_commutative = True  # scalar-valued; the Str/Tuple arguments otherwise leave this undetermined

    @classmethod
    def eval(cls, x, y, z, table):
        return None
    @property
    def table_name(self):
        return str(self.args[3])
    def _eval_derivative(self, s):
        return sp.S.Zero
