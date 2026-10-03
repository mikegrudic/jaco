"""The generated C 1D interpolation helpers: the piecewise-constant one (the derivative of a piecewise-linear
interpolant) returns its last interval's value beyond the last breakpoint, without reading past its n - 1 values."""

import ctypes
import subprocess

from jaco.codegen.printers import _C_INTERP_HELPERS

DRIVER = r"""
struct guarded { double vp[2]; double guard; };
double beyond(double x) {
    static const double xp[3] = {1.0, 2.0, 3.0};
    static const struct guarded v = {{10.0, 20.0}, 999.0};
    return jaco_interp1d_const(x, xp, v.vp, 3, 1);
}
"""


def test_constant_interpolation_stays_in_bounds(tmp_path):
    (tmp_path / "d.c").write_text(_C_INTERP_HELPERS + DRIVER)
    so = tmp_path / "d.so"
    subprocess.run(["gcc", "-O1", "-shared", "-fPIC", "-o", str(so), str(tmp_path / "d.c"), "-lm"], check=True)
    lib = ctypes.CDLL(str(so))
    lib.beyond.restype, lib.beyond.argtypes = ctypes.c_double, [ctypes.c_double]
    assert [lib.beyond(x) for x in (0.5, 1.5, 2.5, 3.0, 7.0)] == [10.0, 10.0, 20.0, 20.0, 20.0]
