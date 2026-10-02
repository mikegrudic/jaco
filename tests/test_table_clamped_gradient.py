"""jaco_table.h: a query clamped to the table edge must return the edge value with a
zero gradient along the clamped axis (the interpolant is constant there), and interior
gradients must match finite differences."""

import ctypes
import subprocess
import tempfile
from pathlib import Path

import numpy as np
import pytest

HEADER = Path(__file__).resolve().parents[1] / "src" / "jaco" / "codegen" / "jaco_table.h"

DRIVER = r"""
#include "jaco_table.h"
double eval2d(double x0, double x1, double *data, int n0, int n1, double *g) {
    JacoTable2D t = {data, 1.0, 100.0, 1.0, 1000.0, n0, n1, 1, 1};
    return jaco_table2d_eval(x0, x1, &t, &g[0], &g[1]);
}
double eval3d(double x0, double x1, double x2, double *data, int n0, int n1, int n2, double *g) {
    JacoTable3D t = {data, 1.0, 100.0, 1.0, 1000.0, 1.0, 10.0, n0, n1, n2, 1, 1, 1};
    return jaco_table3d_eval(x0, x1, x2, &t, &g[0], &g[1], &g[2]);
}
"""


@pytest.fixture(scope="module")
def lib():
    d = Path(tempfile.mkdtemp())
    (d / "driver.c").write_text(DRIVER)
    so = d / "driver.so"
    subprocess.run(["gcc", "-O1", "-shared", "-fPIC", f"-I{HEADER.parent}", "-o", str(so), str(d / "driver.c"), "-lm"], check=True)
    L = ctypes.CDLL(str(so))
    dp = ctypes.POINTER(ctypes.c_double)
    L.eval2d.restype = ctypes.c_double
    L.eval2d.argtypes = [ctypes.c_double, ctypes.c_double, dp, ctypes.c_int, ctypes.c_int, dp]
    L.eval3d.restype = ctypes.c_double
    L.eval3d.argtypes = [ctypes.c_double, ctypes.c_double, ctypes.c_double, dp, ctypes.c_int, ctypes.c_int, ctypes.c_int, dp]
    return L


def _ptr(a):
    return a.ctypes.data_as(ctypes.POINTER(ctypes.c_double))


def test_2d_clamped_axis_has_zero_gradient(lib):
    n0, n1 = 5, 7
    x0 = np.exp(np.linspace(np.log(1), np.log(100), n0))
    x1 = np.exp(np.linspace(np.log(1), np.log(1000), n1))
    data = np.ascontiguousarray(np.log(x0)[:, None] * 2.0 + np.log(x1)[None, :] ** 2)
    g = np.zeros(2)
    # interior: gradient matches central differences
    v = lib.eval2d(10.0, 50.0, _ptr(data), n0, n1, _ptr(g))
    h = 1e-4
    g0 = (lib.eval2d(10.0 * (1 + h), 50.0, _ptr(data), n0, n1, _ptr(np.zeros(2))) - lib.eval2d(10.0 * (1 - h), 50.0, _ptr(data), n0, n1, _ptr(np.zeros(2)))) / (2 * h * 10.0)
    g1 = (lib.eval2d(10.0, 50.0 * (1 + h), _ptr(data), n0, n1, _ptr(np.zeros(2))) - lib.eval2d(10.0, 50.0 * (1 - h), _ptr(data), n0, n1, _ptr(np.zeros(2)))) / (2 * h * 50.0)
    assert g[0] == pytest.approx(g0, rel=1e-6) and g[1] == pytest.approx(g1, rel=1e-6)
    # below the axis-0 floor: edge value, zero gradient on axis 0, finite nonzero on axis 1
    v_edge = lib.eval2d(1.0, 50.0, _ptr(data), n0, n1, _ptr(np.zeros(2)))
    v_out = lib.eval2d(0.2, 50.0, _ptr(data), n0, n1, _ptr(g))
    assert v_out == v_edge and g[0] == 0.0 and g[1] != 0.0 and np.isfinite(g[1])
    # above the axis-1 ceiling
    v_edge = lib.eval2d(10.0, 1000.0, _ptr(data), n0, n1, _ptr(np.zeros(2)))
    v_out = lib.eval2d(10.0, 5000.0, _ptr(data), n0, n1, _ptr(g))
    assert v_out == v_edge and g[1] == 0.0 and g[0] != 0.0
    del v


def test_3d_clamped_axis_has_zero_gradient(lib):
    n0, n1, n2 = 4, 5, 6
    x0 = np.exp(np.linspace(np.log(1), np.log(100), n0))
    x1 = np.exp(np.linspace(np.log(1), np.log(1000), n1))
    x2 = np.exp(np.linspace(np.log(1), np.log(10), n2))
    data = np.ascontiguousarray(np.log(x0)[:, None, None] + 3.0 * np.log(x1)[None, :, None] - np.log(x2)[None, None, :] ** 2)
    g = np.zeros(3)
    lib.eval3d(10.0, 50.0, 3.0, _ptr(data), n0, n1, n2, _ptr(g))
    assert np.all(g != 0.0) and np.all(np.isfinite(g))
    v_edge = lib.eval3d(10.0, 50.0, 10.0, _ptr(data), n0, n1, n2, _ptr(np.zeros(3)))
    v_out = lib.eval3d(10.0, 50.0, 40.0, _ptr(data), n0, n1, n2, _ptr(g))
    assert v_out == v_edge and g[2] == 0.0 and g[0] != 0.0 and g[1] != 0.0
    v_edge = lib.eval3d(1.0, 50.0, 3.0, _ptr(data), n0, n1, n2, _ptr(np.zeros(3)))
    v_out = lib.eval3d(0.01, 50.0, 3.0, _ptr(data), n0, n1, n2, _ptr(g))
    assert v_out == v_edge and g[0] == 0.0 and g[1] != 0.0 and g[2] != 0.0
