"""jaco.check: every rate, its partials, the outputs and the reduced system over a (T, n, x) grid with the generated
C code's float semantics. The shipped models must pass; seeded defects must be found. The shipped models run on a
subsampled grid; JACO_CHECK_FULL=1 runs the standard one (minutes in total)."""

import os
from importlib import import_module

import numpy as np
import pytest
import sympy as sp

import jaco
from jaco.declarations import Parameter, Species
from jaco.model import Model
from jaco.model_check import CheckGrid, _interp1d_const, _table_function
from jaco.processes import CollisionalIonization, GasPhaseRecombination, ThermalTerm

FULL = os.environ.get("JACO_CHECK_FULL") == "1"
T, n_Htot = sp.Symbol("T"), sp.Symbol("n_Htot")
pdv = ThermalTerm(sp.Symbol("pdv_work"), name="PdV work")
GRID = CheckGrid(np.logspace(np.log10(3.0), 9, 6), np.logspace(-3, 9, 3), n_random=2)


def toy(*extra, parameters=()):
    return Model([CollisionalIonization("H"), GasPhaseRecombination("H+"), pdv, *extra], solve_vars=["u", "T", "H+"],
                 time_dependent=["T"], species=[Species("H"), Species("H+"), Species("e-")], parameters=parameters)


@pytest.mark.slow
@pytest.mark.parametrize("name", ["starforge", "starforge_legacy", "wind_comparison"])
def test_shipped_models_pass(name):
    report = jaco.check(import_module(f"jaco.models.{name}").make_model(), "standard" if FULL else GRID, name=name)
    assert report.ok, report.summary()


def test_toy_model_passes():
    report = jaco.check(toy(), "quick")
    assert report.ok, report.summary()
    assert "OK" in report.summary() and report.n_points == 7 * 4 * 4  # floor, trace, at cap, 30% of cap


def test_nonfinite_rate_is_found_with_its_process():
    report = jaco.check(toy(ThermalTerm(-1e-25 * n_Htot**2 * sp.sqrt(T - 100), name="bad sqrt")), "quick")
    assert not report.ok
    found = {(f.what, f.quantity) for f in report.nonfinite}
    assert ("process 'bad sqrt' row heat", "value") in found and ("system row T", "value") in found
    assert all(f.example["T"] < 100 for f in report.nonfinite)


def test_underflow_times_infinity_is_found():
    """y log y -> 0 as y -> 0 in exact arithmetic, but y = x_H+^20 underflows to 0 at a floored x_H+ and the C code
    computes 0 * -inf = nan"""
    y = sp.Symbol("x_H+") ** 20
    report = jaco.check(toy(ThermalTerm(1e-25 * n_Htot**2 * y * sp.log(y), name="underflow")), "quick")
    bad = {f.quantity: f for f in report.nonfinite if f.what == "process 'underflow' row heat"}
    assert "value" in bad and bad["value"].example["x_H+"] < 1e-15


class WrongExp(sp.Function):
    """exp with a wrong derivative"""
    def fdiff(self, argindex=1):
        return 2 * self


def test_wrong_jacobian_is_found(monkeypatch):
    from jaco.model_check import _CSemanticsPrinter
    monkeypatch.setattr(_CSemanticsPrinter, "_print_WrongExp", lambda self, e: f"numpy.exp({self._print(e.args[0])})",
                        raising=False)
    report = jaco.check(toy(ThermalTerm(-1e-25 * n_Htot**2 * WrongExp(-T / 1e4), name="wrong derivative")), "quick")
    assert not report.nonfinite
    found = {(f.what, f.quantity) for f in report.jacobian}
    assert ("process 'wrong derivative' row heat", "d/dT") in found
    assert ("system row T", "Jacobian d/dT") in found


def test_undeclared_symbol_is_flagged():
    report = jaco.check(toy(ThermalTerm(1e-25 * sp.Symbol("G0") * n_Htot, name="typo")), "quick")
    assert report.undeclared == ["G0"] and not report.ok


def test_parameter_without_default_is_reported_unless_given():
    m = toy(ThermalTerm(1e-25 * sp.Symbol("Gamma") * n_Htot, name="heating"), parameters=[Parameter("Gamma")])
    assert jaco.check(m, "quick").missing_values == ["Gamma"]
    assert jaco.check(m, "quick", params={"Gamma": 2.0}).ok


def test_helpers_follow_the_generated_c():
    assert list(_interp1d_const([0.5, 1.5, 2.5, 3.0, 7.0], [1.0, 2.0, 3.0], [10.0, 20.0])) == [10, 10, 20, 20, 20]
    data = np.arange(12.0).reshape(3, 4)
    t = {"t": {"ndim": 2, "shape": (3, 4), "log_axes": [0, 0], "data": data,
               "axis0_min": 0.0, "axis0_max": 2.0, "axis1_min": 0.0, "axis1_max": 3.0}}
    f = _table_function(t)
    assert f((1.5, 2.5), "t", -1) == pytest.approx(4 * 1.5 + 2.5)
    assert f((1.5, 2.5), "t", 0) == pytest.approx(4.0) and f((1.5, 2.5), "t", 1) == pytest.approx(1.0)
    assert f((5.0, 2.5), "t", -1) == pytest.approx(8 + 2.5) and f((5.0, 2.5), "t", 0) == 0.0  # clamped: zero slope
