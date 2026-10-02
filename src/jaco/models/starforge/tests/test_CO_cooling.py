import sympy as sp
import pytest
from jaco.symbols import n_, x_
from ..CO_cooling import CO_cooling, lambda_CO
from ..symbols import T, grad_v, n_Htot, z


@pytest.mark.parametrize("Tval", [10.0, 30.0, 100.0])
@pytest.mark.parametrize("nH", [1e2, 1e4, 1e6])
def test_CO_cooling_removes_energy(Tval, nH):
    x_CO, x_H2 = 1e-4, 0.5
    vals = {T: Tval, n_Htot: nH, grad_v: 3.241e-14, x_("CO"): x_CO, n_("CO"): x_CO * nH, n_("H_2"): x_H2 * nH,
            sp.Symbol("C_2"): 1.0, z: 0.0}
    heat = float(CO_cooling.heat.subs(vals))
    assert heat < 0
    cmb_bath = (Tval - 2.73) / (Tval + 2.73)  # GIZMO's factor on LambdaMol
    assert heat == pytest.approx(-float(lambda_CO.subs(vals)) * x_CO * nH * x_H2 * nH * cmb_bath, rel=1e-12, abs=0)
