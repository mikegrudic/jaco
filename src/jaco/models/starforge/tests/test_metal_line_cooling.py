"""Structure of the tabulated metal-line cooling against GIZMO's LambdaMetal (CoolingRate, FIRE-3 branch)."""

import numpy as np
import pytest
import sympy as sp
from jaco.symbols import n_, x_
from jaco.interpolation import TableInterp2D
from .. import metal_line_cooling as mlc
from ..symbols import T, n_Htot

pytestmark = pytest.mark.skipif(not mlc._hdf5_path().is_file(), reason="spcool_tables.hdf5 not available")

k = sp.Symbol("k")
x_tot, x_Cplus, x_CO = sp.symbols("x_tot x_Cp x_CO")


def table_replaced(process):
    """Process heat with the table lookup replaced by the symbol k."""
    return process.heat.subs({a: k for a in process.heat.atoms(TableInterp2D)})


def gizmo_metal_volumetric(kval, Tv, n_e, nH, x_element):
    """k * n_e * n_X,tot with GIZMO's sub-100 K taper (tables already normalised per element)."""
    logT = np.log10(Tv)
    taper = np.exp(-min((2.0 - logT) ** 2 / 0.1, 40.0)) if logT < 2 else 1.0
    return kval * taper * n_e * x_element * nH


@pytest.mark.parametrize("Tv", [5.0, 30.0, 80.0, 100.0, 1e4, 1e6])
def test_carbon_uses_total_carbon(Tv):
    """Table carbon cooling scales with total C, independent of how C is split among C, C+ and CO."""
    heat = table_replaced(mlc.metal_line_cooling_process("C", "Carbon_cooling"))
    nH, n_e = 10.0, 0.3
    heat = heat.subs({x_("C"): x_tot - x_Cplus - x_CO, x_("C+"): x_Cplus, x_("CO"): x_CO})
    f = sp.lambdify((k, T, n_("e-"), n_Htot, x_tot, x_Cplus, x_CO, sp.Symbol("C_2")), heat, modules="numpy")
    expected = gizmo_metal_volumetric(1e-20, Tv, n_e, nH, 3e-4)
    for xCp, xCO in [(0.0, 0.0), (3e-4, 0.0), (1e-4, 1.5e-4)]:
        assert -f(1e-20, Tv, n_e, nH, 3e-4, xCp, xCO, 1.0) == pytest.approx(expected, rel=1e-12, abs=0)


def test_oxygen_includes_CO_and_other_elements_unchanged():
    heat_O = table_replaced(mlc.metal_line_cooling_process("O", "Oxygen_cooling")).subs(T, 1e4)
    assert sp.simplify(heat_O + k * n_("e-") * n_Htot * (x_("O") + x_("CO")) * sp.Symbol("C_2")) == 0
    heat_Fe = table_replaced(mlc.metal_line_cooling_process("Fe", "Iron_cooling")).subs(T, 1e4)
    assert sp.simplify(heat_Fe + k * n_("e-") * n_Htot * x_("Fe") * sp.Symbol("C_2")) == 0
