"""Tabulated metal-line cooling against GIZMO's LambdaMetal (CoolingRate, GALSF_FB_FIRE_STELLAREVOLUTION > 2 branch)."""

import numpy as np
import pytest
import sympy as sp
from jaco.symbols import n_, x_
from jaco.interpolation import TableInterp2D
from .. import metal_line_cooling as mlc
from ..symbols import T, n_Htot, z

pytestmark = pytest.mark.skipif(not mlc._hdf5_path().is_file(), reason="spcool_tables.hdf5 not available")

k = sp.Symbol("k")
x_tot, x_Cplus, x_CO = sp.symbols("x_tot x_Cp x_CO")


def table_replaced(expr):
    """Expression with every table lookup replaced by the symbol k."""
    return expr.subs({a: k for a in expr.atoms(TableInterp2D)})


def gizmo_metal_volumetric(Lambda_tables, Tv, T_cmb=2.73):
    """GIZMO: taper below 100 K, CMB-bath factor only if the summed rate is positive; negative = heating."""
    logT = np.log10(Tv)
    L = Lambda_tables * (np.exp(-min((2.0 - logT) ** 2 / 0.1, 40.0)) if logT < 2 else 1.0)
    if L > 0:
        L *= (Tv - T_cmb) / (Tv + T_cmb)
    return L


@pytest.mark.parametrize("Tv", [5.0, 30.0, 80.0, 100.0, 1e4, 1e6])
@pytest.mark.parametrize("kval", [1e-20, -1e-21])
def test_carbon_uses_total_carbon(Tv, kval):
    """Table carbon cooling scales with total C, independent of how C is split among C, C+ and CO."""
    heat = table_replaced(mlc.metal_line_cooling_process("C", "Carbon_cooling").heat)
    nH, n_e = 10.0, 0.3
    heat = heat.subs({x_("C"): x_tot - x_Cplus - x_CO, x_("C+"): x_Cplus, x_("CO"): x_CO})
    f = sp.lambdify((k, T, n_("e-"), n_Htot, x_tot, x_Cplus, x_CO, sp.Symbol("C_2"), z), heat, modules="numpy")
    expected = gizmo_metal_volumetric(kval * n_e * 3e-4 * nH, Tv)
    for xCp, xCO in [(0.0, 0.0), (3e-4, 0.0), (1e-4, 1.5e-4)]:
        assert -f(kval, Tv, n_e, nH, 3e-4, xCp, xCO, 1.0, 0.0) == pytest.approx(expected, rel=1e-12, abs=0)


def test_oxygen_includes_CO_and_other_elements_unchanged():
    rate_O = table_replaced(mlc.metal_line_cooling_rate("O", "Oxygen_cooling")).subs(T, 1e4)
    assert sp.simplify(rate_O - k * n_("e-") * n_Htot * (x_("O") + x_("CO"))) == 0
    rate_Fe = table_replaced(mlc.metal_line_cooling_rate("Fe", "Iron_cooling")).subs(T, 1e4)
    assert sp.simplify(rate_Fe - k * n_("e-") * n_Htot * x_("Fe")) == 0


def test_cmb_factor_acts_on_net_sum():
    """Net heating from one element offsets cooling from another before GIZMO's CMB-bath factor is applied."""
    process = mlc.metal_line_cooling()
    tables = sorted(process.heat.atoms(TableInterp2D), key=str)
    assert len(tables) == len(mlc.METAL_SPECIES)
    rest = {sp.Symbol("C_2"): 1.0, z: 0.0, n_("e-"): 1.0, n_Htot: 1.0, x_("C+"): 0.0, x_("CO"): 0.0, sp.Symbol("f_metal"): 1.0}
    rest.update({x_(s): 1.0 for _, s in mlc.METAL_SPECIES})
    for Tv, values in [(20.0, [2e-21] + [-1e-21] * 8), (20.0, [-2e-21] + [1e-22] * 8), (500.0, [1e-21] * 9)]:
        heat = process.heat.xreplace({tab: sp.Float(v) for tab, v in zip(tables, values)})
        expected = gizmo_metal_volumetric(sum(values), Tv)
        assert -float(heat.subs({**rest, T: Tv})) == pytest.approx(expected, rel=1e-12, abs=0)


def test_switched_by_f_metal():
    """GIZMO applies the tables only with a UV background loaded; f_metal carries that switch."""
    heat = mlc.metal_line_cooling().heat
    assert heat.subs(sp.Symbol("f_metal"), 0) == 0
    assert sp.simplify(sp.diff(heat, sp.Symbol("f_metal"), 2)) == 0
