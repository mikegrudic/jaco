import pytest
import sympy as sp
from jaco.symbols import n_
from ..compton import compton_cooling
from ..symbols import T, z

heat = sp.lambdify((T, n_("e-"), z), compton_cooling.heat, modules="numpy")


@pytest.mark.parametrize("Tv,zv", [(2.0, 0.0), (1e4, 0.0), (1e7, 0.0), (1e6, 3.0)])
def test_matches_gizmo_cmb_term(Tv, zv):
    """GIZMO: compton_prefac_eV * n_elec * e_CMB_eV * (T - T_cmb), times nH^2."""
    n_e = 0.1
    expected = 2.16e-35 * n_e * 0.262 * (1 + zv) ** 4 * (Tv - 2.73 * (1 + zv))
    assert -heat(Tv, n_e, zv) == pytest.approx(expected, rel=1e-12, abs=0)
