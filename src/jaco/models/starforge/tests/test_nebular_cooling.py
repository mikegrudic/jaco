"""Nebular cooling against a direct transcription of GIZMO's CoolingRate() (RT_CHEM_PHOTOION && METALS block)."""

import numpy as np
import pytest
import sympy as sp
from jaco.symbols import n_
from ..nebular_cooling import nebular_cooling, f_neb
from ..symbols import T, Z_dust

C_2 = sp.Symbol("C_2")


def gizmo_nebular_volumetric(T, n_e, n_Hp, Z_over_Zsun):
    """GIZMO's LambdaNeb * nH^2 (erg cm^-3 s^-1); n_e and n_Hp in cm^-3."""
    if not (2.0e3 < T < 5.0e4 and n_e > 0 and n_Hp > 0):
        return 0.0
    T4 = T * 1.0e-4
    lnT4 = np.log(T4)
    ne_2 = n_e / 100.0
    log10_fneb = 0.692 + lnT4 * (-0.586 + lnT4 * (0.816 + lnT4 * (-0.505 + lnT4 * (0.118 + lnT4 * (0.00766 - 0.00508 * lnT4)))))
    rate = (n_e * n_Hp * Z_over_Zsun * 3.68e-23 * np.exp(-min(3.86 / T4, 100.0)) / np.sqrt(T4) * 10**log10_fneb
            / (1.0 + 0.12 * ne_2 ** (0.38 - 0.12 * lnT4)))
    return rate / (1.0 + np.exp(min(10.0 * (T - 2.75e4) / 1.5e4, 60.0)))


heat = sp.lambdify((T, n_("e-"), n_("H+"), Z_dust, f_neb, C_2), nebular_cooling.heat, modules="numpy")


@pytest.mark.parametrize("Tval", [3e3, 5e3, 8e3, 1e4, 1.5e4, 2.5e4, 3e4, 4e4, 4.9e4])
@pytest.mark.parametrize("n_e", [1e-3, 1.0, 1e2, 1e4])
@pytest.mark.parametrize("Z", [0.1, 1.0, 3.0])
def test_matches_gizmo_inside_window(Tval, n_e, Z):
    n_Hp = 0.9 * n_e
    expected = gizmo_nebular_volumetric(Tval, n_e, n_Hp, Z)
    assert -heat(Tval, n_e, n_Hp, Z, 1.0, 1.0) == pytest.approx(expected, rel=1e-7, abs=0)


@pytest.mark.parametrize("Tval", [10.0, 1e3, 1.5e3, 5.5e4, 1e5, 1e7])
def test_negligible_outside_window(Tval):
    """GIZMO switches the term off outside 2e3-5e4 K; here it is smoothly ~0 there (edge smoothed within ~25% of 2e3 K)."""
    peak = gizmo_nebular_volumetric(2.0e4, 1.0, 1.0, 1.0)
    assert abs(heat(Tval, 1.0, 1.0, 1.0, 1.0, 1.0)) < 1e-6 * peak


def test_switch_and_finite_derivatives():
    assert heat(1e4, 1.0, 1.0, 1.0, 0.0, 1.0) == 0.0
    for var in (T, n_("e-")):
        d = sp.lambdify((T, n_("e-"), n_("H+"), Z_dust, f_neb, C_2), sp.diff(nebular_cooling.heat, var), modules="numpy")
        for Tval in (1e2, 2e3, 1e4, 2.75e4, 5e4, 1e6):
            for n_e in (0.0, 1e-10, 1.0):
                assert np.isfinite(d(Tval, n_e, 1.0, 1.0, 1.0, 1.0))
