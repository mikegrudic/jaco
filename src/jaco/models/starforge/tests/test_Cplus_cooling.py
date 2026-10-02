"""C+ fine-structure cooling against GIZMO's Lambda_Cplus (COOL_LOW_TEMPERATURES block of CoolingRate)."""

import numpy as np
import pytest
import sympy as sp
from jaco.symbols import n_
from ..line_cooling import LineCoolingSimple
from ..symbols import T, z

SOLAR = {"Z": 0.0142, "He": 0.2703, "C": 2.53e-3}  # GIZMO FIRE-3 SolarAbundances[0, 1, 2]

process = LineCoolingSimple("C+")
heat = sp.lambdify((T, n_("e-"), n_("H"), n_("C+"), sp.Symbol("C_2"), z), process.heat, modules="numpy")


def gizmo_Cplus_volumetric(Tv, nH, x_e, f_Cplus, Z_C):
    """GIZMO's f_Cplus_CCO * Lambda_Cplus [C+ terms only] * nH0 * nH^2 for fully neutral atomic gas (nH0 = 1),
    with the truncation and CMB-bath factors applied to LambdaMol (z = 0)."""
    T_cmb, logT = 2.73, np.log10(Tv)
    truncation = np.exp(-min(((logT - 4.5) / 0.2) ** 2, 40.0)) if logT > 4.5 else 1.0
    return (nH**2 * f_Cplus * Z_C * 4.7e-28 * (Tv**0.15 + 1.04e4 * x_e / np.sqrt(Tv)) * np.exp(-91.211 / Tv)
            * truncation * (Tv - T_cmb) / (Tv + T_cmb))


@pytest.mark.parametrize("Tv", [5.0, 20.0, 100.0, 1000.0, 8000.0, 4e4, 1e5])
@pytest.mark.parametrize("x_e", [1e-4, 1e-2])
def test_matches_gizmo_at_solar(Tv, x_e):
    nH, f_Cplus = 30.0, 0.7
    X_H = 1 - SOLAR["Z"] - SOLAR["He"]
    x_C_tot = SOLAR["C"] / 12.0 / X_H  # as packed by gizmo_to_jaco
    n_Cplus = f_Cplus * x_C_tot * nH
    expected = gizmo_Cplus_volumetric(Tv, nH, x_e, f_Cplus, 1.0)
    # jaco's e- coefficient is 4890e-27 vs GIZMO's 4.7e-28 * 1.04e4 = 4888e-27
    assert -heat(Tv, x_e * nH, nH, n_Cplus, 1.0, 0.0) == pytest.approx(expected, rel=5e-4, abs=0)
