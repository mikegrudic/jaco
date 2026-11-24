"""Photoelectric heating against GIZMO's CoolingRate() (Bakes & Tielens 1994 / Wolfire 2003 form)."""

import numpy as np
import pytest
import sympy as sp
from jaco.symbols import x_
from ..photoelectric_heating import photoelectric_heating
from ..symbols import T, G_0, Z_dust, f_dust, n_Htot

heat = sp.lambdify((T, G_0, Z_dust, f_dust, n_Htot, x_("e-")), photoelectric_heating.heat, modules="numpy")


def gizmo_pe_volumetric(Tv, G0, Z_over_Zsun, nH, x_e):
    """-LambdaPElec * nH^2, with return_dust_to_metals_ratio_vs_solar = 1 (as under SINGLE_STAR_SINK_DYNAMICS)."""
    if Tv >= 1.0e6:
        return 0.0
    x_pe = G0 * np.sqrt(Tv) / (0.5 * (1.0e-12 + x_e) * nH)
    eps = 0.049 / (1 + (x_pe / 1925.0) ** 0.73) + 0.037 * (Tv / 1.0e4) ** 0.7 / (1 + x_pe / 5000.0)
    return 1.3e-24 * G0 / nH * Z_over_Zsun * eps * nH**2


@pytest.mark.parametrize("Tv", [10.0, 100.0, 8e3, 1e5, 7.5e5])
@pytest.mark.parametrize("G0,nH,x_e", [(1.7, 1.0, 1e-3), (100.0, 1e3, 1e-4), (0.1, 0.1, 0.5)])
def test_matches_gizmo_below_cutoff(Tv, G0, nH, x_e):
    expected = gizmo_pe_volumetric(Tv, G0, 0.5, nH, x_e)
    assert heat(Tv, G0, 0.5, 1.0, nH, x_e) == pytest.approx(expected, rel=1e-4, abs=0)


@pytest.mark.parametrize("Tv", [1.3e6, 1e7, 1e9])
def test_off_in_hot_gas(Tv):
    reference = gizmo_pe_volumetric(9.9e5, 1.7, 1.0, 0.01, 1.2)
    assert abs(heat(Tv, 1.7, 1.0, 1.0, 0.01, 1.2)) < 1e-4 * reference
