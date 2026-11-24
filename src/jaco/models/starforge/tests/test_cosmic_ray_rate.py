"""CR ionization rate against GIZMO's Get_CosmicRayIonizationRate_cgs (no CR fluid, RT_ISRF_BACKGROUND, METALS)."""

import numpy as np
import pytest
import sympy as sp
from jaco.symbols import x_
from ..symbols import cosmicray_ionization_rate_H, NH, ISRF, Z_dust, x_solar

PROTONMASS_CGS = 1.67262178e-24
rate = sp.lambdify((NH, ISRF, Z_dust, x_("Fe")), cosmicray_ionization_rate_H, modules="numpy")


def gizmo_zeta(Sigma_cgs, isrf, Z_over_Zsun, Fe_over_Fesun):
    sigma_0 = 2.23e-3
    u_cr = np.sqrt(isrf) * 1.6e-12
    if Sigma_cgs >= sigma_0:
        u_cr *= np.exp(max(-Sigma_cgs / 100.0, -90.0)) * sigma_0 / Sigma_cgs
    return 1.0e-5 * u_cr + 1.0e-21 * Z_over_Zsun + 1.0e-19 * Fe_over_Fesun


@pytest.mark.parametrize("Sigma", [1e-5, 2e-3, 2.23e-3, 3e-3, 0.1, 10.0, 300.0])
@pytest.mark.parametrize("isrf,Z", [(1.0, 1.0), (10.0, 0.1)])
def test_matches_gizmo(Sigma, isrf, Z):
    expected = gizmo_zeta(Sigma, isrf, Z, Z)
    assert rate(Sigma / PROTONMASS_CGS, isrf, Z, Z * x_solar("Fe")) == pytest.approx(expected, rel=1e-10, abs=0)
