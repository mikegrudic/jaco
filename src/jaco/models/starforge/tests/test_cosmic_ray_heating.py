"""CR heating against GIZMO's CR_gas_heating (no CR fluid, COOL_LOW_TEMPERATURES, RT_ISRF_BACKGROUND, non-comoving)."""

import numpy as np
import pytest
import sympy as sp
from jaco.symbols import x_
from ..cosmic_ray_ionization import cosmic_ray_heating, cosmic_ray_ionization
from ..symbols import NH, ISRF, n_Htot

PROTONMASS_CGS = 1.67262178e-24
heat = sp.lambdify((NH, ISRF, n_Htot, x_("e-"), x_("H+")), cosmic_ray_heating.heat, modules="numpy")


def gizmo_cr_heat_volumetric(Sigma_cgs, isrf, nH, n_elec, nH0):
    """CR_gas_heating(...) * nH^2, with Get_CosmicRayEnergyDensity_cgs's RT_ISRF_BACKGROUND attenuation"""
    sigma_0 = 2.23e-3
    u_cr = np.sqrt(isrf) * 1.6e-12
    if Sigma_cgs >= sigma_0:
        u_cr *= np.exp(max(-Sigma_cgs / 100.0, -90.0)) * sigma_0 / Sigma_cgs
    prefac_CR = u_cr / 1.6e-12
    a_hadronic, f_heat_hadronic = 6.37e-16, 1.0 / 6.0
    b_coulomb_ion_per_GeV = 3.09e-16 * (n_elec + 0.57 * nH0) * 0.76
    per_nH2 = (0.87 * f_heat_hadronic * a_hadronic + 0.53 * b_coulomb_ion_per_GeV) * (1.6e-12 * prefac_CR) / (1.0e-2 + nH)
    return per_nH2 * nH**2


@pytest.mark.parametrize("Sigma", [1e-4, 2.23e-3, 0.05, 300.0])
@pytest.mark.parametrize("nH", [1e-3, 1.0, 1e4])
@pytest.mark.parametrize("x_e,x_Hp", [(1e-7, 1e-8), (1e-3, 9e-4), (1.2, 1.0)])
def test_matches_gizmo(Sigma, nH, x_e, x_Hp):
    expected = gizmo_cr_heat_volumetric(Sigma, 4.0, nH, x_e, 1 - x_Hp)
    assert heat(Sigma / PROTONMASS_CGS, 4.0, nH, x_e, x_Hp) == pytest.approx(expected, rel=1e-12, abs=0)


def test_ionization_carries_no_heat():
    assert cosmic_ray_ionization("H").heat == 0
