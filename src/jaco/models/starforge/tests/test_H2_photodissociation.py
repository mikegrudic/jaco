"""H2 photodissociation against GIZMO's update_explicit_molecular_fraction (COOL_MOLECFRAC_NONEQM).

GIZMO evolves the H2 mass fraction per neutral H, f, with dissociation term G_LW * y_ss * f, G_LW = 3.3e-11 urad_G0 / 2;
per molecule that is 3.3e-11 urad_G0 y_ss. Its self-shielding y_ss (Gnedin & Draine 2014) takes the Doppler parameter in
km/s.
"""

import numpy as np
import pytest
import sympy as sp
from jaco.symbols import x_, n_
from ..h2_chemistry.photochemistry import photodissociation, f_selfshield_H2
from ..symbols import T, grad_v, dx, NH, G_LW

PROTONMASS_CGS = 1.67262178e-24


def gizmo_y_ss(Tv, column_cgs, xH0, fH2, gradv_cgs, dx_cm):
    """y_ss as GIZMO computes it, for neutral fraction xH0 and f = 2 n_H2 / n_H0"""
    surface_density_H2_0, x_exp_fac, w0 = 5.0e14 * PROTONMASS_CGS, 0.00085, 0.035
    surface_density_local = xH0 * column_cgs
    v_thermal_rms = 0.111 * np.sqrt(Tv)  # km/s
    dv_turb = gradv_cgs * dx_cm / 1e5  # km/s
    x00 = surface_density_local / surface_density_H2_0
    x01 = x00 / (np.sqrt(1.0 + 3.0 * dv_turb**2 / v_thermal_rms**2) * np.sqrt(2.0) * v_thermal_rms)
    x_ss_1, x_ss_sqrt = 1.0 + fH2 * x01, np.sqrt(1.0 + fH2 * x00)
    return (1.0 - w0) / (x_ss_1 * x_ss_1) + w0 / x_ss_sqrt * np.exp(-min(90.0, x_exp_fac * x_ss_sqrt))


fss = sp.lambdify((T, NH, x_("H_2"), grad_v, dx), f_selfshield_H2(), modules="numpy")
rate = sp.lambdify((T, NH, x_("H_2"), grad_v, dx, G_LW, n_("H_2")), photodissociation("H_2").network["H_2"].rhs,
                   modules="numpy")


@pytest.mark.parametrize("Tv", [10.0, 100.0, 3000.0])
@pytest.mark.parametrize("column", [1e-6, 3e-3, 0.1])
@pytest.mark.parametrize("xH0,xH2", [(1.0, 1e-6), (1.0, 0.3), (0.9, 0.45)])
@pytest.mark.parametrize("gradv", [1e-16, 3e-14])
def test_selfshielding_matches_gizmo(Tv, column, xH0, xH2, gradv):
    dx_cm = 3e18
    expected = gizmo_y_ss(Tv, column, xH0, 2 * xH2 / xH0, gradv, dx_cm)
    assert fss(Tv, column / PROTONMASS_CGS, xH2, gradv, dx_cm) == pytest.approx(expected, rel=1e-10, abs=0)


@pytest.mark.parametrize("G", [1e-3, 1.0, 30.0])
def test_dissociation_rate_matches_gizmo(G):
    Tv, column, xH2, gradv, dx_cm, n_H2 = 20.0, 1e-2, 0.2, 1e-14, 3e18, 50.0
    per_molecule = 3.3e-11 * G * gizmo_y_ss(Tv, column, 1.0, 2 * xH2, gradv, dx_cm)
    assert -rate(Tv, column / PROTONMASS_CGS, xH2, gradv, dx_cm, G, n_H2) == pytest.approx(per_molecule * n_H2, rel=1e-10)


def test_selfshielding_prescription_reaches_the_reaction():
    wg = photodissociation("H_2", self_shielding="Wolcott-Green 2011")
    assert sp.Symbol("N_H_2") in wg.rate.free_symbols
    assert sp.simplify(wg.rate - 3.3e-11 * G_LW * f_selfshield_H2("Wolcott-Green 2011") * n_("H_2")) == 0
