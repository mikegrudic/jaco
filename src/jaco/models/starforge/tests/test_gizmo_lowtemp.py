"""Low-temperature metal and molecular cooling.

STARFORGE_LEGACY against a transcription of GIZMO's COOL_LOW_TEMPERATURES block (cooling.cc, detailed branch: C+ with
[CI] 609 um, HM79 CO with its LVG cap, H2/HD with mass-fraction collider weights), per nHcgs^2; STARFORGE's [CI] on
neutral carbon."""

import numpy as np
import pytest
import sympy as sp
from jaco.symbols import n_
from ..gizmo_lowtemp import gizmo_carbon_cooling, gizmo_H2_cooling
from ..line_cooling import CI_cooling
from ..symbols import T, n_Htot, x_, G_0, NH, grad_v, X_H, z, x_solar

y, C2 = sp.Symbol("y"), sp.Symbol("C_2")
SOLAR_C = 2.53e-3  # GIZMO FIRE-3 SolarAbundances[2]


def gizmo_lowtemp(Tv, nHcgs, nHp, n_elec, f_molec, G0, column, gradv_kms_pc, Z_C, X_Hfrac, Y_Hefrac):
    """(Lambda_Metals_Neutral, Lambda_H2 + Lambda_HD) per nHcgs^2, before the truncation and CMB factors"""
    EXPmax, T3, sqrt_T, logT = 90.0, Tv / 1e3, np.sqrt(Tv), np.log10(Tv)
    nH0 = 1 - nHp
    ncrit_CO, Sigma_crit_CO = 1.9e4 * sqrt_T, 3.0e-5 * Tv / Z_C
    f_Cplus_CCO = 1. / (1. + (nHcgs / (340. * max(0.1, G0)))**2 / sqrt_T)
    Lambda_Cplus = Z_C * (4.7e-28 * (Tv**0.15 + 1.04e4*n_elec/sqrt_T) * np.exp(-min(91.211/Tv, EXPmax)) + 2.08e-29*np.exp(-min(23.6/Tv, EXPmax)))
    Lambda_CCO = Z_C * Tv*sqrt_T * 2.73e-31 / (1. + (nHcgs/ncrit_CO)*(1. + 1.*max(column, 0.017)/Sigma_crit_CO))
    Lambda_CO_HI = 4.42e-28 * gradv_kms_pc * Tv**4 / (nHcgs * nHcgs)
    Lambda_CO = min(Lambda_CO_HI, (1. - f_Cplus_CCO) * Lambda_CCO)
    Lambda_Metals = f_Cplus_CCO * Lambda_Cplus + Lambda_CO
    Lambda_H2_thick = (6.7e-19*np.exp(-min(5.86/T3, EXPmax)) + 1.6e-18*np.exp(-min(11.7/T3, EXPmax)) + 3.e-24*np.exp(-min(0.51/T3, EXPmax))
                       + 9.5e-22*T3**3.76*np.exp(-min(0.0022/(T3*T3*T3), EXPmax))/(1.+0.12*T3**2.1)) / nHcgs
    Lambda_HD_thin = ((1.555e-25 + 1.272e-26*Tv**0.77)*np.exp(-min(128./Tv, EXPmax)) + (2.406e-25 + 1.232e-26*Tv**0.92)*np.exp(-min(255./Tv, EXPmax))) * np.exp(-min(T3*T3/25., EXPmax))
    q = logT - 3.
    poly = lambda c, v: sum(ci * v**i for i, ci in enumerate(c))
    L = max(nH0 - 2.*f_molec, 0) * X_Hfrac * 10**max(poly([-103., 97.59, -48.05, 10.8, -0.9032], logT), -50.)
    L += Y_Hefrac * 10**max(poly([-23.6892, 2.18924, -0.815204, 0.290363, -0.165962, 0.191914], q), -50.)
    L += f_molec * X_Hfrac * 10**max(poly([-23.9621, 2.09434, -0.771514, 0.436934, -0.149132, -0.0336383], q), -50.)
    L += nHp * X_Hfrac * 10**max(poly([-21.7167, 1.38658, -0.379153, 0.114537, -0.232142, 0.0585389], q), -50.)
    le = poly([-22.1903, 1.5729, -0.213351, 0.961498, -0.910232, 0.137497], q) if logT > 2.30103 else poly([-34.2862, -48.5372, -77.1212, -51.3525, -15.1692, -0.981203], q)
    L += n_elec * X_Hfrac * 10**max(le, -50.)
    f_HD = min(0.00126*f_molec, 4.0e-5*nH0)
    r = L / Lambda_H2_thick
    return nH0 * Lambda_Metals, f_molec * L / (1. + r) + f_HD * Lambda_HD_thin / (1. + (f_HD/f_molec)*r)


# T, nHcgs, x_H+, x_e, x_H2, G0, column (g cm^-2), gradv (km/s/pc)
STATES = [(15.0, 300.0, 1e-9, 1e-6, 0.4, 0.05, 0.03, 2.0), (80.0, 30.0, 1e-6, 1.6e-4, 0.05, 1.0, 3e-3, 1.0),
          (3000.0, 1e4, 1e-4, 1e-4, 0.3, 0.01, 1.0, 10.0), (300.0, 1e6, 1e-10, 1e-8, 0.49, 1e-6, 30.0, 50.0),
          (8.0, 1e2, 1e-9, 1e-6, 0.45, 0.3, 1e-3, 0.3)]


@pytest.mark.parametrize("Tv,nH,xHp,xe,xH2,G0,column,gv", STATES)
def test_legacy_matches_gizmo(Tv, nH, xHp, xe, xH2, G0, column, gv):
    Xv, Yv = 0.7155, 0.2703
    vals = {T: Tv, n_Htot: nH, x_("H+"): xHp, x_("e-"): xe, x_("H_2"): xH2, G_0: G0, NH: column / 1.67262178e-24,
            grad_v: gv * 1e5 / 3.085678e18, X_H: Xv, y: 0.25 * Yv / Xv, z: 0.0, sp.Symbol("x_C,tot"): x_solar("C")}
    metals, molecules = gizmo_lowtemp(Tv, nH, xHp, xe, xH2, G0, column, gv, 1.0, Xv, Yv)
    factor = (Tv - 2.73) / (Tv + 2.73)  # truncation is 1 below 10^4.5 K
    assert -float(gizmo_carbon_cooling.heat.subs(vals)) == pytest.approx(nH**2 * metals * factor, rel=1e-10)
    assert -float(gizmo_H2_cooling.heat.subs(vals)) == pytest.approx(nH**2 * molecules * factor, rel=1e-10)


def test_starforge_CI_on_neutral_carbon():
    """Hocuk+16 rate 2.08e-29 exp(-23.6/T) per neutral H nucleus at solar C, on the carbon not in C+ or CO"""
    Tv, n, xHp = 20.0, 300.0, 1e-6
    xC, xCp, xCO = x_solar("C"), 0.2 * x_solar("C"), 0.5 * x_solar("C")
    vals = {T: Tv, n_Htot: n, x_("H+"): xHp, sp.Symbol("x_C,tot"): xC, x_("C+"): xCp, x_("CO"): xCO, C2: 1.0, z: 0.0}
    expected = 2.08e-29 * np.exp(-23.6 / Tv) * (xC - xCp - xCO) / x_solar("C") * n**2 * (1 - xHp) * (Tv - 2.73) / (Tv + 2.73)
    assert -float(CI_cooling.heat.subs(vals)) == pytest.approx(expected, rel=1e-12)
