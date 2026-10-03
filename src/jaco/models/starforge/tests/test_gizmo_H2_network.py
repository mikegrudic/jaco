"""STARFORGE_LEGACY's H2 network against a transcription of GIZMO's update_explicit_molecular_fraction (cooling.cc), and
STARFORGE's grain formation rate R n_H,tot n_HI."""

import numpy as np
import pytest
import sympy as sp
from jaco.symbols import n_
from ..h2_chemistry.gizmo_network import GizmoH2Network
from ..h2_chemistry.grain_formation import grain_formation, H2_dust_formation_rate
from ..symbols import (T, x_, n_Htot, T_dust, Z_dust, f_dust, G_LW, NH, grad_v, dx, ISRF, X_H, cosmicray_ionization_rate_H)

C2, C3, y = sp.Symbol("C_2"), sp.Symbol("C_3"), sp.Symbol("y")


def gizmo_dfdt(Tv, nH_cgs, xH0, x_e, nhp, nHe0, nHep, nHepp, fH2, Tdust, f_dustgas, urad_G0, column_neutral, dv_turb_kms,
               zeta, clump):
    """d fH2/dt of update_explicit_molecular_fraction (the implicit solve's g[f]) without the D/D+ channels (x 1e-10)"""
    EXPmax = 90.0
    log_T, ln_T, sqrt_T = np.log10(Tv), np.log(Tv), np.sqrt(Tv)
    XH, MMF = 0.76, xH0 * fH2  # MolecularMassFraction = xH0 * fH2
    xH2_guess = XH * min(max(MMF, 0), 1)
    xH_guess, xHe_guess = max(XH - xH2_guess, 0), nHe0 + nHep + nHepp
    logT4 = log_T - 4
    ncr_H, ncr_H2 = 10 ** (3.0 - 0.416 * logT4 - 0.327 * logT4**2), 10 ** (4.845 - 1.3 * logT4 + 1.62 * logT4**2)
    ncr_He = 10 ** (5.0792 * (1 - 1.23e-5 * (Tv - 2000)))
    ncrit = 1 / (xH_guess / ncr_H + xH2_guess / ncr_H2 + xHe_guess / ncr_He)
    f_v0 = 1 / (1 + nH_cgs / ncrit)
    f_LTE = 1 - f_v0
    lnTf = ln_T  # T >= 100 K in the states below, where GIZMO's fit is positive and jaco's freeze is inactive
    b_H2Hp = max(0., -3.3232183e-7 + 3.3735382e-7*lnTf - 1.4491368e-7*lnTf**2 + 3.4172805e-8*lnTf**3 - 4.7813720e-9*lnTf**4
                 + 3.9731542e-10*lnTf**5 - 1.8171411e-11*lnTf**6 + 3.5311932e-13*lnTf**7) * np.exp(-min(21237.15/Tv, EXPmax)) * (nhp*nH_cgs) * clump
    b_H2e = 10 ** (f_v0 * np.log10(4.49e-9 * Tv**0.11 * np.exp(-min(101858./Tv, EXPmax))) + f_LTE * np.log10(1.91e-9 * Tv**0.136 * np.exp(-min(53407.1/Tv, EXPmax)))) * (x_e*nH_cgs) * clump
    b_H2HI = 10 ** (f_v0 * np.log10(6.67e-12 * sqrt_T * np.exp(-min(1.+63593./Tv, EXPmax))) + f_LTE * np.log10(3.52e-9 * np.exp(-min(43900./Tv, EXPmax)))) * (xH0*nH_cgs) * clump
    b_H2H2 = 10 ** (f_v0 * np.log10(5.996e-30 * Tv**4.1881 * np.exp(-min(54657.4/Tv, EXPmax)) / (1 + 6.761e-6*Tv)**5.6881) + f_LTE * np.log10(1.3e-9 * np.exp(-min(53300./Tv, EXPmax)))) * (xH0*nH_cgs/2.) * clump
    b_H2He = 10 ** (f_v0 * (-27.029 + 3.801*log_T - 29487./Tv) + f_LTE * (-2.729 - 1.75*log_T - 23474./Tv)) * (nHe0*nH_cgs) * clump
    b_H2Hep = (3.7e-14*np.exp(min(35./Tv, EXPmax)) + 7.2e-15) * ((nHep+nHepp)*nH_cgs) * clump
    b_H2ext = (b_H2Hp + b_H2e + b_H2He + b_H2Hep) / 2
    b_H2HI /= 2
    b_H2H2 /= 4
    nH0 = xH0 * nH_cgs
    a_Z = 3.e-18*sqrt_T / ((1. + 4.e-2*np.sqrt(Tv+Tdust) + 2.e-3*Tv + 8.e-6*Tv*Tv) * (1. + 1.e4/np.exp(min(EXPmax, 600./Tdust)))) * f_dustgas * nH0 * clump
    lnTeV, R51_n = ln_T - 9.35915, 3.62e-17 / nH_cgs
    k1 = -17.845 + 0.762*log_T + 0.1523*log_T**2 - 0.03274*log_T**3 if Tv <= 6000 else -16.420 + 0.1998*log_T**2 - 5.447e-3*log_T**4 + 4.0415e-5*log_T**6
    k1 = 10 ** max(k1, -50.)
    k2 = 1.5e-9 if Tv <= 300 else 4.0e-9 * Tv**-0.17
    k5 = 5.7e-6/sqrt_T + 6.3e-8 - 9.2e-11*sqrt_T + 4.4e-13*Tv
    c15 = [-1.801849334e1, 2.36085220e0, -2.827443e-1, 1.62331664e-2, -3.36501203e-2, 1.17832978e-2, -1.65619470e-3, 1.06827520e-4, -2.63128581e-6]
    k15 = np.exp(max(-EXPmax, sum(c * lnTeV**i for i, c in enumerate(c15))))
    c16 = [-2.0372609e1, 1.13944933e0, -1.4210135e-1, 8.4644554e-3, -1.4327641e-3, 2.0122503e-4, 8.6639632e-5, -2.5850097e-5, 2.4555012e-6, -8.0683825e-8]
    k16 = 1.46629e-16 * Tv**1.78186 if Tv < 1160.45 else np.exp(max(-EXPmax, sum(c * lnTeV**i for i, c in enumerate(c16))))
    k17 = 6.9e-9 * Tv**-0.35 if Tv <= 8000 else 9.6e-7 * Tv**-0.90
    x_p = min(max(nhp, x_e/10.), 2.)
    x_Hminus = k1*xH0*x_e / ((k2+k16)*xH0 + (k5+k17)*x_p + k15*x_e + R51_n)
    a_GP = k2 * x_Hminus * nH0 * clump
    b_3B = (6.0e-32/np.sqrt(sqrt_T) + 2.0e-31/sqrt_T) * nH0 * nH0 * xH0 * clump**3
    G_LW_half, xi_half = 3.3e-11 * urad_G0 / 2, zeta / 2
    x00 = column_neutral / (5.e14 * 1.67262178e-24)
    v_th = 0.111 * sqrt_T
    x01 = x00 / (np.sqrt(1. + 3.*dv_turb_kms**2/v_th**2) * np.sqrt(2.) * v_th)
    w0, x_exp_fac = 0.035, 0.00085
    x_ss_1, x_ss_sqrt = 1. + fH2*x01, np.sqrt(1. + fH2*x00)
    y_ss = (1.-w0)/(x_ss_1*x_ss_1) + w0/x_ss_sqrt*np.exp(-min(EXPmax, x_exp_fac*x_ss_sqrt))
    x_a = b_3B + b_H2HI - b_H2H2
    x_b = a_GP + a_Z + 2.*b_3B + b_H2HI + b_H2ext + xi_half + y_ss * G_LW_half
    x_c = a_GP + a_Z + b_3B
    return x_c - x_b * fH2 + x_a * fH2**2


# T, n_Htot, x_H+, x_e (all electrons), x_He+, x_He++, fH2 (per neutral H), G_LW, N_H, grad_v (s^-1), dx (cm), C_2
STATES = [(100.0, 30.0, 1e-6, 2e-4, 1e-10, 1e-14, 0.1, 1.0, 1e21, 3e-14, 3e18, 1.8),
          (20.0, 1e3, 1e-9, 3e-7, 1e-12, 1e-16, 0.8, 1.0, 1e22, 1e-13, 1e18, 4.0),
          (2000.0, 1e5, 1e-4, 2e-4, 1e-8, 1e-12, 0.5, 3.0, 3e23, 1e-12, 1e17, 1.2),
          (8000.0, 3.0, 0.3, 0.31, 1e-3, 1e-6, 1e-3, 1.0, 1e19, 1e-14, 1e19, 1.0)]


@pytest.mark.parametrize("Tv,nH,xHp,xe,xHep,xHepp,fH2,GLW,NHv,gv,dxv,C2v", STATES)
def test_matches_gizmo(Tv, nH, xHp, xe, xHep, xHepp, fH2, GLW, NHv, gv, dxv, C2v):
    Xv, yv, Td, Zd = 0.7155, 0.0943, 15.0, 1.0
    xH0 = 1 - xHp
    vals = {T: Tv, n_Htot: nH, x_("H+"): xHp, x_("e-"): xe, x_("He+"): xHep, x_("He++"): xHepp, x_("H_2"): 0.5 * fH2 * xH0,
            G_LW: GLW, NH: NHv, grad_v: gv, dx: dxv, C2: C2v, C3: C2v**3, y: yv, X_H: Xv, T_dust: Td, Z_dust: Zd, f_dust: 1.0,
            ISRF: 1.0, x_("Fe"): 3.4e-5}
    zeta = float(cosmicray_ionization_rate_H.subs(vals))
    nH_cgs = nH / Xv  # rho / m_p
    dfdt = gizmo_dfdt(Tv, nH_cgs, xH0, xe, xHp, yv - xHep - xHepp, xHep, xHepp, fH2, Td, Zd, GLW,
                      xH0 * NHv * 1.67262178e-24, gv * dxv / 1e5, zeta, C2v)
    model = float((GizmoH2Network().network["H_2"].rhs / n_Htot).subs(vals))
    assert model == pytest.approx(0.5 * xH0 * dfdt, rel=1e-6)


def test_starforge_grain_formation_uses_all_H_nuclei():
    """R n_H,tot n_HI (Hollenbach & McKee 1979; Glover & Jappsen 2007), not R n_HI^2"""
    p = grain_formation()
    assert sp.simplify(p.rate - H2_dust_formation_rate * n_Htot * n_("H") * C2) == 0
