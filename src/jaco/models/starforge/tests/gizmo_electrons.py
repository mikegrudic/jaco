"""Transcription of GIZMO's neutral-gas free-electron terms (cooling.cc find_abundances_and_rates under
COOL_LOW_TEMPERATURES + SIMPLE_STEADYSTATE_CHEMISTRY; functions in simple_chemistry.cc and cooling.cc), for scalar
inputs in GIZMO's own variables. METALS, COOL_METAL_LINES_BY_SPECIES, RT_ISRF_BACKGROUND, no CR fluid, not comoving,
dust-to-metals ratio 1."""

import numpy as np

PROTONMASS_CGS, ELECTRONMASS_CGS, HYDROGEN_MASSFRAC = 1.6726e-24, 9.10953e-28, 0.76
SOLAR = [0.0142, 0.2703, 2.53e-3, 7.41e-4, 6.13e-3, 1.34e-3, 7.57e-4, 7.12e-4, 3.31e-4, 6.87e-5, 1.38e-3]
WD01 = {"H+": [12.25, 8.074e-6, 1.378, 5.087e2, 1.586e-2, 0.4723, 1.102e-5],
        "C+": [45.58, 6.089e-3, 1.128, 4.331e2, 4.845e-2, 0.8120, 1.333e-4]}


def zeta_cr(column, Z, isrf=1.0):
    """Get_CosmicRayIonizationRate_cgs; column in g cm^-2"""
    u_cr = np.sqrt(isrf) * 1.6e-12
    if column >= 2.23e-3:
        u_cr *= np.exp(max(-column / 100.0, -90.0)) * 2.23e-3 / column
    return 1e-5 * u_cr + 1e-21 * Z[0] / SOLAR[0] + 1e-19 * Z[10] / SOLAR[10]


def get_FUV_G0_mode1(G0_mode0, column, Z, fmol):
    tau_C = min(column * Z[2] / (12 * PROTONMASS_CGS) * 1.6e-17, 100.0)
    r_H2 = min(2.8e-22 * fmol * column * HYDROGEN_MASSFRAC / (2 * PROTONMASS_CGS), 100.0)
    return G0_mode0 * np.exp(-tau_C) * np.exp(-r_H2) / (1 + r_H2)


def alpha_recomb_grain(ion, temp, x_elec, nHcgs, G0, Z):
    psi = G0 * np.sqrt(temp) / (nHcgs * x_elec) + 50
    C = WD01[ion]
    return Z[0] / SOLAR[0] * 1e-14 * C[0] / (1 + C[1] * psi ** C[2] * (1 + C[3] * temp ** C[4] * psi ** (-C[5] - C[6] * np.log(temp))))


def f_Cplus(temp, x_elec, nHcgs, G0, G0_mode1, fmol, zeta, Z):
    ionization_rate = 3.43e-10 * G0_mode1 + 520 * fmol * zeta + 3.85 * zeta
    alpha, beta = np.sqrt(temp / 6.67e-3), np.sqrt(temp / 1.943e6)
    gamma = 0.7849 + 0.1597 * np.exp(-49550 / temp)
    k_rr = 2.995e-9 / (alpha * (1 + alpha) ** (1.0 - gamma) * (1 + beta) ** (1 + gamma))
    k_dr = temp**-1.5 * (6.346e-9 * np.exp(-12.17 / temp) + 9.793e-9 * np.exp(-73.8 / temp) + 1.634e-6 * np.exp(-15230 / temp))
    k_gr = alpha_recomb_grain("C+", temp, x_elec, nHcgs, G0, Z)
    k_cplus_H2 = 2.31e-13 * temp**-1.3 * np.exp(-23 / temp)
    ne, nH2 = nHcgs * x_elec, 0.5 * nHcgs * fmol
    return ionization_rate / (ionization_rate + k_gr * nHcgs + (k_rr + k_dr) * ne + k_cplus_H2 * nH2)


def heavy_ions(temperature, density_cgs, n_elec_HHe, zeta_cr, Z):
    """return_electron_fraction_from_heavy_ions; also returns the regime taken (1-7)"""
    f_dustgas = 0.5 * Z[0] + 1e-15
    n_ion_max = (SOLAR[6] / 24.3) / HYDROGEN_MASSFRAC
    XH = HYDROGEN_MASSFRAC
    if n_elec_HHe > 0.01:
        return n_ion_max, 1
    a_grain_micron, m_ion, mu_eff = 0.1, 24.305 * PROTONMASS_CGS, 2.38
    m_neutrals = mu_eff * PROTONMASS_CGS
    m_grain = 4.189e-12 * 2.4 * a_grain_micron**3
    ngrain_ngas = (m_neutrals / m_grain) * f_dustgas
    k_ei, y0 = 9.77e-8, np.sqrt(m_ion / ELECTRONMASS_CGS)
    y = np.exp(1.0) * y0
    ln_oneplusy = np.log(1.0 + y)
    psi_0 = 1.0 - ln_oneplusy + ln_oneplusy / (1.0 + ln_oneplusy) * np.log(ln_oneplusy * (1.0 + 1.0 / y))
    k_eg_00 = 0.0195 * a_grain_micron**2 * np.sqrt(temperature)
    k_eg_0 = k_eg_00 * np.exp(psi_0)
    n_crit = k_ei * zeta_cr / (k_eg_0 * ngrain_ngas * k_eg_0 * ngrain_ngas)
    n_eff = density_cgs / m_neutrals
    if n_eff < 0.01 * n_crit:
        return min(n_ion_max, np.sqrt(zeta_cr / (k_ei * n_eff)) / (XH * mu_eff)), 2
    if n_eff < 100.0 * n_crit:
        return min(n_ion_max, 0.5 * (np.sqrt(4.0 * n_crit / n_eff + 1.0) - 1.0) * (k_eg_0 * ngrain_ngas) / (k_ei * XH * mu_eff)), 3
    psi_fac = 16.71 / (a_grain_micron * temperature)
    alpha = zeta_cr * psi_fac / (k_eg_00 * ngrain_ngas * ngrain_ngas * n_eff)
    alpha_min, alpha_max = 0.02, 10.0
    if alpha > alpha_max:
        return min(n_ion_max, zeta_cr / (k_eg_0 * ngrain_ngas * XH * mu_eff * n_eff)), 4
    if alpha < 1.0e-4:
        return min(n_ion_max, zeta_cr / (k_eg_00 * ngrain_ngas * XH * mu_eff * n_eff)), 5
    if alpha < alpha_min:
        psi = 0.5 * (1.0 - np.sqrt(1.0 + 4.0 * y0 * alpha))
        return min(n_ion_max, zeta_cr / (k_eg_00 * np.exp(psi) * ngrain_ngas * XH * mu_eff * n_eff)), 6
    psi_xmin = 0.5 * (1.0 - np.sqrt(1.0 + 4.0 * y0 * alpha_min))
    psi = psi_0 + (psi_xmin - psi_0) * 2.0 / (1.0 + alpha / alpha_min)
    return min(n_ion_max, zeta_cr / (k_eg_00 * np.exp(psi) * ngrain_ngas * XH * mu_eff * n_eff)), 7


def alkali(temp, nHcgs, Z):
    if temp < 100:
        return 0.0
    x_K = 1e-7 * Z[0] / SOLAR[0]
    xe = 6.47e-13 * np.sqrt(x_K / 1e-7) * np.sqrt(np.sqrt(temp**3 / 1e9)) * np.sqrt(2.4e15 / nHcgs) * np.exp(-25188 / temp) / 1.15e-11
    return 1.0 / (1 / x_K + 1 / xe)


def Oplus(nHp, Z):
    return Z[4] / SOLAR[4] * 3.2e-4 * nHp


def molecular_ions(temp, nHcgs, zeta, MinGasTemp=0.0):
    return np.sqrt(zeta / (3e-6 / np.sqrt(max(MinGasTemp, temp)) * max(1e2, nHcgs)))


def metal_electrons(temp, nHcgs, density_cgs, n_elec_HHe, nHp, G0, column, fmol, Z, isrf=1.0, iters=4000):
    """Converged x_e - n_elec_HHe of the find_abundances_and_rates loop at fixed ions: (total, C+ part)"""
    zeta = zeta_cr(column, Z, isrf)
    G0_1 = get_FUV_G0_mode1(G0, column, Z, fmol)
    others = heavy_ions(temp, density_cgs, n_elec_HHe, zeta, Z)[0] + alkali(temp, nHcgs, Z) + Oplus(nHp, Z) \
        + molecular_ions(temp, nHcgs, zeta)
    x_C = Z[2] / SOLAR[2] * 1.6e-4
    lo, hi = n_elec_HHe + others, n_elec_HHe + others + x_C  # n_e - base - x_C f(n_e) is increasing: bisect
    for _ in range(iters):
        mid = np.sqrt(lo * hi)
        if mid - n_elec_HHe - others - x_C * f_Cplus(temp, mid, nHcgs, G0, G0_1, fmol, zeta, Z) > 0:
            hi = mid
        else:
            lo = mid
        if hi / lo - 1 < 1e-14:
            break
    ne = np.sqrt(lo * hi)
    Cp = x_C * f_Cplus(temp, ne, nHcgs, G0, G0_1, fmol, zeta, Z)
    return others + Cp, Cp
