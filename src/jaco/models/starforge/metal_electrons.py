"""Free electrons from metals, as GIZMO's cooling module counts them in neutral gas (COOL_LOW_TEMPERATURES +
SIMPLE_STEADYSTATE_CHEMISTRY): find_abundances_and_rates adds C+, heavy ions, thermally ionized alkalis, O+ and
molecular ions to the H and He electrons. The heavy-ion term is GIZMO's cosmic-ray ionization of the neutral gas,
balanced against recombination and grain capture with the charge assumed to end up on Mg+ by charge transfer (hence
its cap at the Mg abundance); its H/He network has no CR term. GIZMO iterates n_e to convergence because the C+
fraction depends on it; here that fixed point is two Newton steps in ln x_C+ from fully ionized carbon (within GIZMO's
own 1% tolerance on n_e), and its discontinuous heavy-ion regimes are blended (_switch). GIZMO's units: n_H is its
nHcgs (an argument here, n_Htot by default), its MolecularMassFraction is 2 x_H2, its gas density is m_p n_Htot / X
and its column in g cm^-2 is m_p N_H."""

import numpy as np
import sympy as sp
from jaco.math import logistic
from .symbols import T, n_Htot, G_0, NH, Z_dust, f_dust, x_, x_solar, cosmicray_ionization_rate_H, SolarAbundances
from .grain_assisted_recombination import alpha_grain

X_H = sp.Symbol("X")  # hydrogen mass fraction
x_C_tot, x_O_tot = sp.Symbol("x_C,tot"), sp.Symbol("x_O,tot")
x_H2 = x_("H_2")
zeta = cosmicray_ionization_rate_H
PROTONMASS_CGS, ELECTRONMASS_CGS = 1.6726e-24, 9.10953e-28  # GIZMO's constants.h
HYDROGEN_MASSFRAC = 0.76
# electrons on the solved H and He ions
x_e_ions = x_("H+") + x_("He+") + 2 * x_("He++")


def G0_carbon():
    """get_FUV_G0 mode 1: the FUV field further shielded by C and H2 (Gong+2017 Eq. 9)"""
    tau_C = sp.Min(NH * X_H * x_C_tot * 1.6e-17, 100.0)
    r_H2 = sp.Min(2.8e-22 * x_H2 * NH * HYDROGEN_MASSFRAC, 100.0)
    return G_0 * sp.exp(-tau_C) * sp.exp(-r_H2) / (1 + r_H2)


def f_Cplus(x_e, n_H=n_Htot, clumping=1):
    """Fraction of gas-phase C in C+ (f_Cplus): photo- and CR ionization against radiative, dielectronic, grain-assisted
    and H2 recombination, the two-body rates times clumping"""
    ionization = 3.43e-10 * G0_carbon() + 520 * 2 * x_H2 * zeta + 3.85 * zeta
    a, b = sp.sqrt(T / 6.67e-3), sp.sqrt(T / 1.943e6)
    g = 0.7849 + 0.1597 * sp.exp(-49550 / T)
    k_rr = 2.995e-9 / (a * (1 + a) ** (1 - g) * (1 + b) ** (1 + g))
    k_dr = T**-1.5 * (6.346e-9 * sp.exp(-12.17 / T) + 9.793e-9 * sp.exp(-73.8 / T) + 1.634e-6 * sp.exp(-15230 / T))
    k_H2 = 2.31e-13 * T**-1.3 * sp.exp(-23 / T)
    recombination = alpha_grain("C+", x_e, n_H) * n_H + (k_rr + k_dr) * n_H * x_e + k_H2 * n_H * x_H2
    return ionization / (ionization + clumping * recombination)


# gas-phase carbon assumed by return_electron_fraction_from_Cplus (Sofia 2004)
x_C_gas = 1.6e-4 * x_C_tot / x_solar("C")


def Cplus_newton_step(x_e_other, y, n_H=n_Htot):
    """(g, h) of a Newton step in ln x_C+ at x_C+ = y: g = ln f_Cplus(x_e_other + y), h = 1 - y d g/d x_e, so that
    ln x_C+ moves by (g - ln(y / x_C_gas)) / h"""
    xe = sp.Dummy("x_e")
    g = sp.log(f_Cplus(xe, n_H))
    return g.subs(xe, x_e_other + y), 1 - sp.diff(g, xe).subs(xe, x_e_other + y) * y


def _switch(a, b, log_ratio):
    """b where log_ratio > 0 and a where it is < 0, blended over ~0.01 dex so the Newton residual stays continuous
    (GIZMO switches regimes discontinuously, by up to ~35%)"""
    u = log_ratio / (0.01 * np.log(10))
    return logistic(-u) * a + logistic(u) * b  # not a + w (b - a): the unselected branch can be huge


def heavy_ion_electrons():
    """return_electron_fraction_from_heavy_ions (METALS + COOL_METAL_LINES_BY_SPECIES; not the IGM cut of comoving
    runs), with its regime switches blended (_switch)"""
    XH, mu_eff, a = HYDROGEN_MASSFRAC, 2.38, 0.1
    f_dustgas = 0.5 * Z_dust * f_dust * SolarAbundances.mass_fraction["Z"] + 1e-15
    n_ion_max = SolarAbundances.mass_fraction["Mg"] / 24.3 / XH
    m_neutrals, m_grain = mu_eff * PROTONMASS_CGS, 4.189e-12 * 2.4 * a**3
    ngr = (m_neutrals / m_grain) * f_dustgas
    k_ei, y0 = 9.77e-8, (24.305 * PROTONMASS_CGS / ELECTRONMASS_CGS) ** 0.5
    y = np.e * y0
    L = np.log(1 + y)
    psi_0 = 1 - L + L / (1 + L) * np.log(L * (1 + 1 / y))
    k_eg_00 = 0.0195 * a * a * sp.sqrt(T)
    k_eg_0 = k_eg_00 * np.exp(psi_0)
    n_crit = k_ei * zeta / (k_eg_0 * ngr) ** 2
    n_eff = n_Htot * PROTONMASS_CGS / X_H / m_neutrals  # GIZMO's density_cgs / m_neutrals
    alpha = zeta * 16.71 / (a * T) / (k_eg_00 * ngr * ngr * n_eff)
    alpha_min, alpha_max = 0.02, 10.0
    psi_small = 0.5 * (1 - sp.sqrt(1 + 4 * y0 * alpha))
    psi_xmin = 0.5 * (1 - np.sqrt(1 + 4 * y0 * alpha_min))
    psi_mid = psi_0 + (psi_xmin - psi_0) * 2 / (1 + alpha / alpha_min)
    to_per_H = 1 / (XH * mu_eff * n_eff)
    # innermost first: dust-dominated recombination, by grain charge regime
    x = _switch(zeta / (k_eg_00 * sp.exp(psi_small) * ngr) * to_per_H, zeta / (k_eg_00 * sp.exp(psi_mid) * ngr) * to_per_H,
                sp.log(alpha / alpha_min))
    x = _switch(zeta / (k_eg_00 * ngr) * to_per_H, x, sp.log(alpha / 1e-4))
    x = _switch(x, zeta / (k_eg_0 * ngr) * to_per_H, sp.log(alpha / alpha_max))
    # gas-phase recombination at low density, interpolated to the dust regime
    x = _switch(0.5 * (sp.sqrt(4 * n_crit / n_eff + 1) - 1) * (k_eg_0 * ngr) / (k_ei * XH * mu_eff), x,
                sp.log(n_eff / (100 * n_crit)))
    x = _switch(sp.sqrt(zeta / (k_ei * n_eff)) / (XH * mu_eff), x, sp.log(n_eff / (0.01 * n_crit)))
    # highly ionized: negligible
    x = _switch(x, n_ion_max, sp.log(x_e_ions / 0.01))
    return sp.Min(n_ion_max, x)


def alkali_electrons(n_H=n_Htot):
    """return_electron_fraction_from_alkali: Saha-limited K, capped at its abundance"""
    x_K = 1e-7 * Z_dust
    K_over_xe = (x_K * 1.15e-11 * sp.exp(25188 / T)
                 / (6.47e-13 * sp.sqrt(x_K / 1e-7) * (T**3 / 1e9) ** 0.25 * sp.sqrt(2.4e15 / n_H)))
    return sp.Piecewise((0, T < 100), (x_K / (1 + K_over_xe), True))


def Oplus_electrons():
    """return_electron_fraction_from_Oplus: O ionized like H (gas-phase O 3.2e-4, Savage & Sembach 1996)"""
    return 3.2e-4 * x_O_tot / x_solar("O") * x_("H+")


def molecular_ion_electrons(n_H=n_Htot):
    """return_electron_fraction_from_molecular_ions (Fromang+2002 recombination, Armitage 2010); GIZMO's floor
    max(MinGasTemp, T) is left to the solver's temperature floor"""
    return sp.sqrt(zeta / (3e-6 / sp.sqrt(T) * sp.Max(1e2, n_H)))


def metal_electrons(n_H=n_Htot):
    """(intermediates, x_e): free electrons per H nucleus on metals (everything but H and He ions), as an expression of
    the model's intermediates, which are evaluated in order; n_H is GIZMO's nHcgs"""
    e_other, g0, h0, y1, g1, h1, e_Cplus = sp.symbols("eMetalOther CpG0 CpH0 CpY1 CpG1 CpH1 eCplus")
    x_e_other = x_e_ions + e_other
    inter = [(e_other, heavy_ion_electrons() + alkali_electrons(n_H) + Oplus_electrons() + molecular_ion_electrons(n_H))]
    inter += list(zip((g0, h0), Cplus_newton_step(x_e_other, x_C_gas, n_H)))
    inter += [(y1, x_C_gas * sp.exp(g0 / h0))]
    inter += list(zip((g1, h1), Cplus_newton_step(x_e_other, y1, n_H)))
    inter += [(e_Cplus, sp.Min(y1 * sp.exp((g1 - g0 / h0) / h1), x_C_gas))]
    return inter, e_Cplus + e_other
