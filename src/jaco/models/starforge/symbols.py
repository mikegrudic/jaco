"""Definition of symbols that are used throughout the model.

NOTE: all quantites are assumed to be in cgs units!
"""

import sympy as sp
from ...symbols import x_, n_
from ...data import SolarAbundances
from ...declarations import Parameter

# Integer mass numbers and solar H mass fraction with which GIZMO's interface packs the x_X parameters
# (x_X = Z_X / (A_X X_H)), so that x_X / x_solar(X) reproduces GIZMO's Z_X / Z_X,sun scalings.
MASS_NUMBER = {"C": 12, "N": 14, "O": 16, "Ne": 20, "Mg": 24, "Si": 28, "S": 32, "Ca": 40, "Fe": 56}
X_H_SOLAR = 1 - SolarAbundances.mass_fraction["Z"] - SolarAbundances.mass_fraction["He"]


def x_solar(element):
    """Solar abundance per H nucleus of an element, in the convention of GIZMO's x_X parameters."""
    return SolarAbundances.mass_fraction[element] / (MASS_NUMBER[element] * X_H_SOLAR)


T = sp.Symbol("T")  # Gas temperature
sqrt_T = sp.sqrt(T)
log_T = sp.log(T, 10.0)
n_Htot = n_("Htot")  # Total number density of H nuclei
X_H = sp.Symbol("X")  # mass fraction of hydrogen
T_dust = sp.Symbol("Td")  # Dust temperature
f_dust = sp.Symbol("f_d")  # Factor accounting for sublimation
Z_dust = sp.Symbol("Z_d")  # Solar-normalized dust abundance. Value of 1 corresponds to Solar neighborhood dust.
G_0 = sp.Symbol("G_0")  # UV radiation field normalized to Habing
G_LW = sp.Symbol("G_LW")  # Lyman-Werner band field in Habing units seen by H2 before its self-shielding
f_shield = sp.Symbol("f_shield")  # Lyman-Werner self-shielding factor
grad_v = sp.Symbol("∇v")  # velocity gradient Frobenius norm in CGS (s^-1)
grad_v_tf = sp.Symbol("∇v_tf")  # Frobenius norm of its trace-free part, |∇v - (∇·v/3) I| (s^-1)
NH = sp.Symbol("N_H")  # column density in nucleons, Sigma/m_p (as GIZMO passes it)
dx = sp.Symbol("Δx")  # effective cell size in cm
ISRF = sp.Symbol("ISRF")  # scaling factor for ISRF and cosmic ray background
EV_CGS = 1.602176634e-12
H2_BINDING_ENERGY = 4.48 * EV_CGS
rho = sp.Symbol("rho")
cs = sp.Symbol("c_s")
T_CMB = sp.Symbol("T_CMB")
A_V = 5.34e-22 * NH * Z_dust * f_dust
z = sp.Symbol("z")  # cosmological redshift
T_cmb = 2.73 * (1 + z)
# GIZMO multiplies molecular, fine-structure and metal-line cooling by this to approximate the CMB bath
cmb_bath_factor = (T - T_cmb) / (T + T_cmb)
# GIZMO's low-temperature cooling block (H2/HD, C+, CO, gas-dust) is cut off where the CIE tables take over
lowtemp_truncation = sp.Piecewise((1, T <= 10**4.5), (sp.exp(-(((log_T - 4.5) / 0.2) ** 2)), True))
# gas-dust coupling additionally falls off where grains are sputtered
dust_sputtering_truncation = sp.Piecewise((1, T <= 3.0e5), (sp.exp(-(((T - 3.0e5) / 2.0e5) ** 2)), True))
# GIZMO's Get_CosmicRayEnergyDensity_cgs: CR energy density falls as Sigma_0/Sigma above Sigma_0 = 2.23e-3 g cm^-2,
# with an exponential cut-off at 100 g cm^-2
PROTONMASS_CGS = 1.67262178e-24
cosmicray_attenuation_fac = sp.Min(1, 2.23e-3 / (PROTONMASS_CGS * NH) * sp.exp(-PROTONMASS_CGS * NH / 100.0))
# CR energy density (erg cm^-3): 1 eV cm^-3 scaled by sqrt(ISRF), as GIZMO assumes under RT_ISRF_BACKGROUND
cosmicray_energy_density = sp.sqrt(ISRF) * 1.6e-12 * cosmicray_attenuation_fac
# 1.6e-17 s^-1 per eV cm^-3 of CRs, plus GIZMO's radioactive-decay floor (K-40 ~ Z; short-lived radionuclides ~ Fe)
cosmicray_ionization_rate_H = 1e-5 * cosmicray_energy_density + 1e-21 * Z_dust + 1e-19 * x_("Fe") / x_solar("Fe")
psi_grain = G_0 * sqrt_T / (0.5 * (1.0e-12 + x_("e-")) * n_Htot)  # grain charging parameter


def clumping_estimator(gradient):
    """GIZMO's clumping estimator (update_explicit_molecular_fraction): <n^2>/<n>^2 = 1 + b^2 M^2 of a lognormal density
    PDF with b = 0.5, M = dv / c_s from a velocity-gradient norm across the cell and the thermal speed of molecular gas,
    c_s = v_th,rms / sqrt(3) with v_th,rms = 0.111 sqrt(T) km/s"""
    return 1 + (0.5 * gradient * dx / (1.11e4 * sqrt_T / sp.sqrt(3))) ** 2


# GIZMO's: on the full gradient norm
clumping_factor = clumping_estimator(grad_v)
# STARFORGE's: on the trace-free part only, so that homologous expansion or compression (e.g. SN ejecta) is not
# counted as sub-grid turbulence
clumping_factor_tracefree = clumping_estimator(grad_v_tf)

# GIZMO's density conventions (STARFORGE_LEGACY): CoolingRate and simple_chemistry.cc take nHcgs = 0.76 rho/m_p as the
# H density; update_explicit_molecular_fraction takes rho/m_p. jaco's n_Htot is the H density, X rho/m_H.
GIZMO_HYDROGEN_MASSFRAC = 0.76
nH_gizmo_cooling = GIZMO_HYDROGEN_MASSFRAC * n_Htot / X_H
nH_gizmo_molecular = n_Htot / X_H

# The inputs of both starforge models besides the core ones (n_Htot, pdv_work, y) and those the species imply
PARAMETERS = [
    Parameter("Δx", "cm", 3e18, "cell size, for the clumping factor and the turbulent line width"),
    Parameter("G_0", "Habing", 1.0, "FUV field, for photoelectric heating, grain charging and the C+ fraction"),
    Parameter("G_LW", "Habing", 1.0, "Lyman-Werner field seen by H2 before its self-shielding"),
    Parameter("ISRF", "", 1.0, "scale of the interstellar radiation field and the cosmic-ray background"),
    Parameter("N_H", "cm^-2", 1e20, "shielding column in nucleons, Sigma/m_p"),
    Parameter("Td", "K", 15.0, "dust temperature"),
    Parameter("X", "", X_H_SOLAR, "H mass fraction"),
    Parameter("Z_d", "", 1.0, "dust abundance relative to the solar neighbourhood"),
    Parameter("f_d", "", 1.0, "dust survival fraction (sublimation)"),
    Parameter("f_metal", "", 1.0, "tabulated metal-line cooling switch, 1 on, 0 off"),
    Parameter("f_neb", "", 0.0, "Kim+23 nebular forbidden-line cooling switch, 1 on, 0 off"),
    Parameter("∇v", "s^-1", 1e-14, "Frobenius norm of the velocity gradient"),
    Parameter("z", "", 0.0, "cosmological redshift"),
]
