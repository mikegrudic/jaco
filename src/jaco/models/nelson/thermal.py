"""Thermal terms from Nelson & Langer (1999, 1997) and Hollenbach et al. (1991).

Rates are in erg cm^-3 s^-1. Positive terms heat the gas and negative terms cool it.
"""

import sympy as sp
from jaco.processes import ThermalTerm
from jaco.symbols import T, n_, n_Htot

NL99, NL97 = "1999ApJ...524..923N", "1997ApJ...482..796N"
GENZEL91 = "1991pssf.conf..155G"
GL78, CP87, HTT91 = "1978ApJ...222..881G", "1987A&A...183...98C", "1991ApJ...377..192H"
eV = 1.602176634e-12  # erg

G_0, A_V, tau_uv, dv_dr = (sp.Symbol(s) for s in ("G_0", "A_V", "tau_uv", "dv_dr"))
n_H2 = n_("H_2")
# Combined H I and H_2 collider density, equal to n_H2 in the fully molecular model.
n_atomic = n_Htot - n_H2
n_crit = 3e3  # cm^-3, C+ collisional de-excitation (Genzel 1991)


def log10(x):
    return sp.log(x) / sp.log(10)


def molecular_line_cooling():
    """NL97 eq. 7, the analytic fit to the Goldsmith & Langer (1978) total molecular line cooling"""
    alpha, beta, delta = 5.0e-30, 0.6 + 0.41 * log10(n_H2), 2.3 - 0.18 * log10(n_H2)
    return ThermalTerm(
        -alpha * (T / 10) ** beta * n_H2**delta,
        name="Molecular line cooling",
        bibliography=[NL97, GL78]
    )


def cplus_cooling():
    """NL97 eq. 8"""
    boltz = sp.exp(-92 / T)
    return ThermalTerm(
        -2.2e-23 * n_atomic * n_("C+") * boltz / (1 + n_H2 / n_crit * (1 + 2 * boltz)),
        name="C+ fine-structure cooling",
        bibliography=[NL97, CP87]
    )


def oxygen_cooling():
    """NL97 eq. 9"""
    return ThermalTerm(
        -1e-26 * n_atomic * n_("O") * (24 * sp.exp(-228 / T) + 7 * sp.exp(-326 / T)) * sp.sqrt(T),
        name="O I fine-structure cooling",
        bibliography=[NL97, CP87]
    )


# Genzel (1991) Table 11 gives (wavelength [cm], A [s^-1], critical density [cm^-3]) for J=1-0 and J=2-1.
# The 609 um line uses the H_2 entry. The 370 um line uses the H I entry as a proxy for H_2.
C_I_LINES = ((609.1354e-4, 7.9e-8, 1e4), (370.415e-4, 2.7e-7, 1.2e3))
HC_OVER_K = 1.438777  # h c / k_B in cm K

# Two passes balance cooling accuracy and symbolic Jacobian size. See the README for the convergence comparison.
C_I_PASSES = 2


def _escape_probability(tau):
    """Return the NL99 escape probability, regularized to stay finite and differentiable at tau = 0."""
    tau = tau + 1e-8
    return (1 - sp.exp(-tau)) / tau


def carbon_cooling(passes=C_I_PASSES):
    """Return the 609 and 370 um C I fine-structure cooling from NL99 eq. 7.

        Lambda = n(C I) sum_lines (h nu) A beta f_u,

    Level populations f follow statistical equilibrium with H_2 collisions and radiative rates A beta
    (Genzel 1991, eqs. 22, 26, and 33). The Sobolev optical depth uses the velocity gradient `dv_dr`.

        tau_ul = A_ul lambda^3 / (8 pi dv/dr) n(C I) (f_l g_u / g_l - f_u).

    Each of `passes` iterations solves the populations in closed form at fixed escape probabilities,
    then updates the escape probabilities from the new optical depths.
    """
    (lam10, A10, nc10), (lam21, A21, nc21) = C_I_LINES
    gamma10, gamma21 = A10 / nc10, A21 / nc21  # Collision rate coefficients in cm^3 s^-1.
    E10, E21 = HC_OVER_K / lam10, HC_OVER_K / lam21  # Line energies h nu / k_B in K.

    # Sobolev prefactors for the optical depths.
    K10 = A10 * lam10**3 / (8 * sp.pi * dv_dr)
    K21 = A21 * lam21**3 / (8 * sp.pi * dv_dr)
    n_C = n_("C")

    # P and Q are upward and downward collision rates per C I atom. Detailed balance uses weights g = (1, 3, 5).
    P1, Q1 = gamma10 * n_H2 * 3 * sp.exp(-E10 / T), gamma10 * n_H2  # 0 <-> 1
    P2, Q2 = gamma21 * n_H2 * sp.Rational(5, 3) * sp.exp(-E21 / T), gamma21 * n_H2  # 1 <-> 2

    def chain(beta10, beta21):
        # Only adjacent transitions are included, so equilibrium balances each pair separately.
        # The ratios s = f0/f1 and r2 = f2/f1 give fractions (s, 1, r2) / Y, where Y = s + 1 + r2.
        s = (Q1 + A10 * beta10) / P1
        r2 = P2 / (Q2 + A21 * beta21)
        return s, r2, s + 1 + r2

    # Start with all C I in the ground state. This makes the 370 um line optically thin.
    beta10 = _escape_probability(3 * K10 * n_C)
    beta21 = sp.S.One
    for _ in range(passes):
        s, r2, Y = chain(beta10, beta21)
        beta10 = _escape_probability(K10 * n_C * (3 * s - 1) / Y)
        beta21 = _escape_probability(K21 * n_C * (sp.Rational(5, 3) - r2) / Y)
    s, r2, Y = chain(beta10, beta21)

    # Convert the line energies from K to erg for the per-atom cooling power.
    k_B = 1.380649e-16
    power = k_B * (E10 * A10 * beta10 + E21 * A21 * beta21 * r2) / Y   # per C I atom, erg s^-1
    return ThermalTerm(
        -n_C * power,
        name="C I fine-structure cooling",
        bibliography=[NL99, GENZEL91]
    )


def surface_cooling():
    """NL97 eq. 10, an artificial cooling that keeps the low-density edge of a cloud below ~100 K"""
    return ThermalTerm(
        -7.31e-25 * n_H2 * sp.exp(-200 / T) * sp.exp(-2 * tau_uv),
        name="Surface cooling",
        bibliography=[NL97]
    )


def cosmic_ray_heating():
    """NL97 eq. 11: zeta_p(H_2) dQ_cr n, with zeta = 2.0e-17 s^-1 and dQ_cr = 20 eV (6.4e-28 n)"""
    return ThermalTerm(
        2.0e-17 * 20 * eV * n_H2,
        name="Cosmic-ray heating",
        bibliography=[NL97, GL78]
    )


def h2_formation_heating():
    """NL97 eq. 12: a third of the cosmic-ray heating, times the 2.0 eV of the 4.48 eV H_2 binding energy that goes to
    the gas"""
    return ThermalTerm(
        2.0e-17 * 20 * eV * n_H2 / 3 * 2.0 / 4.48,
        name="H2 formation heating",
        bibliography=[NL97, GL78]
    )


def dust_temperature():
    """HTT91 eq. 5, with nu_0 = 3e15 s^-1, tau_100 = 0.001 and T_0 = 12.2 G_0^0.2 K, floored at NL99's 10 K"""
    nu_0, tau_100 = 3e15, 1e-3
    T_0 = 12.2 * sp.Max(G_0, 1e-30) ** 0.2
    T_d5 = (8.9e-11 * nu_0 * G_0 * sp.exp(-1.8 * A_V) + 2.75**5
            + 3.4e-2 * (0.42 - sp.log(3.45e-2 * tau_100 * T_0)) * tau_100 * T_0**6)
    return sp.Max(10.0, T_d5**0.2)


def gas_dust_exchange():
    """NL97 eq. 13, cooling when the gas is warmer than the dust, with the HTT91 dust temperature"""
    return ThermalTerm(
        -2.4e-33 * sp.sqrt(T) * (T - dust_temperature()) * n_H2**2,
        name="Gas-dust exchange",
        bibliography=[NL97, NL99, GL78, HTT91]
    )


def grain_charging_factor():
    """f(T, n(e), tau_uv) = (1 - x)^2 / x of NL99 eq. 8: HTT91 B7 with the fit x = 1 - 0.28 / (gamma + 0.5) to the cubic
    B8, and gamma = 2.9e-4 Y G_0 T^1/2 exp(-tau_uv) / n(e) (B6), Y = 0.1"""
    gamma = 2.9e-5 * G_0 * sp.sqrt(T) * sp.exp(-tau_uv) / sp.Max(n_("e-"), 1e-30)
    x = 1 - 0.28 / (gamma + 0.5)
    return (1 - x) ** 2 / x


def photoelectric_heating():
    """NL99 eq. 8, 4.86e-26 G_0 n_H f(T, n(e), tau_uv) exp(-tau_uv), per H nucleus as in HTT91 B7"""
    return ThermalTerm(
        4.86e-26 * G_0 * n_Htot * grain_charging_factor() * sp.exp(-tau_uv),
        name="Photoelectric heating",
        bibliography=[NL99, HTT91]
    )


def thermal_terms():
    return [
        molecular_line_cooling(),
        cplus_cooling(),
        oxygen_cooling(),
        carbon_cooling(),
        surface_cooling(),
        cosmic_ray_heating(),
        h2_formation_heating(),
        gas_dust_exchange(),
        photoelectric_heating()
    ]
