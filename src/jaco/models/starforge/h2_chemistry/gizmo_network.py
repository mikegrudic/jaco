"""GIZMO's single-species H2 network (update_explicit_molecular_fraction in cooling.cc), for STARFORGE_LEGACY.

GIZMO evolves f = 2 n_H2 / n_H0, the molecular mass fraction of the neutral H (n_H0 = neutral nuclei, H2 included),
as df/dt = x_c - x_b f + x_a f^2 with per-molecule rates halved (quartered for H2-H2) for the mass-fraction variable,
at the density rho/m_p, with critical-density weights from the fixed mass fraction 0.76 and the cell's H2, and its
own equilibrium H- abundance. Here x_H2 = x_H0 f / 2 at fixed x_H0. The D and D+ dissociation channels, which GIZMO
scales by 1e-10, are left out. C_2 and C_3 are GIZMO's clumping_factor and its cube, the only clumped terms of the
legacy model.
"""

import sympy as sp
from jaco.process import Process
from jaco.symbols import n_
from ..symbols import (T, log_T, x_, n_Htot, T_dust, Z_dust, f_dust, G_LW, cosmicray_ionization_rate_H,
                       nH_gizmo_molecular, GIZMO_HYDROGEN_MASSFRAC)
from .photochemistry import f_selfshield_H2
from .radiative_association import radiative_association
from .associative_detachment import k2
from .mutual_neutralization import k5
from .collisional_detachment import Hminus_collisional_detachment

C2, C3 = sp.Symbol("C_2"), sp.Symbol("C_3")
ln_T = sp.log(T)
n = nH_gizmo_molecular  # GIZMO's nH_cgs = rho / m_p
xH0 = 1 - x_("H+")
x_e = x_("e-")
x_He0 = sp.Symbol("y") - x_("He+") - x_("He++")


def _lte_interpolated(f_v0, ln_k_v0, ln_k_lte):
    """GA08 2.1.3: 10^(f_v0 log10 k_v0 + (1 - f_v0) log10 k_LTE), from natural logarithms"""
    return sp.exp(f_v0 * ln_k_v0 + (1 - f_v0) * ln_k_lte)


def dissociation_rates():
    """(b_H2HI, b_H2H2, b_H2ext): GIZMO's per-f collisional dissociation rates, halved (quartered for H2-H2)"""
    logT4 = log_T - 4
    ncr_H = 10 ** (3.0 - 0.416 * logT4 - 0.327 * logT4**2)
    ncr_H2 = 10 ** (4.845 - 1.3 * logT4 + 1.62 * logT4**2)
    # GIZMO's 10^(5.0792 (1 - 1.23e-5 (T - 2000))), with the exponent capped where it would overflow
    inv_ncr_He = sp.exp(sp.Min(-sp.log(10) * 5.0792 * (1 - 1.23e-5 * (T - 2000)), 500))
    xH2_guess = GIZMO_HYDROGEN_MASSFRAC * 2 * x_("H_2")  # XH * MolecularMassFraction
    xH_guess = sp.Max(GIZMO_HYDROGEN_MASSFRAC - xH2_guess, 0)
    n_ncrit = n * (xH_guess / ncr_H + xH2_guess / ncr_H2 + sp.Symbol("y") * inv_ncr_He)
    f_v0 = 1 / (1 + n_ncrit)

    # Savin+04 (GA08 A1-7), frozen below 100 K where the fit turns negative
    lnTf = sp.log(sp.Max(T, 100.0))
    k_Hp = (-3.3232183e-7 + 3.3735382e-7 * lnTf - 1.4491368e-7 * lnTf**2 + 3.4172805e-8 * lnTf**3
            - 4.7813720e-9 * lnTf**4 + 3.9731542e-10 * lnTf**5 - 1.8171411e-11 * lnTf**6
            + 3.5311932e-13 * lnTf**7) * sp.exp(-21237.15 / T)
    b_H2Hp = k_Hp * x_("H+") * n * C2
    b_H2e = _lte_interpolated(f_v0, sp.log(4.49e-9) + 0.11 * ln_T - 101858.0 / T,
                              sp.log(1.91e-9) + 0.136 * ln_T - 53407.1 / T) * x_e * n * C2
    b_H2HI = _lte_interpolated(f_v0, sp.log(6.67e-12) + 0.5 * ln_T - (1.0 + 63593.0 / T),
                               sp.log(3.52e-9) - 43900.0 / T) * xH0 * n * C2
    b_H2H2 = _lte_interpolated(f_v0, sp.log(5.996e-30) + 4.1881 * ln_T - 54657.4 / T - 5.6881 * sp.log(1 + 6.761e-6 * T),
                               sp.log(1.3e-9) - 53300.0 / T) * (xH0 * n / 2) * C2
    b_H2He = _lte_interpolated(f_v0, sp.log(10) * (-27.029 + 3.801 * log_T - 29487.0 / T),
                               sp.log(10) * (-2.729 - 1.75 * log_T - 23474.0 / T)) * x_He0 * n * C2
    b_H2Hep = (3.7e-14 * sp.exp(35.0 / T) + 7.2e-15) * (x_("He+") + x_("He++")) * n * C2
    return b_H2HI / 2, b_H2H2 / 4, (b_H2Hp + b_H2e + b_H2He + b_H2Hep) / 2


def x_Hminus():
    """GIZMO's equilibrium H- per H nucleus (Glover & Jappsen 2007 terms), with its x_p = clamp(x_H+, x_e/10, 2) and
    photodetachment floor R51 = 3.62e-17 s^-1"""
    k1 = radiative_association("H").rate_coefficient
    k15 = Hminus_collisional_detachment("e-").rate_coefficient
    k16 = Hminus_collisional_detachment("H").rate_coefficient
    k17 = sp.Min(6.9e-9 * T**-0.35, 9.6e-7 * T**-0.9)
    x_p = sp.Min(sp.Max(x_("H+"), x_e / 10), 2)
    return k1 * xH0 * x_e / ((k2 + k16) * xH0 + (k5 + k17) * x_p + k15 * x_e + 3.62e-17 / n)


def formation_rates():
    """(a_Z + a_GP, b_3B / x_H0): dust and H- formation per (1 - f), and three-body formation per (1 - f)^2 over x_H0"""
    R_dust = (3.0e-18 * sp.sqrt(T) / ((1 + 4.0e-2 * sp.sqrt(T + T_dust) + 2.0e-3 * T + 8.0e-6 * T**2)
                                      * (1 + 1.0e4 / sp.exp(sp.Min(90, 600.0 / T_dust)))))
    a_Z = R_dust * Z_dust * f_dust * xH0 * n * C2
    a_GP = k2 * x_Hminus() * xH0 * n * C2
    b_3B_over_xH0 = (6.0e-32 / T**0.25 + 2.0e-31 / sp.sqrt(T)) * (xH0 * n) ** 2 * C3  # GIZMO: nH0^2 xH0 C_3
    return a_Z + a_GP, b_3B_over_xH0


class GizmoH2Network(Process):
    """d n_H2/dt = n_Htot x_H0/2 (x_c - x_b f + x_a f^2), f = 2 x_H2 / x_H0, written without dividing by x_H0"""

    def __init__(self):
        super().__init__(name="GIZMO H2 network", bibliography=["2008MNRAS.388.1627G", "2007ApJ...666....1G",
                                                                 "2013ApJ...773L..25F", "2014ApJ...795...37G"])
        b_H2HI, b_H2H2, b_H2ext = dissociation_rates()
        a_form, b_3B_over_xH0 = formation_rates()
        b_3B = b_3B_over_xH0 * xH0
        G_LW_half, xi_half = 3.3e-11 * G_LW / 2, cosmicray_ionization_rate_H / 2
        x_H2 = x_("H_2")
        x_c = a_form + b_3B
        x_b = a_form + 2 * b_3B + b_H2HI + b_H2ext + xi_half + f_selfshield_H2() * G_LW_half
        x_a_over_xH0 = b_3B_over_xH0 + b_H2HI / xH0 - b_H2H2 / xH0  # both carry a factor x_H0, which cancels
        rate = n_Htot * (xH0 / 2 * x_c - x_b * x_H2 + 2 * x_a_over_xH0 * x_H2**2)
        self.network["H_2"] += rate
        self.network["H"] -= 2 * rate
