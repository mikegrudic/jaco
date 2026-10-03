"""Implementation of collisional dissociation of H_2 following 2008MNRAS.388.1627G

Rates are assembled in log space, ln k = f_0 ln k_0 + (1 - f_0) ln k_LTE, from closed-form logarithms of the fits:
the generated code evaluates them in double precision, where k_0 and k_LTE underflow to 0 in cold gas and
d(k_0^f_0)/dT = k_0^f_0 (f_0' ln k_0 + ...) then becomes 0 * inf = nan.
"""

from ..symbols import T, log_T, x_, n_Htot
from .chemical_heat import dissociation_heat
from jaco.processes import Reaction
import sympy as sp

LN10 = sp.log(10.0)


def H2_collisional_dissociation(collider, isotopologue="H_2", chemical_heat=True):
    """Symbolic implementation of the rate coefficient for collisional dissociation of H_2, HD, or D_2

    This implements reactions 9, 10, 11, 108, 109, 110, 112, 113, and 114 from Glover & Abel 2008

    Parameters
    ----------
    n:
        total H number density
    collider: str
        Colliding species (implemented: H+, e-, H, H_2)
    isotopologue: str, optional
        Dissociated isotopologue of H_2 (implemented: H_2, HD, D_2)

    Returns
    -------
    Symbolic expression for collisional dissociation rate coefficient in cgs units
    """

    bib = ["2008MNRAS.388.1627G"]  # cite for the implementation of the LTE interpolation
    # use interpolation function from Glover & Abel 2008 [GA08], section 2.1.3, for interpolating between ground state
    # (v=0) and LTE assumptions for states for collisional dissociation rates
    logT4 = log_T - 4.0
    ln_T = sp.ln(T)
    # inverse critical densities 1/n_cr = 10^-(...). The He one grows as 10^(6.2e-5 T) and overflows double precision
    # above ~5e6 K, so its exponent is capped at 500 (where n/n_cr > 1e217 and f_0 ~ 0 regardless).
    inv_ncr_H = sp.exp(-LN10 * (3.0 - 0.416 * logT4 - 0.327 * logT4 * logT4))
    inv_ncr_H2 = sp.exp(-LN10 * (4.845 - 1.3 * logT4 + 1.62 * logT4 * logT4))
    inv_ncr_He = sp.exp(sp.Min(-LN10 * 5.0792 * (1.0 - 1.23e-5 * (T - 2000.0)), 500.0))
    n_ncrit = n_Htot * (x_("H") * inv_ncr_H + x_("H_2") * inv_ncr_H2 + x_("He") * inv_ncr_He)
    if isotopologue == "HD":
        n_ncrit /= 100  # GA08 2.1.7 - HD gets a 100x higher ncrit
    f_0 = 1.0 / (1.0 + n_ncrit)

    k = None  # set directly when k_0 = k_LTE; otherwise built from ln_k_0 and ln_k_LTE
    match collider:
        case "H+":
            # Savin+04 fit, frozen below 100 K: it turns negative below 76 K, where the rate is < 1e-100 anyway
            ln_T_fit = sp.ln(sp.Max(T, 100.0))
            k = (
                -3.3232183e-7
                + 3.3735382e-7 * ln_T_fit
                - 1.4491368e-7 * ln_T_fit**2
                + 3.4172805e-8 * ln_T_fit**3
                - 4.7813720e-9 * ln_T_fit**4
                + 3.9731542e-10 * ln_T_fit**5
                - 1.8171411e-11 * ln_T_fit**6
                + 3.5311932e-13 * ln_T_fit**7
            ) * sp.exp(-21237.15 / T)
            bib.append("2004ApJ...606L.167S")
        case "e-":
            if isotopologue == "H_2":
                ln_k_0 = sp.log(4.49e-9) + 0.11 * ln_T - 101858.0 / T
                ln_k_LTE = sp.log(1.91e-9) + 0.136 * ln_T - 53407.1 / T
            elif isotopologue == "D_2":
                ln_k_0 = sp.log(8.24e-9) + 0.216 * ln_T - 105388 / T
                ln_k_LTE = sp.log(1.91e-9) + 0.163 * ln_T - 53339.7 / T
            elif isotopologue == "HD":
                ln_k_0 = sp.log(5.09e-9) + 0.128 * ln_T - 103258 / T
                ln_k_LTE = sp.log(1.04e-9) + 0.218 * ln_T - 53070.7 / T
                bib.append("2002PPCF...44.2217T")
            bib.append("2002PPCF...44.1263T")
        case "H":
            ln_k_0 = sp.log(6.67e-12) + 0.5 * ln_T - (1.0 + 63593.0 / T)
            ln_k_LTE = sp.log(3.52e-9) - 43900.0 / T
            bib.append("1983ApJ...270..578L")
            bib.append("1986ApJ...302..585M")
            # can update to Martin 1996 following Glover 2015
        case "H_2":
            ln_k_0 = sp.log(5.996e-30) + 4.1881 * ln_T - 54657.4 / T - 5.6881 * sp.log(1.0 + 6.761e-6 * T)
            ln_k_LTE = sp.log(1.3e-9) - 53300.0 / T
            bib.append("1998ApJ...499..793M")
            bib.append("1987ApJ...318...32S")
        case "He":
            ln_k_0 = LN10 * (-27.029 + 3.801 * log_T - 29487.0 / T)
            ln_k_LTE = LN10 * (-2.729 - 1.75 * log_T - 23474.0 / T)
            bib.append("1987ApJ...318..379D")
        case _:
            raise NotImplementedError(f"Collisional dissociation of {isotopologue} by {collider} not implemented.")

    if k is None:
        # interpolating between v=0 and LTE in log space: k = k_0^f_0 * k_LTE^(1 - f_0)
        k = sp.exp(f_0 * ln_k_0 + (1 - f_0) * ln_k_LTE)
    return Reaction(
        f"H_2 + {collider} -> 2H + {collider}",
        rate_coefficient=k,
        heat_per_reaction=dissociation_heat(chemical_heat),
        name=f"Collisional dissociation of {isotopologue} by {collider}",
        bibliography=bib,
    )

