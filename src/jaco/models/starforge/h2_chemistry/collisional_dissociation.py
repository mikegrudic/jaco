"""Implementation of collisional dissociation of H_2 following 2008MNRAS.388.1627G"""

from ..symbols import T, log_T, x_, n_Htot, H2_formation_heat_cgs
from jaco.processes import ChemicalReaction
import sympy as sp


def H2_collisional_dissociation(collider, isotopologue="H_2"):
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
    # define some critical densities for different species
    ncr_H = sp.Pow(10.0, 3.0 - 0.416 * logT4 - 0.327 * logT4 * logT4)
    ncr_H2 = sp.Pow(10.0, 4.845 - 1.3 * logT4 + 1.62 * logT4 * logT4)
    ncr_He = sp.Pow(10.0, 5.0792 * (1.0 - 1.23e-5 * (T - 2000.0)))
    ncrit = 1.0 / (x_("H") / ncr_H + x_("H_2") / ncr_H2 + x_("He") / ncr_He)
    if isotopologue == "HD":
        ncrit *= 100  # GA08 2.1.7 - HD gets a 100x higher ncrit
    n_ncrit = n_Htot / ncrit
    f_0 = 1.0 / (1.0 + n_ncrit)

    match collider:
        case "H+":
            k_0 = k_LTE = (  # can we do a polyval type thing?
                -3.3232183e-7
                + 3.3735382e-7 * ln_T
                - 1.4491368e-7 * ln_T**2
                + 3.4172805e-8 * ln_T**3
                - 4.7813720e-9 * ln_T**4
                + 3.9731542e-10 * ln_T**5
                - 1.8171411e-11 * ln_T**6
                + 3.5311932e-13 * ln_T**7
            ) * sp.exp(-21237.15 / T)
            bib.append("2004ApJ...606L.167S")
        case "e-":
            if isotopologue == "H_2":
                k_0 = 4.49e-9 * sp.Pow(T, 0.11) * sp.exp(-101858.0 / T)
                k_LTE = 1.91e-9 * sp.Pow(T, 0.136) * sp.exp(-53407.1 / T)
            elif isotopologue == "D_2":
                k_0 = 8.24e-9 * sp.Pow(T, 0.216) * sp.exp(-105388 / T)
                k_LTE = 1.91e-9 * sp.Pow(T, 0.163) * sp.exp(-53339.7 / T)
            elif isotopologue == "HD":
                k_0 = 5.09e-9 * sp.Pow(T, 0.128) * sp.exp(-103258 / T)
                k_LTE = 1.04e-9 * sp.Pow(T, 0.218) * sp.exp(-53070.7 / T)
                bib.append("2002PPCF...44.2217T")
            bib.append("2002PPCF...44.1263T")
        case "H":
            k_0 = 6.67e-12 * sp.sqrt(T) * sp.exp(-(1.0 + 63593.0 / T))
            k_LTE = 3.52e-9 * sp.exp(-43900.0 / T)
            bib.append("1983ApJ...270..578L")
            bib.append("1986ApJ...302..585M")
            # can update to Martin 1996 following Glover 2015
        case "H_2":
            k_0 = 5.996e-30 * pow(T, 4.1881) * sp.exp(-(54657.4 / T)) / sp.Pow(1.0 + 6.761e-6 * T, 5.6881)
            k_LTE = 1.3e-9 * sp.exp(-53300.0 / T)
            bib.append("1998ApJ...499..793M")
            bib.append("1987ApJ...318...32S")
        case "He":
            k_0 = 10 ** (-27.029 + 3.801 * log_T - 29487.0 / T)
            k_LTE = 10 ** (-2.729 - 1.75 * log_T - 23474.0 / T)
            bib.append("1987ApJ...318..379D")
        case _:
            raise NotImplementedError(f"Collisional dissociation of {isotopologue} by {collider} not implemented.")

    k_0 = sp.Max(k_0, 0)
    k_LTE = sp.Max(k_LTE, 0)
    # interpolating between v=0 and LTE in log space:
    k = k_0**f_0 * k_LTE ** (1 - f_0)
    return ChemicalReaction(
        f"H_2 + {collider} -> 2H + {collider}",
        rate_coefficient=k,
        heat_per_reaction=-H2_formation_heat_cgs,
        name=f"Collisional dissociation of {isotopologue} by {collider}",
        bibliography=bib,
    )


colliders = "H+", "e-", "H_2", "H", "He"
model_process = sum([H2_collisional_dissociation(c) for c in colliders])
