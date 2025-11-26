"""Implementation of the 3-body H_2 formation channel"""

from jaco.processes import ChemicalReaction
from ..symbols import sp, T, T_dust, f_dust, Z_dust, H2_formation_heat_cgs


def threebody_rate(number: int):
    """Returns the different numbered 3-body rates as implemented in Grackle"""
    match number:
        case 0:
            k = 1e-32 * sp.Min((T / 300) ** -0.38, (T / 300) ** -1.0)
            return k, "Abel 2002"
        case 1:
            # PSS83 three-body rate
            return 5.5e-29 / T, "PSS83"
        case 2:
            # CW83 three-body rate
            return 8.8e-33, "CW83"
        case 3:
            # FH07 three-body rate
            return 1.44e-26 / pow(T, 1.54), "FH07"
        case 4:
            # Rate from Glover (2008), derived by detailed balance from MSM96
            return 7.7e-31 / pow(T, 0.464), "Glover 2008"
        case 5:
            # Rate from Forrey (2013)
            return 6.0e-32 * T**-0.25 + 2.0e-31 * T**-0.5, "2013ApJ...773L..25F"
        case _:
            raise NotImplementedError(f"3-body rate not implemented for n={number}")


rate, bibcode = threebody_rate(5)  # take Forrey 2013 by default

H2_3body_formation = ChemicalReaction(
    "H + H + H -> H_2 + H",
    rate,  # Forrey 2013 rate
    heat_per_reaction=H2_formation_heat_cgs,
    name="3-body formation of H_2",
    bibliography=[bibcode],
)
