"""Implementation of cosmic ray ionization and photodissociation"""

from .symbols import x_, T, cosmicray_ionization_rate_H, n_Htot
from jaco.processes import ChemicalReaction
from jaco.species_strings import remove_electron
from astropy import units as u

# reactions and their rates relative to ionization of H
# Le Teuff 2000
cr_ionization_reaction_rates = {
    "H -> H+ + e-": 1.0,
    "He -> He+ + e-": 1.1,
    "H_2 -> H+ + H + e-": 0.037,
    "H_2 -> H + H": 0.22,
    "H_2 -> H_2+ + e-": 2.0,
    "C -> C+ + e-": 3.8,
    "O -> O+ + e-": 5.7,
    "CO -> CO+ + e-": 6.5,
}


def cosmic_ray_ionization(species):
    """Process describing dissociation of molecules by cosmic rays"""
    equation = f"{species} -> {remove_electron(species)} + e-"
    rate = cr_ionization_reaction_rates[equation] * cosmicray_ionization_rate_H
    if species == "H":
        heat = 20 * u.eV.to(u.erg)
    else:
        heat = 0  # already accounted for heat with H

    return ChemicalReaction(
        equation,
        rate,
        heat_per_reaction=heat,
        name=f"Direct ionization of {species} by cosmic rays",
        bibliography=["2000A&AS..146..157L"],
    )


cr_photoionization_reaction_rates = {
    "C -> C+ + e-": (2800, "1987ApJ...323L.137G"),  # Gredel, Lepp, & Dalgarno 1987
    "CH -> C + H": (4000, "1989ApJ...347..289G"),  # Gredel 1989
    "CH+ -> C+ + H": (960, "1989ApJ...347..289G"),  # Gredel 1989
    "CH_2 -> CH_2+ + e-": (2700, "2000A&AS..146..157L"),  # Le Teuff 2000
    "CH_2 -> CH + H": (2700, "2000A&AS..146..157L"),  # Le Teuff 2000
    "C_2 -> C + C": (130, "1989ApJ...347..289G"),  # Gredel 1989
    "OH -> O + H": (280, "1989ApJ...347..289G"),  # Gredel 1989
    "H_2O -> OH + H": (530, "1989ApJ...347..289G"),  # Gredel 1989
    "O_2 -> O + O": (410, "1989ApJ...347..289G"),  # Gredel 1989
    "O_2 -> O_2+ + e-": (640, "1989ApJ...347..289G"),  # Gredel 1989
    "CO -> C + O": (
        0.21 * T**0.5 * x_("H_2") * x_("CO") ** -0.5,
        "1996ApJ...466..561M",
    ),  # Maloney, Hollenbach & Tielens (1996)
}


def cosmic_ray_photoionization(species):
    """Process describing photoionization/dissociation of a species by cosmic ray-induced photons"""
    equation = f"{species} -> {remove_electron(species)} + e-"
    rate, bib = cr_photoionization_reaction_rates[equation]
    rate *= cosmicray_ionization_rate_H

    return ChemicalReaction(
        equation,
        rate,
        name=f"Cosmic ray photoionization of {species}",
        bibliography=[bib],
    )
