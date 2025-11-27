"""
Routines for parsing and manipulating chemical species' string identifiers. All such behaviours should be defined
here - i.e. this is where all of jaco's conventions for writing down chemical formulae should be.

# C_6H_12O_63+ ambiguous: 6 O's with 3+ charge or 63 O's with 1+ charge?
# assume any contiguous number following a _ is the number of the species in the formula, if you want to
# make sure something gets identified as a charge put a ^

- underscores and spaces are valid for readability and latex compatibility, but are ignored in the spec. All should be
compatible with UMIST.

is C3+ (C_3)+ or C^{3+}?
"""

from string import digits
from .data.atoms import atoms, atomic_weights
from astropy.constants import k_B, m_p, m_e
from astropy import units as u

boltzmann_cgs = k_B.to(u.erg / u.K).value
protonmass_cgs = m_p.to(u.g).value
electronmass_cgs = m_e.to(u.g).value


def strip(species: str) -> str:
    """Strip all extraneous characters from the species that will be ignored"""
    stripped = "".join(species.split())  # whitespace
    stripped = stripped.replace("_", "")
    return stripped


def species_charge(species: str) -> int:
    """Returns the charge number of a species from its name"""
    match species[-1]:
        case "-":
            sign = -1
        case "+":
            sign = 1

    if species[-1] not in ("-", "+"):
        return 0
    elif "++" in species:
        return 2
    elif "--" in species:
        return -2
    elif species[-2] not in digits:
        return sign
    elif "^" in species:  # we have a digit possibly specifying the charge, check for a ^
        # get what comes between the ^ and the sign
        charge = int(species.split("^")[-1].split("+")[0].split("-")[0])
        return sign * charge
    else:
        return sign


def is_an_ion(species: str) -> bool:
    return species_charge(species) != 0 and species != "e-"

def species_mass

def base_species(species: str) -> str:
    """Removes the charge suffix from a species"""
    # FIXME: WILL FAIL FOR GENERAL CHEMICAL FORMULAE
    if species == "e-":
        return "e-"
    base = species.rstrip(digits + "-+")
    return base


def neutralize(species: str) -> str:
    """If an ion, return the neutralized version, otherwise just return the species"""
    if species == "e-":
        raise NotImplementedError("neutralization of e- is undefined")
    return base_species(species)


def charge_suffix(charge: int) -> str:
    """Returns the suffix for an input charge number"""
    match charge:
        case 0:
            return ""
        case 1:
            return "+"
        case 2:
            return "++"
        case -1:
            return "-"
        case -2:
            return "--"
        case _:
            if charge < -2:
                return str(abs(charge)) + "-"
            else:
                return str(abs(charge)) + "+"


def ionize(species: str) -> str:
    """Returns the symbol of the species produced by removing an electron from the input species"""
    charge = species_charge(species)
    return base_species(species) + charge_suffix(charge + 1)


def recombine(species: str) -> str:
    """Returns the symbol of the species produced by adding an electron to the input species"""
    charge = species_charge(species)
    return base_species(species) + charge_suffix(charge - 1)


def species_components(species: str) -> dict:
    """Returns a dict of components (nuclei, electrons) found in a species whose values are the number of that species
    found in it"""

    if species == "e-":  # special behaviour
        return {"e-": 1}

    components = {}
    charge = species_charge(species)
    formula = base_species(strip(species))  # take off charge suffix
    formula = "".join(formula.split())  # remove whitespace

    elements = []
    # scan through looking for elements, adding up the numbers that follow. try 2-letter elements first
    for length in (2, 1):
        if formula[:length] in atoms:
            elements.append(formula[:2])
            formula = formula[length:]

    return components


def species_mass(species: str) -> float:
    """
    Returns the mass of a species in g

    Parameters
    ----------
    species: string
        Identifier for the species

    Returns
    -------
    mass: float
        mass of species in g
    """
    # TODO: IMPLEMENT ME
    return 0

#     # proper way to do this: write function that decomposes a general species into nuclei + electron
#     if species in atomic_weights:  # atom
#         return protonmass_cgs * atomic_weights[species]
#     elif neutralize(species) in atomic_weights:  # atomic ion
#         return protonmass_cgs * atomic_weights[neutralize(species)] - electronmass_cgs * species_charge
#     elif "H_2" in
