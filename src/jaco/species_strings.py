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
from .data.atoms import atoms, atomic_weights, atomic_numbers
from astropy.constants import k_B, m_p, m_e
from astropy import units as u
import sympy as sp

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


def base_species(species: str) -> str:
    """Removes the charge suffix from a species"""
    if species == "e-":
        return "e-"
    if "^" in species:
        return species.split("^")[0]
    else:
        return species.replace("-", "").replace("+", "")


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
                return "^" + str(abs(charge)) + "-"
            else:
                return "^" + str(abs(charge)) + "+"


def remove_electron(species: str) -> str:
    """Returns the symbol of the species produced by removing an electron from the input species"""
    charge = species_charge(species)
    return base_species(species) + charge_suffix(charge + 1)


def add_electron(species: str) -> str:
    """Returns the symbol of the species produced by adding an electron to the input species"""
    charge = species_charge(species)
    return base_species(species) + charge_suffix(charge - 1)


def species_counts(species: str) -> dict:
    """Returns a dict of components (nuclei, electrons) found in a species whose values are the count of that species
    found in it"""

    if species == "e-":  # special behaviour
        return {"e-": 1}

    counts = {}
    formula = base_species(strip(species))  # take off charge suffix
    formula = "".join(formula.split())  # remove whitespace

    # scan through looking for elements, adding up the numbers that follow.
    while formula:
        for num_char in (2, 1, 0):  # find an atomic symbol, starting with the 2-character ones
            if len(formula) < num_char:
                continue
            if formula[:num_char] in atoms:
                atom = formula[:num_char]
                break
        if num_char == 0:
            return {}  # ValueError(f"Could not parse atomic symbols from {species}.")
        formula = formula[num_char:]
        num = 1
        i = 0
        while formula[: i + 1].isnumeric() and i < len(formula):
            num = int(formula[: i + 1])
            i += 1
        formula = formula[i:]
        counts[atom] = num

    charge = species_charge(species)
    num_protons = sum([atomic_numbers[s] * counts[s] for s in counts])
    num_electrons = num_protons - charge
    counts["e-"] = num_electrons
    return counts


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
    mass = 0
    for s, num in species_counts(species).items():
        if s == "e-":
            mass += num * electronmass_cgs
        elif s in atomic_weights:
            mass += num * atomic_weights[s] * (protonmass_cgs + electronmass_cgs)
        else:
            raise ValueError(f"I don't know the mass of species component {s}")
    return mass


def total_atom_abundance(atom: str):
    """Returns symbol representing the abundance of an atom relative to that of H"""
    if atom not in atoms:
        raise (ValueError(f"Cannot define total atomic abundance symbol for {atom}"))
    if atom == "H":
        xtot = 1
    elif atom == "He":
        xtot = sp.Symbol("y")
    else:
        xtot = sp.Symbol("x_" + (atom + ",tot"))
    return xtot


def total_atom_massfrac(atom: str):
    """Returns symbol representing the mass fraction of an atom"""
    if atom not in atoms:
        raise (ValueError(f"Cannot define total atomic mass fraction symbol for {atom}"))
    if atom == "H":
        xtot = sp.Symbol("X")
    elif atom == "He":
        xtot = sp.Symbol("Y")
    else:
        xtot = sp.Symbol("Z_" + atom)
    return xtot


def species_max_abundance(species):
    """Returns symbol representing the maximum possible abundance of a species from atom conservation."""
    if not (isinstance(species, str) or isinstance(species, sp.Symbol)) and hasattr(species, "__iter__"):
        return [species_max_abundance(s) for s in species]
    # if this is an abundance symbol extract the species
    counts = species_counts(str(species))
    if "e-" in counts:
        del counts["e-"]
    print(species, counts)
    return sp.Min(*[total_atom_abundance(a) / counts[a] for a in counts])
