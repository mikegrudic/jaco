"""Tests for manipulation and parsing of string identifiers for chemical species. All conventions should be implemented
in species_strings, and we try a bunch of trickier edge cases here.
"""

from jaco.species_strings import (
    base_species,
    species_charge,
    add_electron,
    remove_electron,
    strip,
    species_counts,
    species_mass,
    species_max_abundance,
    total_atom_abundance,
)
import pytest
import sympy as sp

test_species = [
    "e-",
    "C_6H_12O_6",  # my favorite molecule
    "C_6 H_12 O_63+",
    "C_6 H_12 O_6^3+",
    "C_6H_12O_63+",
    "H_2",
    "H_2+",
    "H2",
    "H2+H-",
    "HCO+",
    "He++",
    "He2+",
    "photon(656nm)",
]


@pytest.mark.parametrize("species", test_species)
def test_species_components(species):
    """Bunch of test cases written by hand to catch edge behaviours"""
    match species:
        case "photon(656nm)":
            assert species_counts(species) == {}
            assert species_mass(species) == 0
        case "e-":
            assert strip(species) == "e-"
            assert species_charge(species) == -1
            assert species_counts(species) == {"e-": 1}
        case "C_6H_12O_6":
            assert base_species(species) == "C_6H_12O_6"
            assert strip(species) == "C6H12O6"
            assert species_charge(species) == 0
            assert remove_electron(species) == species + "+"
            assert add_electron(species) == species + "-"
            assert species_counts(species) == {"C": 6, "H": 12, "O": 6, "e-": 12 + 6 * 6 + 6 * 8}
            assert species_mass(species) == pytest.approx(2.9916e-22)
        case "C_6 H_12 O_63+" | "C_6H_12O_63+":
            assert strip(species) == "C6H12O63+"
            assert species_charge(species) == 1
            assert remove_electron(species) == species + "+"
            assert add_electron(species) == species.replace("+", "")
            assert species_counts(species) == {"C": 6, "H": 12, "O": 63, "e-": 12 + 6 * 6 + 63 * 8 - 1}
        case "C_6 H_12 O_6^3+":
            assert strip(species) == "C6H12O6^3+"
            assert species_charge(species) == 3
            assert remove_electron(species) == "C_6 H_12 O_6^4+"
            assert add_electron(species) == "C_6 H_12 O_6++"
            assert species_counts(species) == {"C": 6, "H": 12, "O": 6, "e-": 12 + 6 * 6 + 6 * 8 - 3}
        case "H_2" | "H2":
            assert strip(species) == "H2"
            assert species_charge(species) == 0
            assert remove_electron(species) == species + "+"
            assert add_electron(species) == species + "-"
            assert species_counts(species) == {"H": 2, "e-": 2}
        case "H-":
            assert strip(species) == "H-"
            assert species_charge(species) == -1
            assert remove_electron(species) == "H"
            assert add_electron(species) == "H--"
            assert species_counts(species) == {"H": 1, "e-": 2}
        case "HCO+":
            assert strip(species) == "HCO+"
            assert species_charge(species) == 1
            assert remove_electron(species) == "HCO++"
            assert add_electron(species) == "HCO"
            assert species_counts(species) == {"H": 1, "C": 1, "O": 1, "e-": 1 + 6 + 8 - 1}
        case "He++":  # HeIII
            assert strip(species) == "He++"
            assert species_charge(species) == 2
            assert remove_electron(species) == "He^3+"  # not actually a thing
            assert add_electron(species) == "He+"
            assert species_counts(species) == {"He": 1, "e-": 0}
        case "He2+":  # this is intended to denote the dihelium ion
            assert strip(species) == "He2+"
            assert species_charge(species) == 1
            assert remove_electron(species) == "He2++"
            assert add_electron(species) == "He2"
            assert species_counts(species) == {"He": 2, "e-": 3}
