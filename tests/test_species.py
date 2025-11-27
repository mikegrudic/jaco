"""Tests for manipulation and parsing of string identifiers for chemical species"""

from jaco.species_strings import species_components, charge_suffix, species_charge, recombine, ionize, strip
from jaco.data.atoms import atomic_weights
import pytest

test_species = ["e-", "C_6H_12O_6", "C_6 H_12 O_63+", "C_6H_12O_63+", "H_2", "H_2+", "H2", "H2+H-", "HCO+"]


@pytest.mark.parametrize("species", test_species)
def test_species_components(species):
    # components = species_components(species)
    # charge = species_charge(species)
    # suffix = charge_suffix(species)
    # recombined = recombine(species)
    # ionized = ionize(species)
    # neutral = neutralize(species)
    stripped = strip(species)
    match species:
        case "e-":
            assert stripped == "e-"
        case "C_6H_12O_6":
            assert stripped == "C6H12O6"
        case "C_6 H_12 O_63+" | "C_6H_12O_63+":
            assert stripped == "C6H12O63+"
        case "H_2" | "H2":
            assert stripped == "H2"
        case "H-":
            assert stripped == "H-"
        case "HCO+":
            assert stripped == "HCO+"
