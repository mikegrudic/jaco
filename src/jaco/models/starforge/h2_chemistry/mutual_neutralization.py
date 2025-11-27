"""Implementation of mutual neutralization of ions of H and D and their ions"""

from ..symbols import T, sqrt_T
from jaco.processes import ChemicalReaction
from jaco.species_strings import base_species as neutralize
import sympy as sp


bib = [
    "1999MNRAS.304..327C",
    "2008MNRAS.388.1627G",
    "Note: Glover & Abel 2008 consider this rate to be order-of-magnitude uncertain, but take the 1999MNRAS.304..327C value fiducially.",
]
k5 = 5.7e-6 / sqrt_T + 6.3e-8 - 9.2e-11 * sqrt_T + 4.4e-13 * T

r5 = ChemicalReaction(
    "H- + H+ -> H + H",
    rate_coefficient=k5,
    name="Mutual neutralization of H− and H+",
    bibliography=bib,
)

r67 = ChemicalReaction(
    "H- + D+ -> H + D",
    rate_coefficient=k5,
    name="Mutual neutralization of H− and D+",
    bibliography=bib,
)

r68 = ChemicalReaction(
    "H+ + D- -> H + D",
    rate_coefficient=k5,
    name="Mutual neutralization of H+ and D-",
    bibliography=bib,
)

r69 = ChemicalReaction(
    "D+ + D- -> D + D",
    rate_coefficient=k5,
    name="Mutual neutralization of D+ and D-",
    bibliography=bib,
)


def mutual_neutralization(species1, species2):
    name = f"Mutual neutralization of {species1} and {species2}"
    valid_reactants = (("H+", "H-"), ("H-", "D+"), ("H+", "D-"), ("D+", "D-"))
    valid_reactants = (set(s) for s in valid_reactants)
    if set((species1, species2)) not in valid_reactants:
        raise NotImplementedError(f"{name} not implemented.")

    return ChemicalReaction(
        f"{species1} + {species2} -> {neutralize(species1)} + {neutralize(species2)}",
        k5,
        name=name,
        bibliography=bib,
    )
