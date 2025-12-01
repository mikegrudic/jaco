"""Implementation of molecular hydrogen chemistry following Hopkins+22"""

import sympy as sp
from . import (
    grain_formation,
    mutual_neutralization,
    three_body,
    associative_detachment,
    radiative_association,
    photochemistry,
    collisional_dissociation,
    collisional_detachment,
    cosmic_ray_dissociation,
)
from jaco.processes import ChemicalReaction
from ..symbols import T

r16 = ChemicalReaction(
    "H+ + H- -> H_2+ + e-", sp.Min(6.9e-9 * T**-0.35, 9.6e-7 * T**-0.9), bibliography=["1978JPhB...11L.671P"]
)  # reaction 16 in Glover & Abel 2008

h2_chemistry_processes = [
    grain_formation.grain_formation,
    three_body.H2_3body_formation,
    collisional_detachment.Hminus_collisional_detachment("H"),
    collisional_detachment.Hminus_collisional_detachment("e-"),
    associative_detachment.associative_detachment("H", "H-"),
    *[collisional_dissociation.H2_collisional_dissociation(c) for c in ("H+", "e-", "H_2", "H", "He")],
    mutual_neutralization.mutual_neutralization("H+", "H-"),
    radiative_association.radiative_association("H"),
    photochemistry.photodissociation("H_2"),
    photochemistry.photodetachment("H-"),
    cosmic_ray_dissociation.cosmic_ray_dissociation("H_2"),
    r16,
]
H2_chemistry = sum(h2_chemistry_processes)
