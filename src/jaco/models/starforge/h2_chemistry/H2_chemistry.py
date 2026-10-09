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
    detailed_balance,
)
from jaco.processes import Reaction
from ..symbols import T
from .chemical_heat import formation_radiation

r16 = Reaction(
    "H+ + H- -> H_2+ + e-", sp.Min(6.9e-9 * T**-0.35, 9.6e-7 * T**-0.9), bibliography=["1978JPhB...11L.671P"]
)  # reaction 16 in Glover & Abel 2008


def h2_chemistry_processes(chemical_heat=True, three_body_detailed_balance=False):
    """The H2/H- reactions; chemical_heat gives the gas the binding energy of the H2 formed and takes that of the H2
    dissociated (h2_formation_radiation takes back what formation radiates); three_body_detailed_balance replaces
    Forrey's three-body rate by detailed_balance.py's"""
    three_body_formation = (detailed_balance.H2_3body_formation_detailed_balance(chemical_heat)
                            if three_body_detailed_balance else [three_body.H2_3body_formation(chemical_heat)])
    return [
        grain_formation.grain_formation(chemical_heat),
        *three_body_formation,
        collisional_detachment.Hminus_collisional_detachment("H"),
        collisional_detachment.Hminus_collisional_detachment("e-"),
        associative_detachment.associative_detachment("H", "H-", chemical_heat),
        *[collisional_dissociation.H2_collisional_dissociation(c, chemical_heat=chemical_heat)
          for c in ("H+", "e-", "H_2", "H", "He")],
        mutual_neutralization.mutual_neutralization("H+", "H-"),
        radiative_association.radiative_association("H"),
        photochemistry.photodissociation("H_2"),
        photochemistry.photodetachment("H-"),
        cosmic_ray_dissociation.cosmic_ray_dissociation("H_2"),
        r16,
    ]


def formation_channel(reaction):
    """chemical_heat's channel ("grain", "H-", "3-body") of an H2 formation reaction, else None"""
    if reaction.name == "Formation of H_2 on dust grains":
        return "grain"
    if reaction.name == "Associative detachment of H with H-":
        return "H-"
    return "3-body" if reaction.name.startswith("3-body formation of H_2") else None


def h2_formation_radiation(chemistry):
    """The cooling terms of the energy radiated by the formation reactions in chemistry that heat the gas
    (chemical_heat.formation_radiation): separate processes, so that an optically-thick limit can cap them and not
    the reactions"""
    return [formation_radiation(p, formation_channel(p)) for p in chemistry if formation_channel(p) and p.heat != 0]
