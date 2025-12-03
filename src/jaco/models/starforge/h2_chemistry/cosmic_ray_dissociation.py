"""Implementation of dissociation by cosmic rays"""

from ..symbols import cosmicray_ionization_rate_H
from jaco.processes import ChemicalReaction


def cosmic_ray_dissociation(species):
    """Process describing dissociation of molecules by cosmic rays"""
    match species:
        case "H_2":
            rate = cosmicray_ionization_rate_H
            return ChemicalReaction(
                "H_2 -> H + H", rate, name="Dissociation of H_2 by cosmic rays", bibliography=["2016ApJ...831...18C"]
            )
