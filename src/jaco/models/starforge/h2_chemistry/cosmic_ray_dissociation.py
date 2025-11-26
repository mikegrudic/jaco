"""Implementation of dissociation by cosmic rays"""

from ..symbols import sp, NH, ISRF
from jaco.processes import ChemicalReaction

cosmicray_attenuation_fac = sp.Min(1, 1e21 / NH * sp.exp(-NH / 1e24))
cosmicray_ionization_rate = sp.sqrt(ISRF) * 1.6e-12 * 1e-5 * cosmicray_attenuation_fac


def cosmic_ray_dissociation(species):
    """Process describing dissociation of molecules by cosmic rays"""
    match species:
        case "H_2":
            rate = cosmicray_ionization_rate
            return ChemicalReaction(
                "H_2 -> H + H", rate, name="Dissociation of H_2 by cosmic rays", bibliography=["2016ApJ...831...18C"]
            )
