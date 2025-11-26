"""Implementation of associative detachment of H and D species"""

from ..symbols import T, H2_formation_heat_cgs
from jaco.processes import ChemicalReaction
import sympy as sp

k2 = sp.Piecewise((1.5e-9, T <= 300.0), (4.0e-9 * T**-0.17, T > 300.0))
# plausible upper and lower_bounds per Glover 2006
k2_upper = 5.0 / 1.3 * k2
k2_lower = 0.5 * k2
bib = [
    "1967ApJ...148L.155S",
    "2008MNRAS.388.1627G",
    "Note: Glover & Abel 2008 consider this rate to be order-of-magnitude uncertain, but take the 1967ApJ...148L.155S value fiducially.",
]  # would want to update to Kreckel 2010 rate 2010Sci...329...69K

k2_new = (
    1.35e-9
    * (sp.Pow(T, 9.8493e-2) + 3.2852e-1 * sp.Pow(T, 5.5610e-1) + 2.771e-7 * sp.Pow(T, 2.1826))
    / (1.0 + 6.191e-3 * sp.Pow(T, 1.0461) + 8.9712e-11 * sp.Pow(T, 3.0424) + 3.2576e-14 * sp.Pow(T, 3.7741))
)  # 2010Sci...329...69K


def associative_detachment(species1, species2):
    species = (species1, species2)
    rate = k2
    product = "H_2"
    if "D-" in species and "H" in species:
        rate *= 0.5
        product = "HD"
    if "D" in species and "H-" in species:
        rate *= 0.5
        product = "HD"
    if "D-" in species and "D" in species:
        product = "D_2"

    return ChemicalReaction(
        f"{species1} + {species2} -> {product} + e-",
        rate,
        heat_per_reaction=H2_formation_heat_cgs,
        name=f"Associative detachment of {species1} with {species2}",
        bibliography=bib,
    )
