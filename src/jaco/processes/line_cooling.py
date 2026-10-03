import sympy as sp
from .thermal_process import ThermalTerm, collisional_thermal_term
from ..symbols import T, T5
from ..data import SolarAbundances

# put analytic fits for cooling efficiencies
line_cooling_coeffs = {
    "H": {"e-": 7.5e-19 * sp.exp(-118348 / T) / (1 + sp.sqrt(T5))},  # 1996ApJS..105...19K
    "He+": {"e-": 5.54e-17 * T**-0.397 * sp.exp(-473638 / T) / (1 + sp.sqrt(T5))},  # 1996ApJS..105...19K
    "C+": {  # 2023MNRAS.519.3154H
        "e-": 1e-27 * 4890 / sp.sqrt(T) * sp.exp(-91.211 / T) / SolarAbundances.x("C"),
        "H": 1e-27 * 0.47 * T**0.15 * sp.exp(-91.211 / T) / SolarAbundances.x("C"),
    },
}
line_cooling_bibliography = {"H": "1996ApJS..105...19K", "He+": "1996ApJS..105...19K", "C+": "2023MNRAS.519.3154H"}


def LineCoolingSimple(emitter: str, collider=None):
    """Cooling by collisional excitation of emitter by collider, well below the critical density with no ambient
    radiation field.

    Parameters
    ----------
    emitter: str
        Emitting excited species
    collider: str, optional
        Exciting colliding species. If None, a list of the terms for every collider with a known rate.

    Returns
    -------
    A two-body ThermalTerm (or a list of them) whose heat is the line cooling rate in erg cm^-3 s^-1
    """
    if emitter not in line_cooling_coeffs:
        raise NotImplementedError(f"Line cooling not implemented for {emitter}")
    if collider is None:
        return [LineCoolingSimple(emitter, c) for c in line_cooling_coeffs[emitter]]
    if collider not in line_cooling_coeffs[emitter]:
        raise NotImplementedError(f"Excitation by collisions with {collider} not implemented for {emitter}")
    return collisional_thermal_term((emitter, collider), -line_cooling_coeffs[emitter][collider],
                                    name=f"{emitter}-{collider} Line Cooling",
                                    bibliography=[line_cooling_bibliography[emitter]])
