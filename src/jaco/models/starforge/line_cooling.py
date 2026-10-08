import sympy as sp
from jaco.processes import collisional_thermal_term
from jaco.symbols import T, T5

# put analytic fits for cooling efficiencies
line_cooling_coeffs = {
    "H": {"e-": 7.5e-19 * sp.exp(-118348 / T) / (1 + sp.sqrt(T5))},  # 1996ApJS..105...19K
    "He+": {"e-": 5.54e-17 * T**-0.397 * sp.exp(-473638 / T) / (1 + sp.sqrt(T5))},  # 1996ApJS..105...19K
}


def LineCoolingSimple(emitter: str, collider=None):
    """Cooling by collisional excitation of emitter by collider, well below the critical density with no ambient
    radiation field: a two-body ThermalTerm (erg cm^-3 s^-1), or a list of them for every known collider if collider is
    None"""
    if emitter not in line_cooling_coeffs:
        raise NotImplementedError(f"Line cooling not implemented for {emitter}")
    if collider is None:
        return [LineCoolingSimple(emitter, c) for c in line_cooling_coeffs[emitter]]
    if collider not in line_cooling_coeffs[emitter]:
        raise NotImplementedError(f"Excitation by collisions with {collider} not implemented for {emitter}")
    return collisional_thermal_term((emitter, collider), -line_cooling_coeffs[emitter][collider],
                                    name=f"{emitter}-{collider} Line Cooling", bibliography=["1996ApJS..105...19K"])
