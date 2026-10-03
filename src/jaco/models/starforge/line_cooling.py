import sympy as sp
from jaco.processes import ThermalTerm, collisional_thermal_term
from jaco.symbols import T, T5, x_, n_Htot
from .symbols import x_solar, cmb_bath_factor, lowtemp_truncation

# put analytic fits for cooling efficiencies
line_cooling_coeffs = {
    "H": {"e-": 7.5e-19 * sp.exp(-118348 / T) / (1 + sp.sqrt(T5))},  # 1996ApJS..105...19K
    "He+": {"e-": 5.54e-17 * T**-0.397 * sp.exp(-473638 / T) / (1 + sp.sqrt(T5))},  # 1996ApJS..105...19K
    "C+": {  # 2023MNRAS.519.3154H: rates per H at solar C, converted to per C+ ion
        "e-": 1e-27 * 4890 / sp.sqrt(T) * sp.exp(-91.211 / T) / x_solar("C") * cmb_bath_factor * lowtemp_truncation,
        "H": 1e-27 * 0.47 * T**0.15 * sp.exp(-91.211 / T) / x_solar("C") * cmb_bath_factor * lowtemp_truncation,
    },
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
    bibliography = ["2023MNRAS.519.3154H"] if emitter == "C+" else ["1996ApJS..105...19K"]
    return collisional_thermal_term((emitter, collider), -line_cooling_coeffs[emitter][collider],
                                    name=f"{emitter}-{collider} Line Cooling", bibliography=bibliography)


# [CI] 609 um fine-structure cooling by collisions with neutral H nuclei (Hocuk+16 rate per H nucleus at solar C, as GIZMO
# carries it), on the neutral carbon the model's C+ and CO leave. GIZMO weights the same rate by the C+ fraction instead.
x_C_neutral = sp.Symbol("x_C,tot") - x_("C+") - x_("CO")
CI_cooling = ThermalTerm(
    -2.08e-29 * sp.exp(-23.6 / T) / x_solar("C") * x_C_neutral * n_Htot**2 * (1 - x_("H+")) * cmb_bath_factor
    * lowtemp_truncation,
    name="[CI] 609 um cooling",
    bibliography=["2016MNRAS.456.2586H"],
    clumping=sp.Symbol("C_2"),
)
