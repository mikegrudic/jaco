"""Inverse Compton cooling off the CMB, as in GIZMO's evaluate_Compton_heating_cooling_rate().

GIZMO's routine also includes the UV background, a fixed Milky-Way ISRF (or the explicit RT field) and
synchrotron losses; those need radiation and magnetic energy densities the model does not take.
"""

from jaco.processes import ThermalTerm
from jaco.symbols import n_
from .symbols import T, z, T_cmb

# 4 sigma_T k_B / (m_e c) = 2.16e-35 erg s^-1 K^-1 per (eV cm^-3); CMB energy density 0.262 (1+z)^4 eV cm^-3
compton_cooling = ThermalTerm(
    -2.16e-35 * 0.262 * (1 + z) ** 4 * n_("e-") * (T - T_cmb),
    name="Inverse Compton cooling (CMB)",
    bibliography=["1986ApJ...301..522I"],
)
