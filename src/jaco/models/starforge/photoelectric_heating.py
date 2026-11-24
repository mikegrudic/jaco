"""Implementation of photoelectric heating"""

from .symbols import G_0, Z_dust, f_dust, n_Htot, T, psi_grain, log_T
from jaco.math import logistic
from jaco.processes import ThermalTerm

efficiency_eps = 0.049 / (1 + pow(psi_grain / 1925.0, 0.73)) + 0.037 * pow(T / 1.0e4, 0.7) / (1 + psi_grain / 5000.0)
# GIZMO applies this heating only for T < 1e6 K (the efficiency fit grows as T^0.7); smoothed over ~0.05 dex
hot_gas_cutoff = logistic((6.0 - log_T) / 0.01)
photoelectric_heating = ThermalTerm(
    1.3e-24 * G_0 * Z_dust * f_dust * efficiency_eps * n_Htot * hot_gas_cutoff,
    name="Photoelectric Heating",
    bibliography=["1994ApJ...427..822B"],
)
