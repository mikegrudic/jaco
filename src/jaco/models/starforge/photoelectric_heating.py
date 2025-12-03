"""Implementation of photoelectric heating"""

from .symbols import G_0, Z_dust, f_dust, n_Htot, T, psi_grain
from jaco.processes import ThermalProcess

efficiency_eps = 0.049 / (1 + pow(psi_grain / 1925.0, 0.73)) + 0.037 * pow(T / 1.0e4, 0.7) / (1 + psi_grain / 5000.0)
photoelectric_heating = ThermalProcess(
    1.3e-24 * G_0 * Z_dust * f_dust * efficiency_eps * n_Htot,
    name="Photoelectric Heating",
    bibliography=["1994ApJ...427..822B"],
)
