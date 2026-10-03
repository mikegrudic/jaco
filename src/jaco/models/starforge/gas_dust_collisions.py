"""Implementation of dust-gas collisions"""

from .symbols import f_dust, Z_dust, sqrt_T, T, T_dust, n_Htot, lowtemp_truncation, dust_sputtering_truncation
import sympy as sp
from jaco.processes import ThermalTerm


def gas_dust_heattransfer_coeff(a_grain_angstrom=10.0):
    """Returns expression for the gas-dust heat transfer coefficient, with GIZMO's high-temperature cut-offs

    Parameters
    ----------
    a_grain_angstrom: assumed minimum grain size in angstrom. Note that
    """
    return (
        1.116e-32 * sqrt_T * (1.0 - 0.8 * sp.exp(-75.0 / T)) * Z_dust * f_dust * (a_grain_angstrom / 10) ** -0.5
        * lowtemp_truncation * dust_sputtering_truncation
    )


def GasDustCollisions():
    """Heat exchange of the gas with the dust in two-body collisions, rate C_2 n_Htot^2 (all gas nuclei with grains in
    proportion to the gas); the dust heat reservoir receives what the gas loses"""
    return ThermalTerm(n_Htot * n_Htot * gas_dust_heattransfer_coeff() * (T_dust - T), name="Gas-dust collisions",
                       bibliography=["1979ApJS...41..555H"], clumping=sp.Symbol("C_2"), reservoir="dust heat")


gas_dust_collisions = GasDustCollisions()
