"""Implementation of dust-gas collisions"""

from .symbols import f_dust, Z_dust, sqrt_T, T, T_dust, n_Htot
import sympy as sp
from jaco.processes import Process
from jaco.equation import Equation
from jaco.symbols import d_dt
import pytest


def gas_dust_heattransfer_coeff(a_grain_angstrom=10.0):
    """Returns expression for the gas-dust heat transfer coefficient

    Parameters
    ----------
    a_grain_angstrom: assumed minimum grain size in angstrom. Note that
    """
    return 1.116e-32 * sqrt_T * (1.0 - 0.8 * sp.exp(-75.0 / T)) * Z_dust * f_dust * (a_grain_angstrom / 10) ** -0.5


# need to make a new kind of process: awkward to do this with NBodyProcess because multiple colliders. Instead,
# we want the rate to go simply as C_2 n_Htot^2
class GasDustCollisions(Process):
    """Process representing gas-dust collisions"""

    def __init__(self):
        super().__init__(name="Gas-dust collisions", bibliography=["1979ApJS...41..555H"])
        self.heat = n_Htot * n_Htot * sp.Symbol("C_2") * gas_dust_heattransfer_coeff() * (T_dust - T)
        self.dust_heat = -self.heat
        self.network["dust heat"] = Equation(d_dt("dust_heat"), self.dust_heat)


gas_dust_collisions = GasDustCollisions()
