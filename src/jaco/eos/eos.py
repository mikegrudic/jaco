"""Routines to generate the symbolic expression for internal energy, mass density, pressure, and c_v, given the species
in the network.
"""

from ..data.atoms import atomic_weights, atoms
from astropy.constants import k_B, m_p, m_e
from ..symbols import T
from astropy import units as u
from ..species_strings import species_charge, base_species, species_mass
import sympy as sp
from dataclasses import dataclass

k_B_cgs = k_B.to(u.erg / u.K).value


def H2_energy_over_kbT(T):
    """Fit for the average energy of a H_2 molecule in units of k_B T"""
    p1 = ((0.0510229560518835 * sp.log(T) - 1.59102802948492) * sp.log(T) + 17.4609732213753) * sp.log(
        T
    ) - 64.5635781542307
    p2 = ((0.102728060235312 * sp.log(T) - 2.42078687298443) * sp.log(T) + 19.4951634944834) * sp.log(
        T
    ) - 51.5605088374882
    return 1.5 + 1 / (1 + sp.exp(-p1)) + 1 / (1 + sp.exp(-p2))


def species_energy(species):
    """Returns the average thermal energy per particle for a species"""
    if species == "H_2":  # currently the only defined special heat capacity, likely the only one that matters
        fac = H2_energy_over_kbT(T)
    else:
        fac = 1.5

    return fac * k_B_cgs * T


def species_heat_capacity(species):
    """Returns the average thermal energy per particle for a species"""
    if species == "H_2":  # currently the only defined special heat capacity, likely the only one that matters
        fac = sp.diff(H2_energy_over_kbT(T), T)
    else:
        fac = 1.5
    return fac * k_B_cgs


class EOS:
    def __init__(self, species):
        self.species = species  # make sure it's all chemical species here?

    @property
    def density(self):
        return sum([species_mass(s) * sp.Symbol(f"n_{s}") for s in self.species])

    @property
    def energy_density(self):
        return sum([species_mass(s) * sp.Symbol(f"n_{s}") for s in self.species])

    @property
    def internal_energy(self):
        u = sum([sp.Symbol(f"n_{s}") * k_B_cgs * T for s in self.species])

    @property
    def heat_capacity(self):
        return

    @property
    def pressure(self):
        """Pressure via ideal gas law"""
        return sum([sp.Symbol(f"n_{s}") for s in self.species]) * k_B * sp.Symbol("T")
