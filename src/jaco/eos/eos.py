"""Routines to generate the symbolic expression for internal energy, mass density, pressure, and c_v, given the species
in the network.
"""

from ..data.atoms import atomic_weights
from astropy.constants import k_B, m_p, m_e
from astropy import units as u
from ..species_strings import species_charge, base_species as neutralize, species_mass
import sympy as sp


# TODO: implement!
class EOS:
    def __init__(self, species):
        self.species = species  # TODO: should make sure it's chemical species

    @property
    def density(self):
        return sum([species_mass(self.species) * sp.Symbol(f"n_{self.species}")])

    @property
    def internal_energy(self):
        return

    @property
    def heat_capacity(self):
        return

    @property
    def pressure(self):
        """Pressure via ideal gas law"""
        # TODO - non-maxwellian/relativistic species?
        return sum([sp.Symbol(f"n_{s}") for s in self.species]) * k_B * sp.Symbol("T")
