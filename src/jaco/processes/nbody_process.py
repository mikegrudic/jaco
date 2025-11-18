"""Class specifying an N-body collisional process with generic methods"""

from ..process import Process
from ..symbols import n_
import sympy as sp


class NBodyProcess(Process):
    """Process implementing special methods specific to n-body processes, whose rates
    all follow the pattern

    rate per volume = k * prod_i(n_i) for i in species

    rate and heat are promoted from attributes to properties implemented to compute this pattern.

    Parameters
    ----------
    colliding_species:
        iterable of strings representing the colliding species.
    rate_coefficient: sympy.core.symbol.Symbol, optional
        Symbol symbol expression for k
    heat_rate_coefficient: sympy.core.symbol.Symbol, optional
        Symbol symbol expression for the heat rate coefficient = average radiated energy * k
    name: str, optional
        Name of the process
    """

    def __init__(self, colliding_species, rate_coefficient=0, heat_rate_coefficient=0, name: str = ""):
        self.name = name
        self.initialize_network()
        self.colliding_species = colliding_species
        self.rate_coefficient = rate_coefficient
        self.heat_rate_coefficient = heat_rate_coefficient
        self.subprocesses = [self]

    @property
    def heat_rate_coefficient(self):
        """Returns the heat rate coefficient of the N-body process"""
        return self.__heat_rate_coefficient

    @heat_rate_coefficient.setter
    def heat_rate_coefficient(self, value):
        """Ensures that the network is always updated when we update the rate coefficient"""
        self.__heat_rate_coefficient = value
        self.heat = value * sp.prod([n_(c) for c in self.colliding_species])

    @property
    def rate(self):
        """Returns the number of events per unit time and volume"""
        if self.rate_coefficient is None:
            return None
        return self.rate_coefficient * sp.prod([n_(c) for c in self.colliding_species])

    # @property
    # def heat(self):
    #     """Returns the energy radiated per unit time and volume"""
    #     return super().heat

    # #        self.__heat = self.heat_rate_coefficient * sp.prod([sp.Symbol("n_" + c) for c in self.colliding_species])
    # #        return self.__heat

    # @heat.setter
    # def heat(self, value):
    #     raise NotImplementedError(
    #         "Cannot directly set the heat value of an N-body process - set the rate coefficient instead."
    #     )

    @property
    def num_colliding_species(self):
        """Number of colliding species"""
        return len(self.colliding_species)
