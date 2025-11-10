"""Implementation of base Process class with methods for managing and solving systems of equations"""

import sympy as sp
import jax.numpy as jnp
import numpy as np
from .numerics import newton_rootsolve
from .symbols import n_, d_dt
from .misc import is_an_ion
from .data import SolarAbundances
from . import eos
from astropy import units
from .equations import EquationSystem, Equation


class Process:
    """
    Top-level class containing a description of a microscopic process

    Most importantly, this implements the procedure for combining processes to build up a network for chemistry
    + conservation equations.
    """

    def __init__(self, name="", bibliography={}):
        """Construct an empty Process instance

        Parameters
        ----------
        name: str, optional
            Name of the process
        """
        self.name = name
        self.initialize_network()
        self.rate = 0
        self.heat = 0
        self.bibliography = bibliography
        self.subprocesses = [self]

    def __repr__(self):
        """Print the name in print()"""
        return self.name

    def __add__(self, other):
        """Sum 2 processes together: define a new process whose rates are the sum of the input process"""
        if other == 0:  # necessary for native sum() routine to work
            return self

        attrs_to_sum = "heat", "subprocesses", "network"  # all rates

        sum_process = Process()
        sum_process.rate = None  # "rate" ceases to be meaningful for composite processes
        for summed_attr in attrs_to_sum:
            attr1, attr2 = getattr(self, summed_attr), getattr(other, summed_attr)
            if attr1 is None or attr2 is None:
                setattr(sum_process, summed_attr, None)
            else:
                setattr(sum_process, summed_attr, attr1 + attr2)

        sum_process.name = f"{self.name} + {other.name}"
        return sum_process

    def __radd__(self, other):
        return self.__add__(other)

    def initialize_network(self):
        self.network = EquationSystem()  # this is a dict for which unknown keys are initialized to 0 by default

    @property
    def heat(self):
        """Energy lost from gas per unit volume and time"""
        return self.__heat

    @heat.setter
    def heat(self, value):
        """Ensures that the network is always updated when we update the heat"""
        self.__heat = value
        self.network["heat"] = Equation(d_dt(n_("heat")), value)

    def solve(
        self,
        known_quantities,
        guess,
        time_dependent=[],
        dt=None,
        model="default",
        verbose=False,
        tol=1e-3,
        careful_steps=10,
    ):
        """
        Solves the equations for a set of desired quantities given a set of known quantities

        Parameters
        ----------
        known_quantities: dict
            Dict of symbolic quantities and their values that will be plugged into the network solve as known quantities.
            Can be arrays if you want to substitute multiple values. If T is included here, we solve for chemical
            equilibrium. If T is not included, solve for thermochemical equilibrium.
        guess: dict, optional
            Dict of symbolic quantities and their values that will be plugged into the network solve as guesses for the
            unknown quantities. Can be arrays if you want to substitute multiple values. Will default to trying sensible
            guesses for recognized quantities.
        normalize_to_H: bool, optional
            Whether to return abundances normalized by the number density of H nucleons (default: True)
        reduce_network: bool, optional
            Whether to solve the reduced version of the network substituting conservation laws (default: True)
        tol: float, optional
            Desired relative error in chemical abundances (default: 1e-3)
        careful_steps: int, optional
            Number of careful initial steps in the Newton solve before full step size is used - try increasing this if
            your solve has trouble converging.

        Returns
        -------
        soldict: dict
            Dict of solved quantities
        """

        return self.network.solve(
            known_quantities,
            guess,
            time_dependent=time_dependent,
            tol=tol,
            careful_steps=careful_steps,
            dt=dt,
            verbose=verbose,
        )
