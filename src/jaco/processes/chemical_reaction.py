"""Reaction: a process given by a chemical equation and its rate."""

import sympy as sp
from warnings import warn
from ..process import Process
from ..symbols import n_

_warned_unreferenced = set()  # equations already reported as lacking a bibliography


class Reaction(Process):
    """A chemical reaction, immutable after construction.

    The events per unit volume and time are ``clumping * rate_coefficient * prod(n_reactant)`` (mass action), or
    ``clumping * rate`` for a rate that is not of mass-action form. Each event removes the reactants, adds the products
    and gives the gas ``heat_per_reaction``.

    Parameters
    ----------
    equation: str
        Syntax: ``<c1>reactant1 + <c2>reactant2 ... -> <c3>product1 + ...``, e.g. ``'3H -> H_2 + H'``, ``'H+ + e- -> H'``.
        Species on both sides (colliders, catalysts) keep their rows.
    rate_coefficient: expression, optional
        k of the mass-action rate.
    heat_per_reaction: expression, optional
        Energy given to the gas per event (erg); negative for cooling.
    rate: expression, optional
        Events per unit volume and time before clumping, instead of a rate coefficient.
    heat_rate_coefficient: expression, optional
        Heat per unit volume and time per prod(n_reactant) (erg cm^3(n-1) s^-1), instead of heat_per_reaction, for fits
        that give the cooling separately from the rate.
    clumping: expression, optional
        Factor <n^k>/<n>^k multiplying the events and heat. Default: C_k for a mass-action rate with k >= 2 reactants
        (a model may exclude reactants it declares to be radiation), 1 otherwise.
    name: str, optional
        Defaults to the equation.
    bibliography: sequence of str, optional
    """

    def __init__(self, equation, rate_coefficient=None, heat_per_reaction=0, *, rate=None, heat_rate_coefficient=None,
                 clumping=None, name="", bibliography=()):
        if rate is not None and rate_coefficient is not None:
            raise ValueError(f"{equation}: give a rate coefficient or a rate, not both")
        if heat_rate_coefficient is not None and (heat_per_reaction != 0 or rate is not None):
            raise ValueError(f"{equation}: heat_rate_coefficient goes with a rate coefficient and replaces heat_per_reaction")
        if rate is None and rate_coefficient is None:
            rate_coefficient = 0.0
        if not bibliography and equation not in _warned_unreferenced:
            _warned_unreferenced.add(equation)
            warn(f"Chemical reaction {equation} does not have a bibliographic reference. Be a lot cooler if it did.",
                 stacklevel=2)
        self.equation = equation
        self.lhs_coeffs, self.rhs_coeffs = self.species_and_coeffs(equation)
        self.reactants = tuple(s for s, c in self.lhs_coeffs.items() for _ in range(c))
        self.rate_coefficient = rate_coefficient
        self.heat_per_reaction = heat_per_reaction
        self.heat_rate_coefficient = heat_rate_coefficient
        self._explicit_rate = rate
        self._clumping_given = clumping
        self.clumping = self.default_clumping() if clumping is None else clumping

        nprod = self.nprod
        if rate is None:
            self.rate = rate_coefficient * nprod * self.clumping
            if heat_rate_coefficient is not None:
                heat = heat_rate_coefficient * nprod * self.clumping
            else:
                heat = rate_coefficient * heat_per_reaction * nprod * self.clumping
        else:
            self.rate = rate * self.clumping
            heat = heat_per_reaction * self.rate

        rows = {}
        for s, coeff in self.lhs_coeffs.items():
            rows[s] = rows.get(s, sp.S.Zero) - self.rate * coeff
        for s, coeff in self.rhs_coeffs.items():
            rows[s] = rows.get(s, sp.S.Zero) + self.rate * coeff
        rows["heat"] = heat
        super().__init__(name or equation, bibliography, rows)

    def default_clumping(self, radiation=()):
        """C_k for a mass-action rate with k >= 2 reactants that are not in radiation, else 1"""
        k = sum(1 for s in self.reactants if s not in radiation)
        return sp.Symbol(f"C_{k}") if self._explicit_rate is None and k > 1 else 1

    @property
    def nprod(self):
        """Product of the reactant number densities"""
        return sp.prod([n_(s) for s in self.reactants])

    def _rebuild(self, **changes):
        args = dict(rate_coefficient=None if self._explicit_rate is not None else self.rate_coefficient,
                    heat_per_reaction=self.heat_per_reaction, rate=self._explicit_rate,
                    heat_rate_coefficient=self.heat_rate_coefficient, clumping=self._clumping_given, name=self.name,
                    bibliography=self.bibliography)
        args.update(changes)
        return Reaction(self.equation, **args)

    def unclumped(self):
        return self if self.clumping == 1 else self._rebuild(clumping=1)

    def with_radiation(self, radiation):
        """Copy whose default clumping excludes reactants in radiation (species a model declares to be radiation)"""
        if self._clumping_given is not None or self.default_clumping(radiation) == self.clumping:
            return self
        return self._rebuild(clumping=self.default_clumping(radiation))

    @staticmethod
    def equation_to_heat(equation):
        """Given an equation, return the enthalpy of the reaction in erg based upon chemical data"""
        raise NotImplementedError("Automatic reaction enthalpy not implemented yet.")

    @staticmethod
    def species_and_coeffs(equation) -> list[dict]:
        """Returns the lists of species and corresponding stoichiometric coefficients from the equation"""
        if "->" not in equation:
            raise ValueError(f"Chemical equation {equation} has no ->")
        # TODO: add check that the equation is balanced
        coeffs_dicts = []
        for idx, side in enumerate(("lhs", "rhs")):
            terms = equation.split("->")[idx].split(" + ")
            terms = [t.strip() for t in terms]
            coefficients = len(terms) * [1]
            species = terms.copy()

            # strip off leading coefficients on the species
            for i, species_string in enumerate(terms):
                for j in range(1, len(species_string)):
                    substr = species_string[:j]
                    if substr.isnumeric():
                        coefficients[i] = int(substr)
                        species[i] = species_string[j:]
                    else:
                        break

            # if we have duplicates, must sum the coefficients, otherwise we lose that info when we make the dict
            coeffdict = {s: 0 for s in species}
            for i, s in enumerate(species):
                coeffdict[s] += coefficients[i]

            coeffs_dicts.append(coeffdict)
        return coeffs_dicts


ChemicalReaction = Reaction
