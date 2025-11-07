"""Implementation of Equation and EquationSystem for representing, manipulationg, and constructing conservation laws"""

import sympy as sp
from .symbols import d_dt, n_, dt, t, BDF
from .data import SolarAbundances


class Equation(sp.core.relational.Equality):
    """Sympy equation where we overload addition/subtraction to apply those operations to the RHS, for summing rate
    equations"""

    def get_summand(self, other):
        """Value-check the operand and return the quantity to be summed in the operation: the expression itself if an expression, or the RHS"""
        if isinstance(other, sp.core.relational.Equality):
            if self.lhs != other.lhs:
                raise ValueError(
                    "Tried to sum incompatible equations. Equation summation only defined for differential equations with the same LHS."
                )
            else:
                return other.rhs
        elif isinstance(other, sp.logic.boolalg.BooleanAtom):
            return 0
        else:
            return other

    def __add__(self, other):
        summand = self.get_summand(other)
        return Equation(self.lhs, self.rhs + summand)

    def __sub__(self, other):
        summand = self.get_summand(other)
        return Equation(self.lhs, self.rhs - summand)

    def __radd__(self, other):
        return self.__add__(other)

    def __iadd__(self, other):
        self = self + other
        return self

    def __isub__(self, other):
        self = self - other
        return self


class EquationSystem(dict):
    """Dict of symbolic expressions with certain superpowers for manipulating sets of conservation equations."""

    def copy(self):
        new = EquationSystem()
        for k in self:
            new[k] = self[k]
        return new

    def __getitem__(self, __key: str):
        """Dict getitem method where we initialize a differential equation for the conservation of a species if the key
        does not exist"""
        if __key not in self:
            self.__setitem__(__key, Equation(d_dt(n_(__key)), 0))  # technically should only be n_ if this is a species
            # need to make sure that d/dt's don't add up when composing equations
        return super().__getitem__(__key)

    def __add__(self, other):
        """Return a dict whose values are the sum of the values of the operands"""
        keys = self.keys() | other.keys()
        new = EquationSystem()
        for k in keys:
            new[k] = self[k] + other[k]
        return new

    @property
    def symbols(self):
        """Returns the set of all symbols in the equations"""
        all = set()
        for e in self.values():
            all.update(e.free_symbols)
        if t in all:  # leave time out
            all.remove(t)
        return all

    @property
    def jacobian(self):
        """Returns a dict of dicts representing the Jacobian of the RHS of the system. Keys are the names of the
        conserved quantities and subkeys are the variable of differentiation.
        """
        return {k: {s: sp.diff(e.rhs, s) for s in self.symbols} for k, e in self.items()}

    def subs(self, expr, replacement):
        """Substitute symbolic expressions throughout the whole network."""
        for k, e in self.items():
            self[k] = e.subs(expr, replacement)

    def reduced_system(self, knowns, time_dependent=[]):
        subsystem = self.copy()
        subsystem.eulerify(time_dependent)
        subsystem.do_conservation_reductions(knowns, time_dependent)
        return subsystem

    def eulerify(self, time_dependent_vars):
        """Insert backward-difference formula for the"""
        # put in backward differences
        for q in self:
            if q in time_dependent_vars:  # insert backward-difference formula
                self[q] = Equation(BDF(q), self[q].rhs)
            else:
                self[q] = Equation(0, self[q].rhs)

    def do_conservation_reductions(self, knowns, time_dependent_vars):
        """Eliminate equations from the system using known conservation laws."""
        # charge neutrality
        if "e-" not in time_dependent_vars:
            self.subs(n_("e-"), n_("H+") + n_("He+") + 2 * n_("He++"))  # general: sum(n_species * ion charge)
            del self["e-"]

        if "n_Htot" in knowns:
            n_Htot = sp.Symbol("n_Htot", positive=True)
            #  general: sum(n_(species containing H) / (number of H in species))  - n_("H_2") / 2 #
            if "H+" not in time_dependent_vars:
                self.subs(n_("H+"), n_Htot - n_("H"))
                if "H+" in self:
                    del self["H+"]

            if "He++" not in time_dependent_vars:
                Y = sp.Symbol("Y")
                y = Y / (1 - Y) / 4
                self.subs(n_("He++"), n_Htot * y - n_("He") - n_("He+"))
                if "He++" in self:
                    del self["He++"]

            # general: substitute highest ionization state with n_Htot * x_element - sum of lower ionization states
            # x_Atot = SolarAbundances.x(A)
            # substitutions[highest ionization state of A] = n_Htot * x_Atot - sum(n_A(i))
            # do we have to write a chemical species class?

    @property
    def rhs(self):
        """Return as dict of rhs-lhs instead of equations"""
        return {k: e.rhs - e.lhs for k, e in self.items()}

    @property
    def ccode(self):
        cse, cseval = sp.cse(self.rhs.values(), order="none")
        return sp.ccode(cse, standard="c99"), sp.ccode(cseval, standard="c99")

    def solve(self, knowns, unknowns, time_dependent=[]):
        subsystem = self.reduced_system(knowns, time_dependent)
        if "T" in knowns:
            del subsystem["T"]

    # @property
    # def steadystate(self, species=None):
    #     """Returns the system with all time derivatives set to 0"""
    #     return {k: Equation(0, e.rhs) for k, e in self.items()}  # steadystate_equations

    # @property
    # def network_ions(self):
    #     """Returns the list of ions involved in a process"""
    #     return [s for s in self.network if is_an_ion(s)]

    # def apply_network_reductions(self, expr):
    #     """Applies the replacements given by network_reduction_replacements to a symbolic expression"""
    #     out = expr
    #     for _ in range(2):  # 2 passes to avoid ordering issues
    #         for n, r in self.network_reduction_replacements.items():
    #             out = out.subs(n, r)
    #     return out

    # @property
    # def reduced_network(self):
    #     """
    #     Returns the chemistry network after substituting known conservation laws:

    #     n_atom = sum(n_{species containing atom} * number of atoms in species)
    #     n_e- = sum(ion charge * n_ion) - want to keep n_e- in the explicit updates, so eliminate the highest ions?

    #     This reduces the network of N rate equations to N - (num_atoms + 1).
    #     """

    #     replacements = self.network_reduction_replacements

    #     reduced_network = {}
    #     for s, rhs in self.network.items():
    #         if n_(s) in replacements:
    #             continue
    #         else:
    #             rhs = self.apply_network_reductions(rhs)
    #         reduced_network[s] = rhs
    #     return reduced_network

    # def get_thermochem_network(self, reduced=True):
    #     """Returns the network including all chemical processes plus the gas heating-cooling equation"""
    #     if reduced:
    #         network = self.reduced_network
    #     else:
    #         network = self.network
    #     return network | {"T": self.apply_network_reductions(self.heat)}  # combine the dicts


# def eulerify(symbol):


# @property
# def network_reduction_replacements(self):
#     """Replacements for reducing the chemistry network with conservation laws"""
#     Y = sp.Symbol("Y")  # general: mass fractions of different atoms other than H. n_i,tot = n_i + sum(n_ions of i)
#     nHtot = sp.Symbol("n_Htot")  # basically always want this

#     substitutions = {
#         n_("e-"): n_("H+") + n_("He+") + 2 * n_("He++"),
#         # n_("He+"): n_("e-") - n_("H+") - 2 * n_("He++"),
#         n_("H+"): nHtot - n_("H"),
#         n_("He++"): Y / (4 - 4 * Y) * nHtot - sp.Symbol("n_He") - sp.Symbol("n_He+"),
#     }
#     return substitutions

# def apply_network_reductions(self, expr):
#     """Applies the replacements given by network_reduction_replacements to a symbolic expression"""
#     out = expr
#     for _ in range(2):  # 2 passes to avoid ordering issues
#         for n, r in self.network_reduction_replacements.items():
#             out = out.subs(n, r)
#     return out

# @property
# def reduced_network(self):
#     """
#     Returns the chemistry network after substituting known conservation laws:

#     n_atom = sum(n_{species containing atom} * number of atoms in species)
#     n_e- = sum(ion charge * n_ion) - want to keep n_e- in the explicit updates, so eliminate the highest ions?

#     This reduces the network of N rate equations to N - (num_atoms + 1).
#     """

#     replacements = self.network_reduction_replacements

#     reduced_network = {}
#     for s, rhs in self.network.items():
#         if n_(s) in replacements:
#             continue
#         else:
#             rhs = self.apply_network_reductions(rhs)
#         reduced_network[s] = rhs
#     return reduced_network

# def get_thermochem_network(self, reduced=True):
#     """Returns the network including all chemical processes plus the gas heating-cooling equation"""
#     if reduced:
#         network = self.reduced_network
#     else:
#         network = self.network
#     return network | {"T": self.apply_network_reductions(self.heat)}  # combine the dicts


#    def backward_eulerfy

# def jacobian

# def numerify

# def numerify_jacobian
