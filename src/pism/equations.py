"""Implementation of Equation and EquationSystem for representing, manipulationg, and constructing conservation laws"""

import sympy as sp
from .symbols import d_dt, n_, x_, dt, t, BDF, n_Htot, T, internal_energy
from .data import SolarAbundances
from jax import numpy as jnp
import numpy as np
from .numerics import newton_rootsolve
from astropy import units


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

    def reduced(self, knowns, time_dependent=[]):
        subsystem = self.copy()
        subsystem.set_time_dependence(time_dependent)
        subsystem.do_conservation_reductions(knowns, time_dependent)
        if "T" in knowns and "T" not in time_dependent:
            del subsystem["heat"]
        return subsystem

    def set_time_dependence(self, time_dependent_vars):
        """Insert backward-difference formulae or set to steady state"""
        # put in backward differences
        for q in self:
            if q in time_dependent_vars:  # insert backward-difference formula
                self[q] = Equation(BDF(q), self[q].rhs)
            else:
                self[q] = Equation(0, self[q].rhs)
        if "T" in time_dependent_vars:  # special behaviour
            self["heat"] = Equation(BDF("T"), self["heat"].rhs)
            if "u" not in self:
                self["u"] = Equation(0, sp.Symbol("u") - internal_energy)

    def do_conservation_reductions(self, knowns, time_dependent_vars):
        """Eliminate equations from the system using known conservation laws."""
        self.substitutions = []

        # since we have n_Htot let's convert all other n's to x's
        for s in self.symbols:
            if "n_" in str(s) and "Htot" not in str(s):
                species = str(s).split("_")[1]
                self.substitutions.append((s, n_Htot * x_(species)))

        # charge neutrality
        if "e-" not in time_dependent_vars:
            self.substitutions.append((x_("e-"), x_("H+") + x_("He+") + 2 * x_("He++")))
            del self["e-"]

        #  general: sum(n_(species containing H) / (number of H in species))  - n_("H_2") / 2 #
        if "H+" not in time_dependent_vars:
            self.substitutions.append((x_("H+"), 1 - x_("H")))
            if "H+" in self:
                del self["H+"]

        if "He++" not in time_dependent_vars:
            y = sp.Symbol("y")
            self.substitutions.append((x_("He++"), y - x_("He") - x_("He+")))
            if "He++" in self:
                del self["He++"]

        for expr, sub in self.substitutions:
            self.subs(expr, sub)

            # general: substitute highest ionization state with n_Htot * x_element - sum of lower ionization states

    @property
    def rhs(self):
        """Return as dict of rhs-lhs instead of equations"""
        return {k: e.rhs - e.lhs for k, e in self.items()}

    @property
    def ccode(self):
        cse, cseval = sp.cse(self.rhs.values(), order="none")
        return sp.ccode(cse, standard="c99"), sp.ccode(cseval, standard="c99")

    @property
    def rhs_scaled(self):
        """Returns a scaled version of the the RHS pulling out the usual factors affecting collision rates"""
        return [r for r in self.rhs.values()]  # / (T**0.5 * n_Htot * n_Htot * 1e-12)

    def solve(
        self,
        knowns,
        guesses,
        time_dependent=[],
        dt=None,
        model="default",
        verbose=False,
        tol=1e-3,
        careful_steps=10,
        symbolic_keys=False,
    ):
        """
        Solves for equilibrium after substituting a set of known quantities, e.g. temperature, metallicity,
        etc.

        Parameters
        ----------
        known_quantities: dict
            Dict of symbolic quantities and their values that will be plugged into the network solve as known quantities.
            Can be arrays if you want to substitute multiple values. If T is included here, we solve for chemical
            equilibrium. If T is not included, solve for thermochemical equilibrium.
        guesses: dict
            Dict of symbolic quantities and their values that will be plugged into the network solve as guesses for the
            unknown quantities. Can be arrays if you want to substitute multiple values. Will default to trying sensible
            guesses for recognized quantities (NOT IMPLEMENTED YET)
        tol: float, optional
            Desired relative error in chemical abundances (default: 1e-3)
        careful_steps: int, optional
            Number of careful initial steps in the Newton solve before full step size is used - try increasing this if
            your solve has trouble converging.

        Returns
        -------
        soldict: dict
            Dict of species and their equilibrium abundances relative to H or raw number densities (depending on
            value of normalize_to_H)
        """

        def printv(*a, **k):
            """Print only if locally verbose=True"""
            if verbose:
                print(*a, **k)

        # first: check knowns and guesses are all same size
        num_params = np.array([len(np.array(guesses[g])) for g in guesses] + [len(np.array(knowns[g])) for g in knowns])
        if not np.all(num_params == num_params[0]):
            raise ValueError("Input parameters and initial guesses must all have the same shape.")
        num_params = num_params[0]

        if dt is not None:
            knowns["Δt"] = np.repeat(dt.to(units.s), num_params)

        if "u" in guesses:
            self["u"] = Equation(0, internal_energy - sp.Symbol("u"))
        subsystem = self.reduced(knowns, time_dependent)
        symbols = subsystem.symbols
        num_equations = len(subsystem)

        # are there any symbols for which we can make a reasonable assumption or directly solve the steady-state approximation?
        prescriptions = {"y": SolarAbundances.x("He"), "Y": SolarAbundances.mass_fraction["He"], "Z": 1.0}
        assumed_values = {}
        if len(symbols) > num_equations + len(knowns):
            undetermined_symbols = symbols.difference(set(sp.Symbol(g) for g in guesses))
            printv(f"Undetermined symbols: {undetermined_symbols}")
            for s in undetermined_symbols:
                # if we have a prescription for this quantity, plug it in here. This should eventually be specified at the model level.
                if str(s) in prescriptions:
                    # case 1: we have given a value, which we should add to the list of knowns
                    assumed_values[str(s)] = np.repeat(prescriptions[str(s)], num_params)
                    printv(f"{s} not specified; assuming {s}={prescriptions[str(s)]}.")
                    symbols = subsystem.symbols
                    # case 2: we have given an expression in terms of the other available quantities: we need to subs it

        # ok now we should have number of symbols unknowns + knowns
        printv(
            f"Free symbols: {symbols}\nKnown values: {list(knowns)}\nAssumed values: {list(assumed_values)}\nEquations solved: {list(subsystem.rhs)}"
        )
        if len(symbols) != len(knowns | assumed_values) + len(subsystem):
            raise ValueError(
                f"Number of free symbols is {len(symbols)} != number of knowns {len(knowns)} + number of assumptions {len(assumed_values)} + number of equations {len(subsystem)}\n"
            )
        else:
            printv(
                f"It's solvin time. Solving for {set(guesses)} based on input {set(knowns)} and assumptions about {set(assumed_values)}"
            )

        guessvals = {}
        paramvals = {}
        for s in subsystem.symbols:
            for g in guesses:
                if g == str(s) or f"x_{g}" == str(s):
                    guessvals[s] = guesses[g]
            for k in knowns | assumed_values:
                if k == str(s) or f"x_{k}" == str(s):
                    paramvals[s] = (knowns | assumed_values)[k]

        lambda_args = [list(guessvals.keys()), list(paramvals.keys())]
        func = sp.lambdify(lambda_args, subsystem.rhs_scaled, modules="jax", cse=True)

        tolerance_vars = [x_("H"), x_("He+") + x_("He"), 1 - x_("H")]
        if "T" in guesses:
            tolerance_vars += [sp.Symbol("T")]
        if "u" in guesses:
            tolerance_vars += [sp.Symbol("u"), subsystem["heat"].rhs]
            # , subsystem["heat"]]  # converge on the internal energy and  cooling rate
        tolfunc = sp.lambdify(lambda_args, tolerance_vars, modules="jax", cse=True)

        def f_numerical(X, *params):
            """JAX function to rootfind"""
            return jnp.array(func(X, params))

        def tolerance_func(X, *params):
            """Solution will terminate if the relative change in this quantity is < tol"""
            return jnp.array(tolfunc(X, params))

        # option to bail here and just provide the RHS

        # jacfunc = sp.lambdify(
        #     lambda_args, [[sp.diff(a, g) for g in guessvals] for a in subsystem.rhs_scaled]
        # )  # , modules="jax", cse=True

        sol, num_iter = newton_rootsolve(
            f_numerical,
            jnp.array([g for g in guessvals.values()]).T,
            jnp.array([p for p in paramvals.values()]).T,
            tolfunc=tolerance_func,
            rtol=tol,
            careful_steps=careful_steps,
            nonnegative=True,
            return_num_iter=True,
        )

        soldict = self.package_solution(sol, guessvals, guesses, paramvals, subsystem, symbolic_keys)

        return soldict

    def package_solution(self, sol, guessvals, guesses, paramvals, subsystem, symbolic_keys):
        # now repack the solution
        soldict = {}
        for i, g in enumerate(guessvals):
            soldict[g] = sol[:, i]
        # do a reverse-pass on the substitutions we made to get all quantities
        values_to_subs = soldict | paramvals
        for expr, sub in reversed(subsystem.substitutions):
            if expr in soldict:
                continue
            if "n_" in str(expr):
                continue
            soldict[expr] = sp.lambdify(list(sub.free_symbols), sub)(
                *[values_to_subs[s] for s in list(sub.free_symbols)]
            )  # should probably make a function of this
            values_to_subs |= soldict
        if not symbolic_keys:
            soldict = {str(k): v for k, v in soldict.items()}
            # if we have a bunch of x_'s, should also link up keys in the original input format, e.g. H->x_H
            if np.any(["x_" in k for k in guesses]):  # if we specified abundances with x_ notation, return same
                return soldict
            soldict2 = {}  # otherwise return with input format where keys are simple species strings
            for k in soldict:
                if "x_" in str(k):
                    soldict2[str(k).replace("x_", "")] = soldict[k]
                else:
                    soldict2[k] = soldict[k]
            soldict = soldict2

        return soldict
