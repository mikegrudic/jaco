"""Process: the immutable atom of a microphysics model. Processes compose with + into composites."""

import sympy as sp
from .symbols import n_, d_dt
from .equation import Equation
from .equation_system import EquationSystem


def _row_lhs(key):
    """d/dt of the conserved quantity of a row"""
    return d_dt(n_(key))


class Process:
    """A set of rate equations: one row per species (number density per unit time) and the gas heat row (energy per
    unit volume and time), immutable after construction.

    Physics is built from :class:`~jaco.processes.Reaction` and :class:`~jaco.processes.ThermalTerm`; ``a + b`` is a
    composite process whose rows are the sums and whose ``subprocesses`` are the atoms of both.

    Parameters
    ----------
    name: str, optional
        Name of the process; a model uses it as the process's id.
    bibliography: sequence of str, optional
        Bibcodes (or notes) for the rates.
    rows: dict, optional
        Maps species (or "heat", or an energy reservoir such as "dust heat") to its rate expression; the heat row is 0
        if not given.
    """

    def __init__(self, name="", bibliography=(), rows=None):
        rows = rows or {}
        network = EquationSystem()
        if "heat" not in rows:
            network["heat"] = Equation(_row_lhs("heat"), 0)
        for key, rhs in rows.items():
            network[key] = Equation(_row_lhs(key), rhs)
        self._init(name, bibliography, network, None)

    def _init(self, name, bibliography, network, subprocesses):
        """Set the state and freeze; subclasses call this last"""
        self.name = name
        self.bibliography = list(bibliography)
        self._network = network
        self._subprocesses = [self] if subprocesses is None else list(subprocesses)
        self._frozen = True

    @staticmethod
    def _composite(name, bibliography, network, subprocesses):
        obj = Process.__new__(Process)
        obj._init(name, bibliography, network, subprocesses)
        return obj

    def __setattr__(self, attr, value):
        if getattr(self, "_frozen", False):
            raise AttributeError(f"{type(self).__name__} '{self.name}' is immutable; build a new process instead")
        object.__setattr__(self, attr, value)

    def __repr__(self):
        return self.name

    @property
    def network(self):
        """The rate equations, as a fresh EquationSystem the caller may modify"""
        return self._network.copy()

    @property
    def heat(self):
        """Energy given to the gas per unit volume and time"""
        return dict.__getitem__(self._network, "heat").rhs

    @property
    def subprocesses(self):
        """The atomic processes this process is the sum of"""
        return list(self._subprocesses)

    def __add__(self, other):
        """Process whose rates are the sums of the operands' rates; the operands are unchanged"""
        if isinstance(other, (int, float)) and other == 0:  # sum() starts from 0
            return self
        if not isinstance(other, Process):
            return NotImplemented
        return Process._composite(f"{self.name} +\n{other.name}", self.bibliography + other.bibliography,
                                  self._network + other._network, self._subprocesses + other._subprocesses)

    def __radd__(self, other):
        if isinstance(other, (int, float)) and other == 0:
            return self
        return NotImplemented

    def transformed(self, rule):
        """Copy with rule (a function of an expression) applied to every rate and heat expression"""
        if len(self._subprocesses) > 1:
            return sum(p.transformed(rule) for p in self._subprocesses)
        return Process(self.name, self.bibliography, {k: rule(sp.sympify(e.rhs)) for k, e in self._network.items()})

    def with_rows(self, rows):
        """Copy with rows (species or reservoir -> rate expression) added to its own, e.g. the share of its cooling a
        model deposits in a radiation band"""
        if len(self._subprocesses) > 1:
            raise ValueError("with_rows applies to an atomic process")
        new = {k: e.rhs for k, e in self._network.items()}
        for k, v in rows.items():
            new[k] = new.get(k, sp.S.Zero) + v
        return Process(self.name, self.bibliography, new)

    def unclumped(self):
        """Copy without the declared clumping factors (processes that declare none are returned unchanged)"""
        if len(self._subprocesses) > 1:
            return sum(p.unclumped() for p in self._subprocesses)
        return self

    def with_radiation(self, radiation):
        """Copy whose default clumping leaves out reactants in radiation (the species a model declares radiation)"""
        if len(self._subprocesses) > 1:
            return sum(p.with_radiation(radiation) for p in self._subprocesses)
        return self

    def solve(self, known_quantities, guess, time_dependent=[], dt=None, verbose=False, tol=1e-3, careful_steps=10):
        """Solve the network for the guessed quantities given the known ones; see :meth:`EquationSystem.solve`"""
        return self.network.solve(known_quantities, guess, time_dependent=time_dependent, tol=tol,
                                  careful_steps=careful_steps, dt=dt, verbose=verbose)

    def solver_functions(self, solve_vars, time_dependent=[], return_jac=False, return_dict=False):
        """The RHS of the reduced system and its Jacobian; see :meth:`EquationSystem.solver_functions`"""
        return self.network.solver_functions(solve_vars, time_dependent, return_jac, return_dict)
