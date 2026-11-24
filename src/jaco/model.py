"""Model: a keyed collection of processes and the declarations that make their rate equations a solvable system."""

from dataclasses import dataclass, field
from typing import Callable

import sympy as sp

from .process import Process
from .declarations import Parameter, Species, Output, Variable, merge_declarations


@dataclass(frozen=True)
class Rule:
    """A model-level rewrite of the model's processes, applied when the network is assembled, so processes added to the
    model later are rewritten too.

    Parameters
    ----------
    name: str
    apply: callable
        Maps a process to the process the model uses in its place.
    exempt: collection of str
        Ids of processes the rule leaves alone; each must be in the model.
    only: collection of str, optional
        Ids of the only processes the rule rewrites (None: all but the exempt); each must be in the model.
    """

    name: str
    apply: Callable[[Process], Process]
    exempt: frozenset = field(default_factory=frozenset)
    only: frozenset = None

    def __post_init__(self):
        object.__setattr__(self, "exempt", frozenset(self.exempt))
        if self.only is not None:
            object.__setattr__(self, "only", frozenset(self.only))

    def applies_to(self, process_id):
        return process_id not in self.exempt and (self.only is None or process_id in self.only)


def _merge_dicts(what, a, b):
    out = dict(a)
    for k, v in b.items():
        if k in out and sp.sympify(out[k]) != sp.sympify(v):
            raise ValueError(f"conflicting {what} for {k}: {out[k]} vs {v}")
        out[k] = v
    return out


def _merge_sequences(what, a, b):
    if a and b and tuple(a) != tuple(b):
        raise ValueError(f"conflicting {what}: {list(a)} vs {list(b)}")
    return tuple(a or b)


class Model:
    """A microphysics model: processes keyed by name plus everything the reduction needs, immutable after construction.

    The steady-state closures, the rules and the reductions are all worked out from the current processes when the
    network is assembled, so ``model + process`` is the model with that process in every respect.

    Parameters
    ----------
    processes: iterable of Process
        Composite processes are split into their atoms; each atom's name is its id and must be unique.
    solve_vars: sequence of str
        Solve variables in index order ("u", "T", then species).
    time_dependent: sequence of str
        Solve variables that get a backward-Euler term; the others are in steady state.
    steady_state: sequence of str
        Species eliminated by their own steady state, worked out from their rate equation when the network is
        assembled; the equation must be linear in the species' density. The abundance is clamped at 0.
    fixed: dict
        Species whose abundance is prescribed: name -> expression of the solve variables and parameters.
    derived: dict
        Parameter name -> expression of the solve variables, substituted before differentiation.
    intermediates: sequence of (Symbol, expression)
        Named quantities evaluated in order before the equations; see :class:`EquationSystem`.
    fixed_electrons: expression, optional
        Free electrons per H nucleus on species outside the network.
    rules: sequence of Rule
        Rewrites of the processes, applied in order when the network is assembled.
    parameters: sequence of Parameter
        The model's inputs besides the core ones (jaco.declarations.CORE_PARAMETERS) and those its species imply.
    species: sequence of Species
        Every species of the model with its kind. With species declared, every row of the network must be declared,
        the EOS and the conservation sums take exactly the material species, reactants declared as radiation do not
        count towards the default clumping order, and code generation refuses any undeclared symbol. The radiation and
        energy species that are not solved for are outputs.
    outputs: sequence of Output
        Further quantities the generated code evaluates at the converged state, e.g. the emission into one radiation
        band as a sum of named cooling terms.
    variables: sequence of Variable
        Solve variables that are neither u, T nor species (e.g. the dust temperature), each determined by the steady
        state of its row; each must be in solve_vars and not time-dependent.
    """

    def __init__(self, processes=(), *, solve_vars=(), time_dependent=(), steady_state=(), fixed=None, derived=None,
                 intermediates=(), fixed_electrons=None, rules=(), parameters=(), species=(), outputs=(),
                 variables=()):
        atoms = {}
        for p in processes:
            if not isinstance(p, Process):
                raise TypeError(f"not a process: {p!r}")
            for atom in p.subprocesses:
                if not atom.name:
                    raise ValueError("every process in a model needs a name, its id")
                if atom.name in atoms:
                    raise ValueError(f"duplicate process id '{atom.name}'")
                atoms[atom.name] = atom
        self._processes = atoms
        self.solve_vars = tuple(solve_vars)
        self.time_dependent = tuple(time_dependent)
        self.steady_state = tuple(steady_state)
        self.fixed = dict(fixed or {})
        self.derived = dict(derived or {})
        self.intermediates = tuple((sp.Symbol(w) if isinstance(w, str) else w, W) for w, W in intermediates)
        self.fixed_electrons = fixed_electrons
        self.rules = tuple(rules)
        self.parameters = tuple(merge_declarations("parameter", parameters).values())
        self.species = tuple(merge_declarations("species", species).values())
        self.outputs = tuple(merge_declarations("output", outputs).values())
        self.variables = tuple(merge_declarations("variable", variables).values())
        self._check()
        self._assembled = None
        self._frozen = True

    def _check(self):
        if not set(self.time_dependent) <= set(self.solve_vars):
            raise ValueError(f"time-dependent {set(self.time_dependent) - set(self.solve_vars)} not in solve_vars")
        for what, names in (("steady-state", self.steady_state), ("fixed", self.fixed)):
            if set(names) & set(self.solve_vars):
                raise ValueError(f"{what} species {set(names) & set(self.solve_vars)} are solve variables")
        if set(self.steady_state) & set(self.fixed):
            raise ValueError(f"species {set(self.steady_state) & set(self.fixed)} are both steady-state and fixed")
        if len({w for w, _ in self.intermediates}) != len(self.intermediates):
            raise ValueError("an intermediate is defined twice")
        if len({r.name for r in self.rules}) != len(self.rules):
            raise ValueError("two rules share a name")
        if not all(isinstance(p, Parameter) for p in self.parameters):
            raise TypeError("parameters must be Parameter declarations")
        if not all(isinstance(s, Species) for s in self.species):
            raise TypeError("species must be Species declarations")
        if not all(isinstance(o, Output) for o in self.outputs):
            raise TypeError("outputs must be Output declarations")
        if not all(isinstance(v, Variable) for v in self.variables):
            raise TypeError("variables must be Variable declarations")
        names = {v.name for v in self.variables}
        if names - set(self.solve_vars):
            raise ValueError(f"variables {sorted(names - set(self.solve_vars))} are not solve variables")
        if names & set(self.time_dependent):
            raise ValueError(f"variables {sorted(names & set(self.time_dependent))} are determined by a steady state")
        if names & {s.name for s in self.species}:
            raise ValueError(f"{sorted(names & {s.name for s in self.species})} declared as both species and variables")
        declared = {s.name for s in self.species}
        if declared:
            undeclared = (set(self.solve_vars) - {"u", "T"} - names | set(self.steady_state) | set(self.fixed)) - declared
            if undeclared:
                raise ValueError(f"species {sorted(undeclared)} are not declared")

    def __setattr__(self, attr, value):
        if getattr(self, "_frozen", False) and attr != "_assembled":
            raise AttributeError("Model is immutable; use +, without, replace or evolve")
        object.__setattr__(self, attr, value)

    def __repr__(self):
        return f"Model({list(self._processes)})"

    def _declarations(self):
        return dict(solve_vars=self.solve_vars, time_dependent=self.time_dependent, steady_state=self.steady_state,
                    fixed=self.fixed, derived=self.derived, intermediates=self.intermediates,
                    fixed_electrons=self.fixed_electrons, rules=self.rules, parameters=self.parameters,
                    species=self.species, outputs=self.outputs, variables=self.variables)

    def evolve(self, processes=None, **declarations):
        """Copy with the processes and/or declarations replaced"""
        unknown = set(declarations) - set(self._declarations())
        if unknown:
            raise TypeError(f"unknown declarations {unknown}")
        procs = list(self._processes.values()) if processes is None else processes
        return Model(procs, **{**self._declarations(), **declarations})

    # --- the keyed collection ---

    @property
    def processes(self):
        """id -> process, as declared (before the rules)"""
        return dict(self._processes)

    def __contains__(self, process_id):
        return process_id in self._processes

    def __add__(self, other):
        """The model with the processes of other (a Process or a Model) added; a Model's declarations are merged,
        and conflicting declarations raise"""
        if isinstance(other, (int, float)) and other == 0:
            return self
        if isinstance(other, Process):
            return self.evolve(list(self._processes.values()) + [other])
        if not isinstance(other, Model):
            return NotImplemented
        rules = {r.name: r for r in self.rules}
        for r in other.rules:
            if r.name in rules and rules[r.name] != r:
                raise ValueError(f"conflicting rules named '{r.name}'")
            rules[r.name] = r
        fe = [e for e in (self.fixed_electrons, other.fixed_electrons) if e is not None]
        if len(fe) == 2 and sp.sympify(fe[0]) != sp.sympify(fe[1]):
            raise ValueError("both models prescribe fixed_electrons")
        inter = dict(self.intermediates)
        for w, W in other.intermediates:
            if w in inter and sp.sympify(inter[w]) != sp.sympify(W):
                raise ValueError(f"conflicting definitions of intermediate {w}")
            inter[w] = W
        return Model(list(self._processes.values()) + list(other._processes.values()),
                     solve_vars=_merge_sequences("solve_vars", self.solve_vars, other.solve_vars),
                     time_dependent=_merge_sequences("time_dependent", self.time_dependent, other.time_dependent),
                     steady_state=tuple(dict.fromkeys(self.steady_state + other.steady_state)),
                     fixed=_merge_dicts("fixed abundances", self.fixed, other.fixed),
                     derived=_merge_dicts("derived parameters", self.derived, other.derived),
                     intermediates=list(inter.items()), fixed_electrons=fe[0] if fe else None,
                     rules=list(rules.values()),
                     parameters=merge_declarations("parameter", self.parameters, other.parameters).values(),
                     species=merge_declarations("species", self.species, other.species).values(),
                     outputs=merge_declarations("output", self.outputs, other.outputs).values(),
                     variables=merge_declarations("variable", self.variables, other.variables).values())

    def __radd__(self, other):
        if isinstance(other, (int, float)) and other == 0:
            return self
        if isinstance(other, Process):
            return self.evolve([other] + list(self._processes.values()))
        return NotImplemented

    def without(self, *process_ids):
        """The model without the given processes; each must be present"""
        missing = [i for i in process_ids if i not in self._processes]
        if missing:
            raise KeyError(f"no processes {missing} in the model")
        return self.evolve([p for i, p in self._processes.items() if i not in process_ids])

    def replace(self, process_id, new):
        """The model with process process_id replaced by new (a Process or a list of them), in the same place"""
        if process_id not in self._processes:
            raise KeyError(f"no process '{process_id}' in the model")
        new = list(new) if isinstance(new, (list, tuple)) else [new]
        procs = []
        for i, p in self._processes.items():
            procs += new if i == process_id else [p]
        return self.evolve(procs)

    # --- assembly ---

    @property
    def subprocesses(self):
        """The processes as the network sees them, with the rules applied"""
        return list(self._assemble()[0])

    @property
    def network(self):
        """The summed rate equations with the model's declarations attached, as a fresh EquationSystem"""
        return self._assemble()[1].copy()

    @property
    def heat(self):
        return dict.__getitem__(self._assemble()[1], "heat").rhs

    @property
    def bibliography(self):
        """id -> bibliography"""
        return {i: p.bibliography for i, p in self._processes.items()}

    def _assemble(self):
        if self._assembled is not None:
            return self._assembled
        for rule in self.rules:
            unknown = rule.exempt - set(self._processes)
            if unknown:
                raise ValueError(f"rule '{rule.name}' exempts processes not in the model: {sorted(unknown)}")
            unknown = (rule.only or frozenset()) - set(self._processes)
            if unknown:
                raise ValueError(f"rule '{rule.name}' applies to processes not in the model: {sorted(unknown)}")
        radiation = {s.name for s in self.species if s.kind == "radiation"}
        effective = []
        for i, p in self._processes.items():
            p = p.with_radiation(radiation) if radiation else p
            for rule in self.rules:
                if rule.applies_to(i):
                    p = rule.apply(p)
            effective.append(p)
        network = sum(effective).network if effective else Process().network
        if self.species:
            undeclared = set(network) - {"heat"} - {s.name for s in self.species}
            if undeclared:
                raise ValueError(f"the processes have rows for undeclared species {sorted(undeclared)}")
        network.species_kinds = {s.name: s.kind for s in self.species}
        network.species_declarations = {s.name: s for s in self.species}
        network.parameters = {p.name: p for p in self.parameters}
        for v in self.variables:
            if v.row not in network:
                raise ValueError(f"variable {v.name}: no process has the row {v.row!r} that determines it")
        network.variables = {v.name: v for v in self.variables}
        network.outputs = self._resolved_outputs(effective)
        closures = {s: network.steady_state_closure(s) for s in self.steady_state}
        network.fixed_species = {**self.fixed, **closures}
        network.steady_state = self.steady_state
        network.derived_params = dict(self.derived)
        network.intermediates = list(self.intermediates)
        network.fixed_electrons = self.fixed_electrons
        self._assembled = (effective, network)
        return self._assembled

    def _resolved_outputs(self, effective):
        """The declared outputs with heat_of resolved against the processes as assembled (rules applied)"""
        procs = dict(zip(self._processes, effective))

        def row(key):
            i, r = (key, "heat") if isinstance(key, str) else key
            net = procs[i]._network
            return dict.__getitem__(net, r).rhs if r in net else sp.S.Zero

        def total(terms):
            missing = [k for k, _ in terms if (k if isinstance(k, str) else k[0]) not in procs]
            if missing:
                raise ValueError(f"output {o.name} takes terms of processes not in the model: {missing}")
            return sum((w * row(k) for k, w in terms), sp.S.Zero)

        resolved = []
        for o in self.outputs:
            expr = o.expr.xreplace({sp.Symbol(k): total(v) for k, v in o.sums}) + total(o.heat_of)
            resolved.append(Output(o.name, expr, o.units, o.doc))
        return resolved

    # --- solving and code generation ---

    def solver_functions(self, return_jac=False, return_dict=False):
        """RHS of the reduced system in the solve variables (and its Jacobian); see EquationSystem.solver_functions"""
        return self.network.solver_functions(self.solve_vars, self.time_dependent, return_jac, return_dict)

    def solve(self, known_quantities, guess, time_dependent=None, dt=None, verbose=False, tol=1e-3, careful_steps=10):
        """Solve for the guessed quantities given the known ones (time_dependent defaults to the model's); see
        EquationSystem.solve"""
        td = self.time_dependent if time_dependent is None else time_dependent
        return self.network.solve(known_quantities, guess, time_dependent=list(td), tol=tol,
                                  careful_steps=careful_steps, dt=dt, verbose=verbose)

    def generate_code(self, output_dir=".", **kwargs):
        """Write the generated sources for the model's solve variables; see jaco.codegen.gizmo.generate_funcjac_code"""
        from .codegen.gizmo import generate_funcjac_code
        return generate_funcjac_code(self, output_dir=output_dir, **kwargs)
