"""Declared parameters, species and outputs: the contract between a model and the code that hosts it."""

import re
from dataclasses import dataclass

import sympy as sp

from .data import SolarAbundances
from .species_strings import species_counts, total_atom_abundance
from .symbols import sanitize_name

SPECIES_KINDS = ("material", "trace", "radiation", "energy")


@dataclass(frozen=True)
class Parameter:
    """An input of the model: a quantity the host provides (or the solver assumes at its default).

    Parameters
    ----------
    name: str
        Symbol name as the expressions use it (e.g. "G_0", "∇v").
    units: str
        Units, CGS unless stated.
    default: float, optional
        Value the Python solver assumes when none is given; None if the quantity must be given.
    doc: str
    """

    name: str
    units: str = ""
    default: float = None
    doc: str = ""

    @property
    def symbol(self):
        return sp.Symbol(self.name)


@dataclass(frozen=True)
class Species:
    """A species of the model, with the kind that decides how the reduction treats it.

    Kinds:

    - "material": has an element composition (parsed from its name), enters the EOS, the charge and element
      conservation sums and the clumping order of the reactions it takes part in.
    - "trace": material too dilute to count (e.g. HD, whose D the element table lacks): enters the clumping order but
      neither the EOS nor the conservation sums.
    - "radiation": photons of a band; number density per H nucleus like a species, no mass, no charge, no clumping.
    - "energy": an energy reservoir (e.g. "dust heat") whose row is an energy density.

    All but "energy" are written as abundances per H nucleus, n_s = n_Htot x_s. A radiation or energy species that is
    not solved for is an output of the generated code: its row, evaluated at the converged state.
    """

    name: str
    kind: str = "material"
    doc: str = ""

    def __post_init__(self):
        if self.kind not in SPECIES_KINDS:
            raise ValueError(f"species {self.name}: kind must be one of {SPECIES_KINDS}, not {self.kind!r}")
        if self.kind == "material" and self.name != "e-" and not species_counts(self.name):
            raise ValueError(f"material species {self.name} has no element composition")


def output_identifier(name):
    """C identifier of an output named name (a species such as "dust heat" or "photon_assoc,H")"""
    return re.sub(r"\W", "_", sanitize_name(name))


@dataclass(frozen=True)
class Output:
    """A quantity the generated code evaluates at the converged state, without derivatives, for the host to apply over
    the step: under backward Euler, Delta_t times the value is its integral over the step.

    The rows of the radiation and energy species the model does not solve for are outputs already (net production per
    unit volume and time). An Output adds a quantity that is not a row, e.g. the emission into one radiation band as
    the sum of named cooling terms.

    Parameters
    ----------
    name: str
        C identifier: the field Outputs.<name> and the index IDX_OUT_<name> of the generated code.
    expr: expression, optional
        In the model's symbols (solve variables, parameters, number densities n_X, abundances x_X); reduced like the
        rate equations.
    units: str
    doc: str
    heat_of: dict or sequence, optional
        Process id -> weight (a sequence of ids: weight 1). Adds the sum of weight * heat of those processes as the model
        assembles them, i.e. after its rules. Heat is energy given to the gas, so a cooling luminosity takes weight -1.
    """

    name: str
    expr: object = 0
    units: str = ""
    doc: str = ""
    heat_of: tuple = ()

    def __post_init__(self):
        if not re.fullmatch(r"[A-Za-z_]\w*", self.name):
            raise ValueError(f"output name {self.name!r} is not a C identifier")
        h = self.heat_of
        items = h.items() if isinstance(h, dict) else [(i, 1) if isinstance(i, str) else tuple(i) for i in h]
        object.__setattr__(self, "heat_of", tuple((str(i), sp.sympify(w)) for i, w in items))
        object.__setattr__(self, "expr", sp.sympify(self.expr))


CORE_PARAMETERS = (
    Parameter("n_Htot", "cm^-3", None, "number density of H nuclei"),
    Parameter("T", "K", None, "gas temperature, where it is not solved for"),
    Parameter("pdv_work", "erg cm^-3 s^-1", 0.0, "compressional (PdV) heating rate"),
    Parameter("y", "", SolarAbundances.x("He"), "He nuclei per H nucleus"),
    Parameter("Y", "", SolarAbundances.mass_fraction["He"], "He mass fraction"),
    Parameter("Z", "Z_sun", 1.0, "metallicity"),
    Parameter("C_2", "", 1.0, "clumping factor <n^2>/<n>^2 of two-body rates"),
    Parameter("C_3", "", 1.0, "clumping factor <n^3>/<n>^3 of three-body rates"),
)


def merge_declarations(what, *groups):
    """name -> declaration over groups of declarations; the same name declared differently raises"""
    out = {}
    for group in groups:
        for d in group:
            if d.name in out and out[d.name] != d:
                raise ValueError(f"conflicting declarations of {what} {d.name}: {out[d.name]} vs {d}")
            out[d.name] = d
    return out


def implied_parameters(species, time_dependent=()):
    """Parameters the declarations imply: the abundance of each species (an input where it is neither solved nor
    eliminated), the total abundance of each element of the material species, and the start-of-step values of the
    time-dependent variables with the step Δt"""
    out = []
    for s in species.values():
        if s.kind != "energy":
            out.append(Parameter(f"x_{s.name}", "", None, s.doc or f"{s.name} per H nucleus"))
    elements = {el for s in species.values() if s.kind == "material" for el in species_counts(s.name) if el != "e-"}
    for el in sorted(elements):
        xtot = total_atom_abundance(el)
        if isinstance(xtot, sp.Symbol):
            out.append(Parameter(str(xtot), "", None, f"{el} nuclei per H nucleus, all species"))
    if time_dependent:
        out.append(Parameter("Δt", "s", None, "time step"))
    for v in time_dependent:
        if v == "T":
            out.append(Parameter("u_initial", "erg g^-1", None, "specific internal energy at the start of the step"))
        else:
            out.append(Parameter(f"x_{v}_initial", "", None, f"{v} per H nucleus at the start of the step"))
    return out


def declared_parameter_names(parameters):
    """sanitized (C identifier) name -> Parameter"""
    return {sanitize_name(p.name): p for p in parameters.values()}
