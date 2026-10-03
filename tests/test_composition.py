"""Regression tests for process composition: + must not change its operands, rates must not accumulate when
reassigned, and closures derived from a model must follow later additions to it."""

import pytest
import sympy as sp

from jaco.equation import Equation
from jaco.equation_system import EquationSystem
from jaco.processes import ChemicalReaction, CollisionalIonization, GasPhaseRecombination, ThermalProcess
from jaco.symbols import d_dt, n_


D1 = pytest.mark.xfail(strict=True, reason="EquationSystem.__add__ creates the missing keys in its operands")


def snapshot(process):
    """Everything about a process's network that the reduction and code generation read"""
    net = process.network
    return {k: sp.srepr(net[k]) for k in sorted(net)}, sorted(net.chemical_species)


@D1
def test_equation_system_add_leaves_operands_unchanged():
    a = EquationSystem()
    a["heat"] = Equation(d_dt(n_("heat")), sp.Symbol("Q"))
    b = EquationSystem()
    b["H"] = Equation(d_dt(n_("H")), -sp.Symbol("k") * n_("H"))
    total = a + b
    assert set(total) == {"heat", "H"}
    assert set(a) == {"heat"} and set(b) == {"H"}


@D1
def test_process_add_leaves_operands_unchanged():
    heat = ThermalProcess(sp.Symbol("Q"), name="heating")
    recombination = GasPhaseRecombination("H+")
    before = snapshot(heat), snapshot(recombination)
    recombination + heat
    heat + recombination
    assert (snapshot(heat), snapshot(recombination)) == before
    assert set(heat.network) == {"heat"}


@D1
def test_build_A_then_B_equals_fresh_B():
    """A process shared by two models (e.g. a module-level object) must not carry model A's species into model B"""
    shared = ThermalProcess(sp.Symbol("Q"), name="shared heating")

    def build_B():
        return shared + ThermalProcess(sp.Symbol("pdv_work"), name="PdV work")

    fresh = snapshot(build_B())
    GasPhaseRecombination("H+") + CollisionalIonization("H") + shared  # model A
    assert snapshot(build_B()) == fresh


@D1
def test_model_build_leaves_its_processes_unchanged():
    """After a full model build every pure thermal term still has only the heat equation"""
    from jaco.models.starforge import make_model

    model = make_model()
    thermal = [p for p in model.subprocesses if type(p) is ThermalProcess]
    assert thermal
    assert {p.name: sorted(p.network) for p in thermal} == {p.name: ["heat"] for p in thermal}


def _reassign(process, attr, value):
    """Set attr if the process allows it; False if the process is immutable"""
    try:
        setattr(process, attr, value)
    except AttributeError:
        return False
    return True


D2_REASON = "rate setters resync the network incrementally; fixed by immutable processes (R3)"


@pytest.mark.xfail(strict=True, reason=D2_REASON)
def test_reassigning_recombination_rate_does_not_double_it():
    process = GasPhaseRecombination("H+")
    before = snapshot(process)
    if _reassign(process, "rate_coefficient", process.rate_coefficient):
        assert snapshot(process) == before


@pytest.mark.xfail(strict=True, reason=D2_REASON)
def test_reassigning_ionization_rate_does_not_double_it():
    process = CollisionalIonization("H")
    before = snapshot(process)
    if _reassign(process, "rate", process.rate):
        assert snapshot(process) == before


@pytest.mark.xfail(strict=True, reason=D2_REASON)
def test_chemical_reaction_network_follows_its_rate():
    k1, k2 = sp.symbols("k1 k2")
    reaction = ChemicalReaction("H + H -> H_2", k1, bibliography=["test"])
    if _reassign(reaction, "rate_coefficient", k2):
        assert reaction.network["H_2"].rhs == reaction.rate


@pytest.mark.slow
@pytest.mark.xfail(strict=True, reason="the H- steady state is computed once in make_model(); fixed by closure rules "
                                       "evaluated at reduction time (R1)")
def test_added_Hminus_reaction_reaches_reduced_system():
    from jaco.models.starforge import make_model, SOLVE_VARS, TIME_DEPENDENT

    k = sp.Symbol("k_test")
    model = make_model() + ChemicalReaction("H + e- -> H-", k, bibliography=["test"])
    rhs, _ = model.network.solver_functions(SOLVE_VARS, TIME_DEPENDENT)
    assert any(k in sp.sympify(r).free_symbols for r in rhs)
