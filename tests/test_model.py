"""Model: a keyed process collection whose declarations are worked out when the network is assembled."""

import pytest
import sympy as sp

from jaco.model import Model, Rule
from jaco.processes import Reaction, ThermalTerm, GasPhaseRecombination, CollisionalIonization
from jaco.symbols import n_, x_

k1, k2, k3, k4 = sp.symbols("k1 k2 k3 k4")
T, n_Htot = sp.Symbol("T"), sp.Symbol("n_Htot")


def toy(**declarations):
    """H/H+/H- with H- in steady state"""
    processes = [CollisionalIonization("H"), GasPhaseRecombination("H+"),
                 Reaction("H + e- -> H-", k1, name="attachment", bibliography=["t"]),
                 Reaction("H- + H -> H + H + e-", k2, name="detachment", bibliography=["t"]),
                 Reaction("H- + H+ -> 2H", k4, name="neutralization", bibliography=["t"]),
                 ThermalTerm(sp.Symbol("pdv_work"), name="PdV work")]
    return Model(processes, **{"solve_vars": ["H+"], "steady_state": ["H-"], **declarations})


def test_ids_are_unique():
    with pytest.raises(ValueError, match="duplicate"):
        Model([ThermalTerm(1, name="a"), ThermalTerm(2, name="a")])
    with pytest.raises(ValueError, match="duplicate"):
        toy() + Reaction("H + e- -> H-", k3, name="attachment", bibliography=["t"])
    with pytest.raises(ValueError, match="name"):
        Model([ThermalTerm(1)])


def test_composites_are_split_into_atoms():
    m = Model([CollisionalIonization("H") + GasPhaseRecombination("H+")])
    assert list(m.processes) == ["Collisional Ionization of H", "Gas-phase recombination of H+"]


def test_without_and_replace():
    m = toy()
    assert "detachment" not in m.without("detachment")
    with pytest.raises(KeyError):
        m.without("no such process")
    new = Reaction("H- + H -> H + H + e-", k3, name="detachment (new rate)", bibliography=["t"])
    replaced = m.replace("detachment", new)
    assert list(replaced.processes) == [*list(m.processes)[:3], "detachment (new rate)", *list(m.processes)[4:]]
    with pytest.raises(KeyError):
        m.replace("no such process", new)
    with pytest.raises(AttributeError):
        m.solve_vars = ("H",)


def test_steady_state_follows_the_processes():
    """The H- closure is worked out from the processes the model has when its network is assembled (D4)"""
    closure = toy().network.fixed_species["H-"]
    assert k3 not in closure.free_symbols
    more = toy() + Reaction("H- -> H + e-", k3, name="photodetachment", bibliography=["t"])
    assert k3 in more.network.fixed_species["H-"].free_symbols
    rhs, _ = more.solver_functions()
    assert any(k3 in sp.sympify(r).free_symbols for r in rhs)
    assert k2 not in toy().without("detachment").network.fixed_species["H-"].free_symbols


def test_nonlinear_steady_state_raises():
    m = toy() + Reaction("H- + H- -> H_2 + e- + e-", k3, name="H- pairs", bibliography=["t"])
    with pytest.raises(ValueError, match="not linear"):
        m.network


def test_merging_models():
    a = toy(fixed={"C+": 1e-4})
    b = Model([ThermalTerm(sp.Symbol("Q"), name="heating")], fixed={"C+": 1e-4}, derived={"C_2": 1 + T})
    merged = a + b
    assert "heating" in merged and merged.derived == {"C_2": 1 + T} and merged.solve_vars == ("H+",)
    with pytest.raises(ValueError, match="conflicting fixed"):
        a + Model(fixed={"C+": 2e-4})
    with pytest.raises(ValueError, match="conflicting solve_vars"):
        a + Model(solve_vars=["H"])
    with pytest.raises(ValueError, match="conflicting derived"):
        b + Model(derived={"C_2": 2})


def test_rules_rewrite_processes_at_assembly():
    halve = Rule("halve the rates", lambda p: p.transformed(lambda e: e / 2), exempt={"PdV work"})
    m = toy(rules=[halve])
    effective = {p.name: p for p in m.subprocesses}
    assert effective["attachment"].network["H-"].rhs == m.processes["attachment"].rate / 2
    assert effective["PdV work"].heat == sp.Symbol("pdv_work")
    with pytest.raises(ValueError, match="exempts processes not in the model"):
        m.without("PdV work").network


def test_reduction_reports_discarded_equations():
    red = toy().network.reduced({"T", "n_Htot"}, [])
    assert red.discarded == {"e-": "charge neutrality", "H": "H conservation", "H-": "steady state", "heat": "T is known"}


def test_solver_functions_and_solve_leave_the_network_alone():
    net = toy(solve_vars=["u", "T", "H+"], time_dependent=["T"]).network
    keys = set(net)
    net.solver_functions(["u", "T", "H+"], ["T"])
    assert set(net) == keys
    knowns = {"T": [1e4], "n_Htot": [1.0]}
    Model([CollisionalIonization("H"), GasPhaseRecombination("H+")]).solve(knowns, {"H+": [0.5]}, time_dependent=[])
    assert set(knowns) == {"T", "n_Htot"}
