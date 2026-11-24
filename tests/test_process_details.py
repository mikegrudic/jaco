"""Small invariants of individual processes and helpers."""

import warnings

import sympy as sp

from jaco.processes.chemical_reaction import ChemicalReaction
from jaco.processes.recombination import Recombination
from jaco.symbols import n_, sanitize_symbols


def test_recombination_colliders_are_ordered():
    process = Recombination("H+")
    assert process.reactants == ("H+", "e-")
    assert process.nprod == n_("H+") * n_("e-")
    assert process.clumping == sp.Symbol("C_2")


def test_sanitize_symbols_renames_in_one_pass():
    xHp, ne, dt = sp.Symbol("x_H+"), sp.Symbol("n_e-"), sp.Symbol("Δt")
    expr = xHp * ne + sp.exp(-xHp / dt) + sp.Symbol("T")
    clean = {"x_H+": sp.Symbol("x_Hplus"), "n_e-": sp.Symbol("n_eminus"), "Δt": sp.Symbol("Delta_t")}
    assert sanitize_symbols(expr) == expr.subs({sp.Symbol(k): v for k, v in clean.items()})
    assert sanitize_symbols(sp.Matrix([[xHp, ne]])) == sp.Matrix([[clean["x_H+"], clean["n_e-"]]])
    assert sanitize_symbols([xHp, (ne, dt)]) == [clean["x_H+"], [clean["n_e-"], clean["Δt"]]]


def test_missing_bibliography_warns_once_per_reaction():
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        for k in (1.0, 2.0):
            ChemicalReaction("Uu + Vv -> UuVv", k)
    assert len([w for w in caught if "bibliographic reference" in str(w.message)]) == 1


def test_default_clumping_counts_reactants():
    k = sp.Symbol("k")
    assert ChemicalReaction("H_2 -> 2H", k, bibliography=["t"]).clumping == 1
    assert ChemicalReaction("H + e- -> H-", k, bibliography=["t"]).clumping == sp.Symbol("C_2")
    assert ChemicalReaction("3H -> H_2 + H", k, bibliography=["t"]).clumping == sp.Symbol("C_3")
    explicit = ChemicalReaction("H+ -> H", rate=k * n_("H+"), bibliography=["t"])
    assert explicit.clumping == 1 and explicit.rate == k * n_("H+")
    photo = ChemicalReaction("H + photon_EUV -> H+ + e-", k, bibliography=["t"])
    assert photo.clumping == sp.Symbol("C_2")  # until a model declares photon_EUV radiation
    assert photo.with_radiation({"photon_EUV"}).clumping == 1
    assert photo.with_radiation({"photon_EUV"}).rate == k * n_("H") * n_("photon_EUV")


def test_collisional_ionization_is_clumped_by_default():
    """Every two-body rate between material species carries C_2, collisional ionization included"""
    from jaco.processes import CollisionalIonization, GasPhaseRecombination

    for process in (CollisionalIonization("H"), GasPhaseRecombination("H+")):
        assert process.clumping == sp.Symbol("C_2")
        unclumped = process.unclumped()
        assert unclumped.clumping == 1 and sp.Symbol("C_2") not in unclumped.rate.free_symbols
        assert unclumped.rate == process.rate.subs(sp.Symbol("C_2"), 1)
        assert unclumped.heat == process.heat.subs(sp.Symbol("C_2"), 1)


def test_bundles_are_lists():
    from jaco.processes import LineCoolingSimple, CollisionalIonization, GasPhaseRecombination

    for bundle in (LineCoolingSimple("C+"), CollisionalIonization(), GasPhaseRecombination()):
        assert isinstance(bundle, list) and len({p.name for p in bundle}) == len(bundle) > 1
