"""Small invariants of individual processes and helpers."""

import warnings

import sympy as sp

from jaco.processes.chemical_reaction import ChemicalReaction
from jaco.processes.recombination import Recombination
from jaco.symbols import n_, sanitize_symbols


def test_recombination_colliders_are_ordered():
    process = Recombination("H+")
    assert process.colliding_species == ("H+", "e-")
    assert process.nprod == n_("H+") * n_("e-")
    assert process.clumping_factor == sp.Symbol("C_2")


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
