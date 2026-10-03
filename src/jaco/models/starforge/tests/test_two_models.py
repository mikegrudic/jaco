"""The two models built from the starforge process library differ exactly as starforge_legacy's docstring lists."""

import pytest
import sympy as sp
from jaco.symbols import n_, x_
from jaco.processes import Reaction
from ..starforge import make_model
from ..symbols import n_Htot, X_H, T, G_0, Z_dust, f_dust
from ...starforge_legacy import make_model as make_legacy, GIZMO_CLUMPING, GIZMO_DENSITY

SF_ONLY = {"[CI] 609 um cooling", "CO Cooling", "C+-e- Line Cooling", "C+-H Line Cooling", "H2 + HD Line Cooling",
           "Grain-assisted recombination of H+", "Charge transfer of H+ to Mg", "Direct ionization of H by cosmic rays",
           "Formation of H_2 on dust grains"}
LEGACY_ONLY = {"GIZMO C+, [CI] and CO cooling", "GIZMO H2 + HD cooling", "GIZMO H2 network"}


@pytest.fixture(scope="module")
def models():
    return make_model(), make_legacy()


def names(model):
    return {p.name for p in model.subprocesses}


def test_process_sets(models):
    sf, legacy = models
    assert SF_ONLY <= names(sf) and not SF_ONLY & names(legacy)
    assert LEGACY_ONLY <= names(legacy) and not LEGACY_ONLY & names(sf)
    shared = names(sf) & names(legacy)
    assert {"Gas-dust collisions", "Cosmic ray heating", "Photoelectric Heating", "Metal line cooling",
            "Nebular forbidden-line cooling", "Inverse Compton cooling (CMB)", "Gas-phase recombination of H+"} <= shared


def test_legacy_rates_at_gizmo_nHcgs(models):
    """rule GIZMO_DENSITY: legacy rates see GIZMO's nHcgs = 0.76 rho/m_p = (0.76/X) n_Htot"""
    sf, legacy = models
    lam = 0.76 / X_H
    for name, power in [("Photoelectric Heating", 1), ("Gas-phase recombination of H+", 2), ("Gas-dust collisions", 2)]:
        h_sf = next(p for p in sf.subprocesses if p.name == name).network["heat"].rhs.subs(sp.Symbol("C_2"), 1)
        h_lg = next(p for p in legacy.subprocesses if p.name == name).network["heat"].rhs
        rep = {s: lam * s for s in h_sf.free_symbols if str(s).startswith("n_")}
        vals = {T: 80.0, n_Htot: 30.0, X_H: 0.7155, G_0: 1.0, Z_dust: 1.0, f_dust: 1.0, x_("e-"): 2e-4,
                sp.Symbol("Td"): 15.0, n_("H+"): 3e-4, n_("e-"): 6e-3}
        assert float(h_lg.subs(vals)) == pytest.approx(float(h_sf.xreplace(rep).subs(vals)), rel=1e-12)
        assert float(h_lg.subs(vals)) == pytest.approx(float(h_sf.subs(vals)) * (0.76 / 0.7155) ** power, rel=0.05)
    assert X_H in set().union(*[sp.sympify(W).free_symbols for _, W in legacy.network.intermediates])


def test_declarations(models):
    sf, legacy = models
    assert sf.solve_vars == legacy.solve_vars == ("u", "T", "H+", "He+", "He++", "H_2")
    assert sf.time_dependent == legacy.time_dependent == ("T", "H_2")
    assert sf.steady_state == ("H-",) and not legacy.steady_state
    assert set(sf.fixed) == {"H_2+", "HD", "C+", "CO"} and legacy.fixed == {"C+": 0, "CO": 0}
    assert not sf.rules and legacy.rules == (GIZMO_CLUMPING, GIZMO_DENSITY)


def test_legacy_rules_apply_to_added_processes(models):
    """The legacy rules are applied when the network is assembled, so a reaction added later is unclumped and sees
    GIZMO's density like the rest"""
    _, legacy = models
    k = sp.Symbol("k_test")
    added = legacy + Reaction("H + e- -> H+ + 2e-", k, name="test ionization", bibliography=["test"])
    effective = next(p for p in added.subprocesses if p.name == "test ionization")
    rate = effective.network["H+"].rhs
    lam = 0.76 / X_H
    assert sp.simplify(rate - k * lam**2 * n_("H") * n_("e-")) == 0
    assert added.processes["test ionization"].clumping == sp.Symbol("C_2")  # the process itself is unchanged
