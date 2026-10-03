"""The two models built from the starforge process library differ exactly by the switches in switches.py."""

import pytest
import sympy as sp
from jaco.symbols import n_, x_
from ..starforge import make_model
from ..switches import STARFORGE, STARFORGE_LEGACY
from ..symbols import n_Htot, X_H, T, G_0, Z_dust, f_dust
from ...starforge_legacy import make_model as make_legacy

SF_ONLY = {"[CI] 609 um cooling", "CO Cooling", "C+-e- Line Cooling", "C+-H Line Cooling", "H2 + HD Line Cooling",
           "Grain-assisted recombination of H+", "Charge transfer of H+ to Mg", "Direct ionization of H by cosmic rays",
           "Formation of H_2 on dust grains"}
LEGACY_ONLY = {"GIZMO C+, [CI] and CO cooling", "GIZMO H2 + HD cooling", "GIZMO H2 network"}


@pytest.fixture(scope="module")
def models():
    return make_model(STARFORGE), make_legacy()


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
    """switch rate_density: legacy rates see GIZMO's nHcgs = 0.76 rho/m_p = (0.76/X) n_Htot"""
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


def test_switch_values():
    assert STARFORGE.clumping == "all" and STARFORGE_LEGACY.clumping == "h2_chemistry"
    assert STARFORGE.h2_chemical_heat and not STARFORGE_LEGACY.h2_chemical_heat
    assert (STARFORGE.electrons, STARFORGE_LEGACY.electrons) == ("solved", "gizmo")
    with pytest.raises(ValueError):
        make_model(STARFORGE.__class__(**{**STARFORGE_LEGACY.__dict__, "h2_chemical_heat": True}))
