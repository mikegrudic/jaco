"""Sources and sinks of H+ in the two models.

STARFORGE_LEGACY: GIZMO's H balance (find_abundances_and_rates, eqn 33: nH0 = aHp / (aHp + geH0 + gJH0ne)) has collisional
and UV-background/RT photoionization only; its cosmic-ray ionization of the neutrals is the heavy-ion electron term.
STARFORGE: cosmic rays ionize atomic H at the attenuated rate GIZMO uses for its CR heating and heavy ions, and H+ also
recombines on grains (WD01) and by charge transfer to Mg."""

import pytest
import sympy as sp
from jaco.symbols import n_
from ..starforge import make_model
from ...starforge_legacy import make_model as make_legacy
from ..symbols import cosmicray_ionization_rate_H, x_


def H_plus_terms(model):
    sources, sinks = {}, {}
    for p in model.subprocesses:
        if "H+" not in p.network or p.network["H+"].rhs == 0:
            continue
        rhs = p.network["H+"].rhs
        made = sp.simplify(rhs.subs(n_("H+"), 0))
        (sources if made != 0 else sinks)[p.name] = rhs
    return sources, sinks


def test_legacy_only_collisions_make_Hplus():
    sources, sinks = H_plus_terms(make_legacy())
    assert list(sources) == ["Collisional Ionization of H"]
    assert list(sinks) == ["Gas-phase recombination of H+"]


def test_starforge_cosmic_rays_ionize_H():
    sources, sinks = H_plus_terms(make_model())
    assert sp.simplify(sources.pop("Direct ionization of H by cosmic rays") - cosmicray_ionization_rate_H * n_("H")) == 0
    assert list(sources) == ["Collisional Ionization of H"]
    for name in ("Gas-phase recombination of H+", "Grain-assisted recombination of H+", "Charge transfer of H+ to Mg"):
        assert name in sinks
