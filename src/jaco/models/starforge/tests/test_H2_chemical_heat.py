"""Heat of H2 formation and collisional dissociation.

STARFORGE gives the gas the formation heat partitioned by the critical density of H2, as compiled by Nickerson, Teyssier
& Rosdahl 2018 (Eqs. 46-47: Hollenbach & McKee 1979 for grains and H-, Omukai 2000 for three-body formation), and takes
4.48 eV per collisional dissociation. The rest of the binding energy is radiated, under the optically-thick limit.
GIZMO's cooling module (STARFORGE_LEGACY) carries neither."""

import numpy as np
import pytest
import sympy as sp
from jaco.symbols import n_
from ..h2_chemistry.H2_chemistry import h2_chemistry_processes, h2_formation_radiation
from ..h2_chemistry.chemical_heat import RADIATION_OF
from ..h2_chemistry.gizmo_network import GizmoH2Network
from ..symbols import T, n_Htot, x_, T_dust, f_dust, Z_dust

EV = 1.602176634e-12
ROUNDOFF = 1e-14 * 4.48 * EV  # of the heat per H2 formed: the binding energy less its radiated part


def n_cr_HM79(Tv, xH, xH2):
    return 1e6 / np.sqrt(Tv) / (1.6 * xH * np.exp(-((400 / Tv) ** 2)) + 1.4 * xH2 * np.exp(-12000 / (Tv + 1200)))


def by_name(processes, start):
    return next(p for p in processes if p.name.startswith(start))


STATES = [(30.0, 1e2, 0.4, 0.3), (300.0, 1e4, 0.1, 0.45), (3000.0, 1e9, 0.9, 0.05), (100.0, 1.0, 1.0, 1e-6)]


@pytest.mark.parametrize("Tv,nH,xH,xH2", STATES)
@pytest.mark.parametrize("start,channel_eV", [("Formation of H_2 on dust grains", lambda f: 0.2 + 4.2 * f),
                                              ("Associative detachment of H with H-", lambda f: 3.53 * f),
                                              ("3-body formation of H_2", lambda f: 4.48 * f)])
def test_formation_heat_per_event(Tv, nH, xH, xH2, start, channel_eV):
    """the reaction gives the gas the binding energy and its radiation takes back all but HM79's or Omukai's heat"""
    chemistry = h2_chemistry_processes(True)
    p = by_name(chemistry, start)
    rad = by_name(h2_formation_radiation(chemistry), RADIATION_OF + start)
    f = 1 / (1 + n_cr_HM79(Tv, xH, xH2) / nH)
    vals = {T: Tv, n_Htot: nH, x_("H"): xH, x_("H_2"): xH2, T_dust: 15.0, f_dust: 1.0, Z_dust: 1.0, sp.Symbol("C_2"): 1.3,
            sp.Symbol("C_3"): 1.3**3, n_("H"): xH * nH, n_("H-"): 1e-9 * nH, n_("e-"): 1e-4 * nH}
    assert float((p.heat / p.rate).subs(vals)) == pytest.approx(4.48 * EV, rel=1e-12, abs=0)
    heat_per_event = float(((p.heat + rad.heat) / p.rate).subs(vals))
    assert heat_per_event == pytest.approx(channel_eV(f) * EV, rel=1e-10, abs=ROUNDOFF)


@pytest.mark.parametrize("collider", ["H+", "e-", "H_2", "H", "He"])
def test_dissociation_takes_binding_energy(collider):
    p = by_name(h2_chemistry_processes(True), f"Collisional dissociation of H_2 by {collider}")
    assert sp.simplify(p.heat / p.rate + 4.48 * EV) == 0


def test_no_chemical_heat_in_legacy():
    assert all(p.heat == 0 for p in h2_chemistry_processes(False))
    assert h2_formation_radiation(h2_chemistry_processes(False)) == []
    assert GizmoH2Network().network["heat"].rhs == 0


def test_three_body_rate_fit_is_a_keyword():
    from ..h2_chemistry.three_body import H2_3body_formation, threebody_rate
    for fit in (0, 4, 5):
        p = H2_3body_formation(False, fit=fit)
        assert sp.simplify(p.rate_coefficient - threebody_rate(fit)[0]) == 0
        assert p.bibliography == [threebody_rate(fit)[1]]


def test_one_reaction_can_be_left_out():
    """The H2 chemistry is a list of reactions with unique names, not a pre-summed bundle"""
    processes = h2_chemistry_processes(True)
    names = [p.name for p in processes]
    assert len(set(names)) == len(names)
    rest = sum(p for p in processes if p.name != "H+ + H- -> H_2+ + e-")
    assert "H_2+" not in rest.network


# --- the formation radiation under the optically-thick limit (starforge.make_model) ---

FORMATION = ("Formation of H_2 on dust grains", "Associative detachment of H with H-", "3-body formation of H_2")


@pytest.fixture(scope="module")
def assembled():
    """the starforge model's processes as it assembles them (the optically-thick limit's factor f_thick on the capped)"""
    from ..starforge import make_model
    return {p.name: p for p in make_model().subprocesses}


def formation_cycle(processes, start, f_thick, Tv, nH, xH, xH2):
    """(heat per H2 formed by every process booked to the formation, the radiation's share of it), f_thick given"""
    from ..optically_thick import F_THICK
    r = processes[start]
    terms = [p for name, p in processes.items() if name in (start, RADIATION_OF + start)]
    vals = {T: Tv, n_Htot: nH, x_("H"): xH, x_("H_2"): xH2, T_dust: Tv, f_dust: 1.0, Z_dust: 1.0, sp.Symbol("C_2"): 1.0,
            sp.Symbol("C_3"): 1.0, n_("H"): xH * nH, n_("H-"): 1e-12 * nH, n_("e-"): 1e-8 * nH, F_THICK: f_thick}
    rad = processes.get(RADIATION_OF + start)
    return (float((sum(p.heat for p in terms) / r.rate).subs(vals)),
            float((rad.heat / r.rate).subs(vals)) if rad is not None else 0.0)


@pytest.mark.parametrize("Tv,nH,xH,xH2", [(1800.0, 7e13, 1e-3, 0.4995), (300.0, 1e4, 0.1, 0.45)])
@pytest.mark.parametrize("start", FORMATION)
def test_cycle_conserves_energy_where_the_limit_binds(assembled, start, Tv, nH, xH, xH2):
    """f_thick = 1e-6 (the cooling 1e6 x the limit): a formation then a collisional dissociation leave the gas with
    1e-6 of the radiated remainder"""
    heat, radiated = formation_cycle(assembled, start, 1e-6, Tv, nH, xH, xH2)
    assert abs(heat - 4.48 * EV) <= 1e-6 * 4.48 * EV
    assert radiated <= 0


@pytest.mark.parametrize("Tv,nH,xH,xH2", STATES)
@pytest.mark.parametrize("start,channel_eV", [("Formation of H_2 on dust grains", lambda f: 0.2 + 4.2 * f),
                                              ("Associative detachment of H with H-", lambda f: 3.53 * f),
                                              ("3-body formation of H_2", lambda f: 4.48 * f)])
def test_remainder_radiates_where_the_limit_does_not_bind(assembled, start, channel_eV, Tv, nH, xH, xH2):
    """f_thick = 1: the gas gets HM79's or Omukai's heat per H2 formed, and the rest of 4.48 eV is radiated"""
    f = 1 / (1 + n_cr_HM79(Tv, xH, xH2) / nH)
    heat, radiated = formation_cycle(assembled, start, 1.0, Tv, nH, xH, xH2)
    assert heat == pytest.approx(channel_eV(f) * EV, rel=1e-10, abs=ROUNDOFF)
    assert -radiated == pytest.approx((4.48 - channel_eV(f)) * EV, rel=1e-10, abs=0)
