"""Heat of H2 formation and collisional dissociation.

STARFORGE gives the gas the formation heat partitioned by the critical density of H2, as compiled by Nickerson, Teyssier
& Rosdahl 2018 (Eqs. 46-47: Hollenbach & McKee 1979 for grains and H-, Omukai 2000 for three-body formation), and takes
4.48 eV per collisional dissociation. GIZMO's cooling module (STARFORGE_LEGACY) carries neither."""

import numpy as np
import pytest
import sympy as sp
from jaco.symbols import n_
from ..h2_chemistry.H2_chemistry import h2_chemistry_processes
from ..h2_chemistry.gizmo_network import GizmoH2Network
from ..symbols import T, n_Htot, x_, T_dust, f_dust, Z_dust

EV = 1.602176634e-12


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
    p = by_name(h2_chemistry_processes(True), start)
    f = 1 / (1 + n_cr_HM79(Tv, xH, xH2) / nH)
    vals = {T: Tv, n_Htot: nH, x_("H"): xH, x_("H_2"): xH2, T_dust: 15.0, f_dust: 1.0, Z_dust: 1.0, sp.Symbol("C_2"): 1.3,
            sp.Symbol("C_3"): 1.3**3, n_("H"): xH * nH, n_("H-"): 1e-9 * nH, n_("e-"): 1e-4 * nH}
    heat_per_event = float((p.heat / p.rate).subs(vals))
    assert heat_per_event == pytest.approx(channel_eV(f) * EV, rel=1e-10)


@pytest.mark.parametrize("collider", ["H+", "e-", "H_2", "H", "He"])
def test_dissociation_takes_binding_energy(collider):
    p = by_name(h2_chemistry_processes(True), f"Collisional dissociation of H_2 by {collider}")
    assert sp.simplify(p.heat / p.rate + 4.48 * EV) == 0


def test_no_chemical_heat_in_legacy():
    assert all(p.heat == 0 for p in h2_chemistry_processes(False))
    assert GizmoH2Network().network["heat"].rhs == 0
