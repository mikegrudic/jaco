"""The dust temperature of the non-RT models against GIZMO's per-iterate estimate (SINGLE_STAR_SINK_DYNAMICS):
get_equilibrium_dust_temperature_estimate -> rt_eqm_dust_temp -> dust_dEdt (cooling.cc 2335-2380, rt_utilities.cc
1580-1690, rt_dust_opacity.cc 66-100, eos.cc 284-290 at gizmo_jaco_collapse 1d6267ca), transcribed below.

The model's Td is the root of its dust heat row at the gas temperature of the solve. GIZMO's is the root of dust_dEdt at
each temperature its cooling evaluates; the reference here is that function's exact root (GIZMO's own iteration stops
at 1e-3 of T - Td)."""

import math

import numpy as np
import pytest
import sympy as sp
from scipy.optimize import brentq

from jaco.symbols import x_, n_, sanitize_symbols
from jaco.model_check import _c_lambdify
from jaco.interpolation import tables_in
from ..symbols import T, n_Htot, X_H
from .. import dust_opacity as do
from ..starforge import make_model as make_starforge
from ...starforge_legacy import make_model as make_legacy
from ...starforge_legacy_RT_EUV import make_model as make_euv
from .test_legacy_rt import gizmo_planck_mean

S = sp.Symbol
MP, X_GIZMO, SIGMA = 1.6726e-24, 0.76, 5.67e-5
PROTONMASS_JACO = 1.67262178e-24  # the host's Sigma -> N_H_cool conversion


@pytest.fixture(scope="module")
def models():
    return {"starforge": make_starforge(), "starforge_legacy": make_legacy(), "starforge_legacy_RT_EUV": make_euv()}


# --- GIZMO transcription (non-RT, SINGLE_STAR_SINK_DYNAMICS) ------------------------------------------------------

def gizmo_kappa_dust(Td, radtransfer=False):
    """rt_kappa_adaptive_IR_band(Td, Td, 1, 1) per solar metallicity [cm^2/g]: dust_planck_mean_opacity, which
    without RT_INFRARED returns ~0 (MIN_REAL_NUMBER) at or above 1500 K, times 1 - albedo under RADTRANSFER"""
    if Td >= 1500.0:
        return 0.0
    k = gizmo_planck_mean(Td, Td)
    return k * (1 - 0.5 / (1 + 725.0**2 / (1 + Td**2))) if radtransfer else k


def gizmo_dust_dEdt(Td, Tg, nHcgs, rho, Zsol, A, Sigma, radtransfer=False):
    """dust_dEdt [erg cm^-3 s^-1]: gas_dust_heating_coeff (untruncated, dust-to-metals 1) nHcgs^2 (T - Td) + A
    - 4 sigma rho kappa Td^4 / (1 + tau^2), tau = Sigma kappa"""
    coupling = 1.116e-32 * math.sqrt(Tg) * (1 - 0.8 * math.exp(-75.0 / Tg)) * Zsol * nHcgs**2 * (Tg - Td)
    kappa = gizmo_kappa_dust(Td, radtransfer) * Zsol
    return coupling + A - 4 * SIGMA * rho * kappa * Td**4 / (1 + (Sigma * kappa) ** 2)


def gizmo_dust_temperature(Tg, nHcgs, rho, Zsol, A, Sigma, radtransfer=False, lo=2.73, hi=1e4):
    """The root of dust_dEdt in [2.73, 1e4] K that a walk from the gas temperature reaches (at a bound, the bound)"""
    f = lambda t: gizmo_dust_dEdt(t, Tg, nHcgs, rho, Zsol, A, Sigma, radtransfer)
    return _walk_root(f, min(max(Tg, lo), hi), lo, hi)


def _walk_root(f, guess, lo=2.73, hi=1e4):
    t, ft = guess, f(guess)
    if ft == 0:
        return t
    fac, up = 1.1, ft > 0
    while True:
        t_new = min(hi, t * fac) if up else max(lo, t / fac)
        f_new = f(t_new)
        if (f_new > 0) != up or f_new == 0:
            a, b = sorted((t, t_new))
            return math.exp(brentq(lambda y: f(math.exp(y)), math.log(a), math.log(b), xtol=1e-14))
        if t_new in (lo, hi):
            return t_new
        t, fac = t_new, fac * 1.1


# --- the model's dust temperature ---------------------------------------------------------------------------------

def values(Tg, n, Zsol, A, Sigma, f_thick=1.0):
    """Symbol -> value of a neutral, partly molecular state; X = 0.7155, so nHcgs = 0.76 n_Htot / X"""
    X = 0.7155
    v = {T: Tg, n_Htot: n, X_H: X, S("Z_d"): Zsol, S("f_d"): 1.0, S("A_dust"): A, S("rho"): n * MP / X,
         S("N_H_cool"): Sigma / PROTONMASS_JACO, S("C_2"): 1.0, S("C_3"): 1.0, S("f_thick"): f_thick}
    for s, x in {"H": 0.5, "H_2": 0.25, "H+": 1e-6, "e-": 1e-6, "He": 0.094}.items():
        v[x_(s)], v[n_(s)] = x, x * n
    return v


_ROWS = {}


def dust_row(model, vals):
    """Td -> the model's dust heat row (the sum of every process's) at vals, rules applied"""
    if id(model) not in _ROWS:
        e = sp.sympify(sum(p.network["dust heat"].rhs for p in model.subprocesses if "dust heat" in p.network))
        syms = [do.T_dust] + sorted(e.free_symbols - {do.T_dust}, key=str)
        _ROWS[id(model)] = (syms, _c_lambdify([sanitize_symbols(s) for s in syms], [sanitize_symbols(e)], tables_in([e])))
    syms, f = _ROWS[id(model)]
    missing = set(syms[1:]) - set(vals)
    assert not missing, missing
    rest = [np.array([float(vals[s])]) for s in syms[1:]]
    return lambda t: float(f(np.array([t]), *rest)[0][0])


def model_dust_temperature(model, vals, frozen_Td):
    """The model's dust temperature at the gas temperature vals[T]: the root of its dust heat row where Td is solved,
    else the value the host packed at the start of the step"""
    if "Td" not in model.solve_vars:
        return frozen_Td
    return _walk_root(dust_row(model, vals), min(max(vals[T], 2.73), 1e4))


def _near_switch(Td):
    edges = (*do.ZONE_BOUNDARIES, 1500.0)
    return any(abs(Td - b) < 5.0 for b in edges)


# n_Htot, gas T, Sigma [g cm^-2], A [erg cm^-3 s^-1], Z: diffuse to first core
GRID = [(n, Tg, Sig, A, Z) for n in (1e2, 1e6, 1e10, 1e14) for Tg in (8.0, 40.0, 300.0, 1200.0, 2500.0)
        for Sig in (1e-3, 1.0, 30.0) for A in (1e-34, 1e-27) for Z in (1.0, 0.1)]


def _gizmo_args(n, Tg, Sig, A, Z, nH_convention="gizmo"):
    """dust_dEdt's arguments at a GRID state: GIZMO's nHcgs = 0.76 rho/m_p, or starforge's H density n_Htot, which its
    gas-dust collisions (and so the dust balance, which takes what the gas gives) use instead"""
    X = 0.7155
    return Tg, (X_GIZMO * n / X if nH_convention == "gizmo" else n), n * MP / X, Z, A * n, Sig


CONVENTION = {"starforge": "starforge", "starforge_legacy": "gizmo", "starforge_legacy_RT_EUV": "gizmo"}


@pytest.mark.parametrize("name", ["starforge", "starforge_legacy"])
def test_dust_temperature_is_gizmos_per_iterate_estimate(models, name):
    """The model's dust temperature at a gas temperature against GIZMO's root there. The start-of-step value the host
    packed (the root at 0.7 T, as a cell that has since been heated) is the frozen path, which this rejects where
    the gas-dust coupling sets Td"""
    model, conv = models[name], CONVENTION[name]
    for n, Tg, Sig, A, Z in GRID:
        ref = gizmo_dust_temperature(*_gizmo_args(n, Tg, Sig, A, Z, conv))
        frozen = gizmo_dust_temperature(*_gizmo_args(n, 0.7 * Tg, Sig, A, Z, conv))
        got = model_dust_temperature(model, values(Tg, n, Z, A * n, Sig), frozen)
        if _near_switch(ref) or _near_switch(got):
            continue
        assert got == pytest.approx(ref, rel=1e-8), (name, n, Tg, Sig, A, Z, got, ref)


@pytest.mark.parametrize("name,radtransfer", [("starforge_legacy", False), ("starforge_legacy_RT_EUV", True)])
def test_legacy_dust_balance_is_dust_dEdt(models, name, radtransfer):
    """starforge_legacy's dust heat row is dust_dEdt term for term (away from the smoothed switches), the coupling
    without CoolingRate's high-temperature truncations; under RADTRANSFER with the emission's albedo factor"""
    model = models[name]
    for n, Tg, Sig, A, Z in GRID + [(1e8, 5e4, 1.0, 1e-27, 1.0), (1e6, 4e5, 1.0, 1e-27, 1.0)]:
        row = dust_row(model, values(Tg, n, Z, A * n, Sig))
        for Td in (3.0, 20.0, 100.0, 600.0, 1400.0, 2000.0):
            if _near_switch(Td):
                continue
            ref = gizmo_dust_dEdt(Td, *_gizmo_args(n, Tg, Sig, A, Z), radtransfer=radtransfer)
            assert row(Td) == pytest.approx(ref, rel=1e-9, abs=1e-12 * abs(A * n)), (n, Tg, Sig, A, Z, Td)


def test_starforge_caps_the_dust_radiation_not_its_coupling(models):
    """starforge: the optically-thick limit's factor scales the dust's absorption and emission, not the gas-dust
    collisions, so the dust stays coupled where the limit binds; at f = 1 the row is GIZMO's below the truncations"""
    model = models["starforge"]
    for n, Tg, Sig, A, Z in [(1e14, 1500.0, 30.0, 1e-27, 1.0), (1e6, 40.0, 1.0, 1e-27, 1.0)]:
        args = _gizmo_args(n, Tg, Sig, A, Z, "starforge")
        for Td in (20.0, 300.0, 1400.0):
            full = dust_row(model, values(Tg, n, Z, A * n, Sig))(Td)
            capped = dust_row(model, values(Tg, n, Z, A * n, Sig, f_thick=1e-3))(Td)
            coupling = gizmo_dust_dEdt(Td, *args[:3], args[3], 0.0, 1e100)  # no emission (tau huge), no absorption
            radiation = gizmo_dust_dEdt(Td, *args) - coupling
            assert full == pytest.approx(coupling + radiation, rel=1e-9)
            assert capped == pytest.approx(coupling + 1e-3 * radiation, rel=1e-9)


def test_dust_temperature_is_a_steady_state_variable(models):
    """Td is a solve variable of every non-RT model, in steady state, on the dust heat row; no Td parameter remains"""
    for name, model in models.items():
        assert "Td" in model.solve_vars and "Td" not in model.time_dependent, name
        assert [v.row for v in model.variables if v.name == "Td"] == ["dust heat"], name
        assert "Td" not in {p.name for p in model.parameters}, name
