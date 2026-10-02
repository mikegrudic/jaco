"""Metal free electrons against a transcription of GIZMO's neutral-gas electron budget (gizmo_electrons.py).

Inputs are mapped as GIZMO's interface packs them: n_Htot is GIZMO's nHcgs, MolecularMassFraction = 2 x_H2,
column = m_p N_H, x_X = Z_X / (A_X X); X is solar so that x_X / x_solar(X) = Z_X / Z_X,sun exactly."""

import numpy as np
import pytest
import sympy as sp
from jaco.symbols import x_
from ..metal_electrons import (metal_electrons, heavy_ion_electrons, alkali_electrons, Oplus_electrons,
                               molecular_ion_electrons, G0_carbon, f_Cplus, X_H)
from ..symbols import T, n_Htot, G_0, NH, Z_dust, f_dust, ISRF, x_solar
from .c_semantics import c_lambdify
from . import gizmo_electrons as G

XSUN = 1 - G.SOLAR[0] - G.SOLAR[1]
x_C_tot, x_O_tot, x_e = sp.Symbol("x_C,tot"), sp.Symbol("x_O,tot"), sp.Symbol("x_e")
ARGS = (T, n_Htot, G_0, NH, Z_dust, f_dust, ISRF, X_H, x_C_tot, x_O_tot, x_("Fe"), x_("H+"), x_("He+"), x_("He++"),
        x_("H_2"))


def state(Tv, nH, G0, column, fmol, xHp, Zsol=1.0, isrf=1.0):
    """(jaco arguments, GIZMO inputs) for one cell; He ions at 1e-3 and 1e-6 of x_H+"""
    Z = [Zsol * z for z in G.SOLAR]
    Z[1] = G.SOLAR[1]
    xHep, xHepp = 1e-3 * xHp, 1e-6 * xHp
    args = (Tv, nH, G0, column / G.PROTONMASS_CGS, Zsol, 1.0, isrf, XSUN, Z[2] / (12 * XSUN), Z[4] / (16 * XSUN),
            Z[10] / (56 * XSUN), xHp, xHep, xHepp, 0.5 * fmol)
    gz = dict(temp=Tv, nHcgs=nH, density_cgs=nH * G.PROTONMASS_CGS / XSUN, n_elec_HHe=xHp + xHep + 2 * xHepp,
              nHp=xHp, G0=G0, column=column, fmol=fmol, Z=Z, isrf=isrf)
    return args, gz


def evaluate(expr, args):
    return float(c_lambdify(ARGS, expr)(*args))


def zeta_of(gz):
    return G.zeta_cr(gz["column"], gz["Z"], gz["isrf"])


# diffuse atomic, CNM, transition, shielded molecular, dense core, protostellar, warm ionized, hot dense
STATES = [state(*s) for s in [
    (6000.0, 0.3, 1.7, 1e-5, 0.0, 3e-3), (100.0, 20.0, 1.0, 1e-3, 0.05, 1e-5), (40.0, 80.0, 0.3, 1e-2, 0.3, 1e-6),
    (15.0, 500.0, 0.05, 0.05, 0.9, 1e-8), (10.0, 1e5, 1e-6, 1.0, 1.0, 1e-10), (30.0, 1e9, 0.0, 30.0, 1.0, 1e-14),
    (8000.0, 0.1, 1.7, 1e-6, 0.0, 0.05), (1500.0, 1e8, 1e-3, 10.0, 0.5, 1e-9), (200.0, 3.0, 0.5, 3e-3, 0.01, 2e-4, 0.1, 10.0),
]]


@pytest.mark.parametrize("args,gz", STATES)
def test_closed_form_terms(args, gz):
    zeta = zeta_of(gz)
    assert evaluate(molecular_ion_electrons(), args) == pytest.approx(G.molecular_ions(gz["temp"], gz["nHcgs"], zeta), rel=1e-4)
    assert evaluate(alkali_electrons(), args) == pytest.approx(G.alkali(gz["temp"], gz["nHcgs"], gz["Z"]), rel=1e-10, abs=1e-300)
    assert evaluate(Oplus_electrons(), args) == pytest.approx(G.Oplus(gz["nHp"], gz["Z"]), rel=1e-10, abs=1e-300)
    G0c = G.get_FUV_G0_mode1(gz["G0"], gz["column"], gz["Z"], gz["fmol"])
    assert evaluate(G0_carbon(), args) == pytest.approx(G0c, rel=1e-10, abs=1e-300)
    for xe in (1e-7, 1e-4, 1e-2):
        expected = G.f_Cplus(gz["temp"], xe, gz["nHcgs"], gz["G0"], G0c, gz["fmol"], zeta, gz["Z"])
        assert evaluate(f_Cplus(x_e).subs(x_e, xe), args) == pytest.approx(expected, rel=1e-4)


# every branch of return_electron_fraction_from_heavy_ions: (T, n_H, column, x_H+) -> regime. Regime 4 (alpha > 10
# above 100 n_crit) needs T < 0.03 K.
HEAVY = [((100.0, 1.0, 1e-5, 0.02), 1), ((10.0, 1.0, 1e-5, 1e-6), 2), ((10.0, 1e6, 1.0, 1e-12), 3),
         ((0.003, 1e13, 1.0, 1e-14), 4), ((10.0, 1e13, 1.0, 1e-14), 5), ((10.0, 1e11, 1.0, 1e-14), 6),
         ((10.0, 1e10, 1.0, 1e-14), 7)]


@pytest.mark.parametrize("cell,regime", HEAVY)
def test_heavy_ions_every_regime(cell, regime):
    Tv, nH, column, xHp = cell
    args, gz = state(Tv, nH, 0.0, column, 1.0, xHp)
    expected, taken = G.heavy_ions(gz["temp"], gz["density_cgs"], gz["n_elec_HHe"], zeta_of(gz), gz["Z"])
    assert taken == regime
    assert evaluate(heavy_ion_electrons(), args) == pytest.approx(expected, rel=1e-4)


def evaluate_metal_electrons(args):
    inter, total = metal_electrons()
    vals = {}
    for w, W in inter:
        vals[w] = float(c_lambdify(ARGS + tuple(vals), W)(*args, *vals.values()))
    return float(c_lambdify(ARGS + tuple(vals), total)(*args, *vals.values())), vals


@pytest.mark.parametrize("args,gz", STATES)
def test_total_matches_converged_gizmo_loop(args, gz):
    """Two Newton steps reach GIZMO's converged n_e within its own 1% tolerance"""
    expected, _ = G.metal_electrons(**gz)
    total, _ = evaluate_metal_electrons(args)
    assert total == pytest.approx(expected, rel=1e-2)


def test_Cplus_plateau():
    """Photoionized C+ in diffuse neutral gas: GIZMO's 1.6e-4 x Z_C"""
    args, gz = STATES[1]
    _, vals = evaluate_metal_electrons(args)
    assert vals[sp.Symbol("eCplus")] == pytest.approx(1.6e-4, rel=0.02)


@pytest.mark.parametrize("Tv,lg_lo,lg_hi,xHp", [(30.0, 4.4, 4.9, 1e-12),   # gas-phase -> interpolated
                                               (100.0, 7.9, 8.4, 1e-14),  # interpolated -> dust (35% jump in GIZMO)
                                               (100.0, 0.0, 0.0, None)])  # ionized cut at x_e(H, He) = 0.01
def test_heavy_ions_continuous_across_switches(Tv, lg_lo, lg_hi, xHp):
    """GIZMO's regime switches are blended over ~0.01 dex: on a 0.001 dex grid no step carries more than a tenth of the
    change across the window"""
    f = c_lambdify(ARGS, heavy_ion_electrons())
    if xHp is None:  # scan x_H+ across 0.01 at n_H = 1
        grid = 10 ** np.arange(-2.05, -1.95, 1e-3)
        vals = np.array([f(*state(Tv, 1.0, 0.0, 1e-5, 0.0, x / 1.001002)[0]) for x in grid])
    else:
        grid = 10 ** np.arange(lg_lo, lg_hi, 1e-3)
        vals = np.array([f(*state(Tv, n, 0.0, 1.0, 1.0, xHp)[0]) for n in grid])
    assert np.all(np.isfinite(vals))
    steps = np.abs(np.diff(np.log(vals)))
    assert steps.max() < 0.1 * steps.sum()
