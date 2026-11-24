"""STARFORGE's explicit ionization balance (ionization_balance.py): the Newton steps of the C+/Mg+/molecular-ion fixed point,
charge balance in the cosmic-ray-dominated regime, and that regime against GIZMO's heavy-ion term.

GIZMO's return_electron_fraction_from_heavy_ions is its CR ionization of the neutrals, with the charge on Mg+ but
recombining at k_ei = 9.77e-8 cm^3 s^-1 (a molecular-ion rate) and captured by 0.1 um grains only. STARFORGE's Mg+
recombines radiatively (2.8e-11 at 10 K) and on grains at the WD01 rate (~2e-14 n_H per ion, small grains and PAHs
included, ~500x GIZMO's grain capture), and molecular ions recombine at Fromang's 3e-6/sqrt(T). With charge transfer
H+ (and molecular ions) -> Mg the two are within 2x at n_H ~ 1e2 cm^-3; denser, STARFORGE's grains take over (x_e ~ zeta/n
against GIZMO's sqrt(zeta/n)) and it falls below GIZMO, by 4x at 1e3 and 13x at 1e4 cm^-3 for these states."""

import itertools
import numpy as np
import pytest
import sympy as sp
from ..ionization_balance import (ion_abundances, solved_electrons, K_CT_HPLUS_MG, alpha_rr_Mgplus, beta_molion, x_Mg,
                                  Oplus_electrons)
from ..metal_electrons import x_e_ions, alkali_electrons, metal_electrons
from ..grain_assisted_recombination import alpha_grain
from ..symbols import T, n_Htot, G_0, NH, Z_dust, f_dust, ISRF, x_, cosmicray_ionization_rate_H as zeta
from jaco.processes.recombination import hydrogenic_recombination_rate
from .c_semantics import c_lambdify
from . import gizmo_electrons as G

X, C2, xe = sp.Symbol("X"), sp.Symbol("C_2"), sp.Symbol("xe")
xC, xO, xFe = sp.Symbol("x_C,tot"), sp.Symbol("x_O,tot"), sp.Symbol("x_Fe")
ARGS = (T, n_Htot, G_0, NH, Z_dust, f_dust, ISRF, X, xC, xO, x_Mg, xFe, x_("H+"), x_("He+"), x_("He++"), x_("H_2"), C2)
SOLAR_X = 1 - G.SOLAR[0] - G.SOLAR[1]
METALS = (G.SOLAR[2] / 12 / SOLAR_X, G.SOLAR[4] / 16 / SOLAR_X, G.SOLAR[6] / 24 / SOLAR_X, G.SOLAR[10] / 56 / SOLAR_X)
F = c_lambdify(ARGS + (xe,), sum(ion_abundances(xe)))
other = c_lambdify(ARGS, alkali_electrons() + Oplus_electrons())
ions = c_lambdify(ARGS, x_e_ions)
INTER, TOTAL, _ = solved_electrons()
FUNCS = [c_lambdify(ARGS + tuple(w for w, _ in INTER[:k]), W) for k, (w, W) in enumerate(INTER)]
FTOTAL = c_lambdify(ARGS + tuple(w for w, _ in INTER), TOTAL)


def args_of(Tv, nH, G0, column, fmol, xHp, C2v=1.0):
    xH2 = 0.5 * fmol * (1 - xHp)
    return (Tv, nH, G0, column / G.PROTONMASS_CGS, 1.0, 1.0, 1.0, SOLAR_X, *METALS, xHp, 1e-3 * xHp, 1e-6 * xHp, xH2, C2v)


def newton_electrons(a):
    vals = []
    for f in FUNCS:
        vals.append(float(f(*a, *vals)))
    return float(FTOTAL(*a, *vals))


def root_y(a, x_o):
    """bisection on ln F(y) = ln y"""
    lo, hi = -60.0, 0.0
    for _ in range(200):
        m = 0.5 * (lo + hi)
        lo, hi = (m, hi) if np.log(float(F(*a, x_o + np.exp(m)))) > m else (lo, m)
    return np.exp(0.5 * (lo + hi))


GRID = list(itertools.product([10.0, 100.0, 8e3], [0.3, 1e3, 1e8], [1.7, 1e-4, 0.0], [1e-4, 1.0], [0.0, 1.0], [1e-12, 1e-3],
                              [1.0, 10.0]))


@pytest.mark.parametrize("Tv,nH,G0,column,fmol,xHp,C2v", GRID[::7])
def test_newton_steps_reach_the_root(Tv, nH, G0, column, fmol, xHp, C2v):
    a = args_of(Tv, nH, G0, column, fmol, xHp, C2v)
    x_o = float(ions(*a)) + float(other(*a))
    y = newton_electrons(a) - float(other(*a))
    assert y == pytest.approx(root_y(a, x_o), rel=1e-2)


ZETA = c_lambdify(ARGS, zeta)
H_SINK = c_lambdify(ARGS + (xe,), hydrogenic_recombination_rate(1) * xe + alpha_grain("H+", xe)
                    + K_CT_HPLUS_MG * (x_Mg - ion_abundances(xe)[1]))  # per H+ and n_H


def cr_dominated_state(Tv, nH, column, fmol):
    """x_H+ from CR ionization of H against its radiative, grain and charge-transfer sinks, with the other free
    electrons solved alongside; no FUV, clumping 1"""
    def at(xHp):
        a = args_of(Tv, nH, 0.0, column, fmol, xHp)
        x_o = float(ions(*a)) + float(other(*a))
        return a, x_o + root_y(a, x_o)
    lo, hi = -45.0, np.log(0.1)
    for _ in range(100):
        m = 0.5 * (lo + hi)
        a, x_e = at(np.exp(m))
        xHI = 1 - np.exp(m) - fmol * (1 - np.exp(m))
        lo, hi = (m, hi) if float(ZETA(*a)) * xHI > np.exp(m) * nH * float(H_SINK(*a, x_e)) else (lo, m)
    return at(np.exp(0.5 * (lo + hi)))


DENSE = [(20.0, 1e2, 0.1, 0.5), (15.0, 1e3, 0.3, 1.0), (10.0, 1e4, 1.0, 1.0)]


@pytest.mark.parametrize("Tv,nH,column,fmol", DENSE)
def test_cr_ionization_balanced_by_recombination(Tv, nH, column, fmol):
    """Every CR ionization of H or H2 ends in a recombination: radiative or grain (H+, Mg+), dissociative (molecular ions)"""
    a, x_e = cr_dominated_state(Tv, nH, column, fmol)
    xHp, xH2 = a[12], a[15]
    xC, xMg, xmol = [float(c_lambdify(ARGS + (xe,), v)(*a, x_e)) for v in ion_abundances(xe)]
    rate = lambda e: float(c_lambdify(ARGS + (xe,), e)(*a, x_e))
    production = rate(zeta) * (1 - xHp - 2 * xH2 + 2 * xH2)
    recombination = nH * (xHp * rate(hydrogenic_recombination_rate(1) * xe + alpha_grain("H+", xe))
                          + xMg * rate(alpha_rr_Mgplus * xe + alpha_grain("Mg+", xe)) + xmol * rate(beta_molion * xe))
    assert xC < 0.05 * x_e  # carbon is neither photo- nor much CR-ionized here
    assert recombination == pytest.approx(production, rel=0.05)


@pytest.mark.parametrize("Tv,nH,column,fmol,ratio_lo,ratio_hi", [(20.0, 1e2, 0.1, 0.5, 0.3, 0.7),
                                                                 (15.0, 1e3, 0.3, 1.0, 0.15, 0.4),
                                                                 (10.0, 1e4, 1.0, 1.0, 0.04, 0.12)])
def test_cr_dominated_electrons_vs_gizmo(Tv, nH, column, fmol, ratio_lo, ratio_hi):
    """GIZMO is in its gas-phase-recombination regime (II) here; the ratio falls as the grains take over in STARFORGE"""
    a, x_e = cr_dominated_state(Tv, nH, column, fmol)
    Z = list(G.SOLAR)
    density = nH * G.PROTONMASS_CGS / SOLAR_X
    _, regime = G.heavy_ions(Tv, density, 0.0, G.zeta_cr(column, Z), Z)
    assert regime == 2
    legacy_args = list(args_of(Tv, nH, 0.0, column, fmol, 1e-14))
    inter, total = metal_electrons()
    vals = []
    for k, (w, W) in enumerate(inter):
        vals.append(float(c_lambdify(ARGS + tuple(v for v, _ in inter[:k]), W)(*legacy_args, *vals)))
    legacy = float(c_lambdify(ARGS + tuple(v for v, _ in inter), total)(*legacy_args, *vals))
    assert ratio_lo < x_e / legacy < ratio_hi
