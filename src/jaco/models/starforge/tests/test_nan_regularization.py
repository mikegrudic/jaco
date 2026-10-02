"""H2 collisional dissociation stays finite, with finite partials, at the states where the generated C code used to
produce nan (cold gas; hot gas), and is unchanged elsewhere. Evaluated with C floating-point
semantics (see c_semantics.py); the derivatives are taken in the solve variables after the conservation reductions."""

import numpy as np
import pytest
import sympy as sp
from jaco.symbols import x_
from ..h2_chemistry.collisional_dissociation import H2_collisional_dissociation
from ..symbols import T, n_Htot, NH, grad_v, dx
from .c_semantics import c_lambdify

y = sp.Symbol("y")
x_Hp, x_H2, x_Hep, x_Hepp = x_("H+"), x_("H_2"), x_("He+"), x_("He++")
REDUCTIONS = {x_("H"): 1 - x_Hp - 2 * x_H2, x_("He"): y - x_Hep - x_Hepp}
SOLVE_VARS = (T, x_H2, x_Hp)
ARGS = (T, n_Htot, x_Hp, x_H2, x_Hep, x_Hepp, y, NH, grad_v, dx)
DISSOCIATIONS = [("H_2", c) for c in ("H+", "e-", "H", "H_2", "He")] + [("HD", "e-"), ("D_2", "e-")]


def value_and_partials(expr):
    expr = expr.subs(REDUCTIONS)
    return c_lambdify(ARGS, expr), [c_lambdify(ARGS, sp.diff(expr, v)) for v in SOLVE_VARS]


def state(Tv, nH, xHp, xH2, NHv=1e21):
    """ARGS values; He mostly neutral, column and velocity gradient typical of a GMC cell."""
    return (Tv, nH, xHp, xH2, 1e-4 * xHp, 1e-5 * xHp, 0.0944, NHv, 1e-13, 3e17)


PATHOLOGICAL = {
    "cold molecular": state(10.0, 1e4, 1e-7, 0.49),
    "very cold": state(3.0, 1e6, 1e-9, 0.4999),
    "fully ionized, no H2": state(2e4, 1.0, 1.0, 0.0),
    "fully ionized, trace H2": state(1e6, 1e-2, 1.0, 1e-15),
    "hot": state(1e8, 1e-3, 1.0, 0.0),
}
REGULAR = [state(*s) for s in [(100.0, 1e2, 1e-4, 0.3), (300.0, 1e4, 1e-6, 0.45), (1e3, 1.0, 1e-3, 0.1),
                               (3e3, 1e6, 1e-2, 0.2), (1e4, 1e2, 0.5, 1e-3), (3e4, 1.0, 0.9, 1e-6), (1e5, 1e8, 0.1, 1e-2)]]


def old_dissociation_rate(iso, collider, Tv, nH, xHp, xH2, xHep, xHepp, yv, *_):
    """The previous k_0^f_0 k_LTE^(1-f_0) form (as in GIZMO), transcribed directly."""
    lgT, lnT = np.log10(Tv), np.log(Tv)
    xH, xHe = 1 - xHp - 2 * xH2, yv - xHep - xHepp
    lt4 = lgT - 4.0
    ncrit = 1.0 / (xH / 10 ** (3.0 - 0.416 * lt4 - 0.327 * lt4**2) + xH2 / 10 ** (4.845 - 1.3 * lt4 + 1.62 * lt4**2)
                   + xHe / 10 ** (5.0792 * (1.0 - 1.23e-5 * (Tv - 2000.0))))
    f0 = 1.0 / (1.0 + nH / (100 * ncrit if iso == "HD" else ncrit))
    if collider == "H+":
        k0 = kLTE = (-3.3232183e-7 + 3.3735382e-7 * lnT - 1.4491368e-7 * lnT**2 + 3.4172805e-8 * lnT**3
                     - 4.7813720e-9 * lnT**4 + 3.9731542e-10 * lnT**5 - 1.8171411e-11 * lnT**6
                     + 3.5311932e-13 * lnT**7) * np.exp(-21237.15 / Tv)
    elif collider == "e-":
        a0, b0, c0, aL, bL, cL = {"H_2": (4.49e-9, 0.11, 101858.0, 1.91e-9, 0.136, 53407.1),
                                  "D_2": (8.24e-9, 0.216, 105388, 1.91e-9, 0.163, 53339.7),
                                  "HD": (5.09e-9, 0.128, 103258, 1.04e-9, 0.218, 53070.7)}[iso]
        k0, kLTE = a0 * Tv**b0 * np.exp(-c0 / Tv), aL * Tv**bL * np.exp(-cL / Tv)
    elif collider == "H":
        k0, kLTE = 6.67e-12 * np.sqrt(Tv) * np.exp(-(1.0 + 63593.0 / Tv)), 3.52e-9 * np.exp(-43900.0 / Tv)
    elif collider == "H_2":
        k0 = 5.996e-30 * Tv**4.1881 * np.exp(-54657.4 / Tv) / (1.0 + 6.761e-6 * Tv) ** 5.6881
        kLTE = 1.3e-9 * np.exp(-53300.0 / Tv)
    else:
        k0, kLTE = 10 ** (-27.029 + 3.801 * lgT - 29487.0 / Tv), 10 ** (-2.729 - 1.75 * lgT - 23474.0 / Tv)
    return max(k0, 0) ** f0 * max(kLTE, 0) ** (1 - f0)


@pytest.mark.parametrize("iso,collider", DISSOCIATIONS)
@pytest.mark.parametrize("point", list(PATHOLOGICAL))
def test_dissociation_finite(iso, collider, point):
    f, partials = value_and_partials(H2_collisional_dissociation(collider, iso).rate_coefficient)
    vals = PATHOLOGICAL[point]
    assert np.isfinite(f(*vals)) and f(*vals) >= 0
    for d, v in zip(partials, SOLVE_VARS):
        assert np.isfinite(d(*vals)), f"d k / d {v} not finite"


@pytest.mark.parametrize("iso,collider", DISSOCIATIONS)
def test_dissociation_unchanged_in_regular_regime(iso, collider):
    f, _ = value_and_partials(H2_collisional_dissociation(collider, iso).rate_coefficient)
    checked = 0
    for vals in REGULAR:
        expected = old_dissociation_rate(iso, collider, *vals)
        if expected < 1e-200:  # below the generated code's exp(-500) clamp; neither form is meaningful there
            continue
        assert f(*vals) == pytest.approx(expected, rel=1e-10, abs=0)
        checked += 1
    assert checked >= 5
