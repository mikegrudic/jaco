"""Sub-grid clumping: STARFORGE multiplies every two-body rate by C_2 and every three-body rate by C_3; STARFORGE_LEGACY
(rule GIZMO_CLUMPING) clumps only the H2 terms GIZMO's update_explicit_molecular_fraction multiplies by its
clumping_factor.

Walks every process of each model: each rate and heat expression must scale as C_2^1 (C_3^1) if its process is two-body
(three-body) and not at all if it is one-body, at a sub-critical state. New processes must be classified here."""

import numpy as np
import pytest
import sympy as sp
from jaco.symbols import n_, x_
from ..starforge import make_model
from ...starforge_legacy import make_model as make_legacy
from ..ionization_balance import ion_abundances
from ..symbols import clumping_factor, T, grad_v, grad_v_tf, dx, x_solar

C2, C3 = sp.Symbol("C_2"), sp.Symbol("C_3")
ONE_BODY = {"Cosmic ray heating", "Photoelectric Heating", "Inverse Compton cooling (CMB)", "Photodissociation of H_2",
            "Photodetachment of H-", "Dissociation of H_2 by cosmic rays", "Direct ionization of H by cosmic rays",
            "PdV work"}
THREE_BODY = {"3-body formation of H_2"}

ABUNDANCES = {"H": 0.6, "H+": 1e-4, "H_2": 0.2, "He": 0.094, "He+": 1e-8, "He++": 1e-12, "e-": 2e-4, "H-": 1e-10,
              "H_2+": 0.0, "HD": 1e-5, "C+": 1e-4, "CO": 5e-5, "C": 1.5e-4, "O": 5e-4, "N": 7e-5, "Ne": 9e-5, "Mg": 4e-5,
              "Si": 3e-5, "S": 1e-5, "Ca": 2e-6, "Fe": 3e-5}
PARAMS = {"T": 100.0, "n_Htot": 0.1, "G_0": 1.0, "G_LW": 1.0, "N_H": 1e21, "∇v": 1e-14, "Δx": 1e18, "Td": 15.0, "X": 0.7155,
          "y": 0.094, "z": 0.0, "f_metal": 1.0, "f_neb": 1.0, "Z_d": 1.0, "f_d": 1.0, "ISRF": 1.0, "x_C,tot": x_solar("C"),
          "x_O,tot": x_solar("O"), "xIbMgplus": 1e-5, "pdv_work": 1e-25}


def evaluate(expr, c2=1.0, c3=1.0):
    vals = {sp.Symbol(k): v for k, v in PARAMS.items()}
    for s, x in ABUNDANCES.items():
        vals[x_(s)] = x
        vals[n_(s)] = x * PARAMS["n_Htot"]
    vals.update({C2: c2, C3: c3})
    missing = expr.free_symbols - set(vals)
    assert not missing, missing
    return float(sp.N(expr.xreplace(vals)))


def clumping_exponents(expr):
    f0 = evaluate(expr)
    if f0 == 0:
        return None
    return (np.log(evaluate(expr, c2=1.01) / f0) / np.log(1.01), np.log(evaluate(expr, c3=1.01) / f0) / np.log(1.01))


def expressions(process):
    return {k: eq.rhs for k, eq in process.network.items() if eq.rhs != 0}


@pytest.mark.parametrize("process", make_model().subprocesses, ids=lambda p: p.name)
def test_starforge_clumps_every_collision(process):
    expected = (0, 0) if process.name in ONE_BODY else (0, 1) if process.name in THREE_BODY else (1, 0)
    for key, expr in expressions(process).items():
        exps = clumping_exponents(expr)
        if exps is not None:
            assert exps == pytest.approx(expected, abs=0.02), key


@pytest.mark.parametrize("process", make_legacy().subprocesses, ids=lambda p: p.name)
def test_legacy_clumps_only_H2_chemistry(process):
    for key, expr in expressions(process).items():
        if process.name == "GIZMO H2 network":
            assert {C2, C3} <= expr.free_symbols, key
        else:
            assert not {C2, C3} & expr.free_symbols, key


def by_name(model, name):
    return next(p for p in model.subprocesses if p.name == name)


def test_cooling_and_H2_terms_each_model():
    sf, legacy = make_model(), make_legacy()
    assert clumping_exponents(by_name(sf, "Gas-dust collisions").network["heat"].rhs)[0] == pytest.approx(1)
    assert clumping_exponents(by_name(sf, "Formation of H_2 on dust grains").network["H_2"].rhs)[0] == pytest.approx(1)
    assert C2 not in by_name(legacy, "Gas-dust collisions").network["heat"].rhs.free_symbols
    # with no H2 yet the GIZMO network is all formation, which GIZMO clumps (a_Z, a_GP by C_2, b_3B by C_3)
    h2 = by_name(legacy, "GIZMO H2 network").network["H_2"].rhs
    saved = ABUNDANCES["H_2"]
    ABUNDANCES["H_2"] = 1e-12
    try:
        assert clumping_exponents(h2)[0] == pytest.approx(1, abs=0.01)
    finally:
        ABUNDANCES["H_2"] = saved


def test_starforge_ionization_balance_clumps_recombination():
    """C+, Mg+ and molecular ions at fixed x_e: their two-body sinks scale with C_2, so clumping lowers them"""
    x_e = sp.Symbol("x_e")
    abundances = ion_abundances(x_e)
    lo = [evaluate(a.subs(x_e, 2e-4), c2=1.0) for a in abundances]
    hi = [evaluate(a.subs(x_e, 2e-4), c2=4.0) for a in abundances]
    assert all(h < l for h, l in zip(hi, lo))


def test_clumping_estimator_is_gizmos():
    """1 + (b dv / c_s)^2, b = 0.5, dv = |grad v| dx, c_s = v_th,rms / sqrt(3) with v_th,rms = 0.111 sqrt(T) km/s"""
    Tv, gv, dxv = 50.0, 2e-14, 3e18
    dv_kms = gv * dxv / 1e5
    expected = 1 + (0.5 * dv_kms / (0.111 * np.sqrt(Tv) / np.sqrt(3))) ** 2
    assert float(clumping_factor.subs({T: Tv, grad_v: gv, dx: dxv})) == pytest.approx(expected, rel=1e-12)


def test_clumping_gradient_each_model():
    """STARFORGE estimates C_2 from the trace-free velocity gradient, so homologous flow is not sub-grid turbulence;
    STARFORGE_LEGACY keeps GIZMO's full norm"""
    sf, legacy = make_model().derived["C_2"].free_symbols, make_legacy().derived["C_2"].free_symbols
    assert grad_v_tf in sf and grad_v not in sf
    assert grad_v in legacy and grad_v_tf not in legacy
