"""starforge_legacy_RT against direct transcriptions of GIZMO's legacy RT coupling (cooling/cooling.cc, rt_chem.cc,
rt_utilities.cc at gizmo_jaco_dev bd3500d7)."""

import numpy as np
import pytest
import sympy as sp
from jaco.symbols import n_, x_, dt
from jaco.processes.recombination import gasphase_recombination_rates
from jaco.processes.ionization import collisional_ionization_rates
from ..symbols import T, n_Htot, X_H, z
from ...starforge_legacy import make_model as make_legacy
from ...starforge_legacy_RT import (make_model, slab_average, GIZMO_ABUNDANCES, GIZMO_EMISSION_BACKGROUND, Gamma_HI,
                                    sigma_HI, eps_HI, c_tilde, T_bg)

STATE = {T: 8.0e3, n_Htot: 100.0, X_H: 0.7, z: 0.0}
F_IR = sp.Symbol("f_IR_selfabs")


@pytest.fixture(scope="module")
def models():
    return make_model(), make_legacy()


def process(model, name):
    return next(p for p in model.subprocesses if p.name == name)


def per_H(expr, xs):
    """expr with number densities n_s = n_Htot x_s, at STATE"""
    rep = {n_(s): n_Htot * x for s, x in xs.items()}
    return float(expr.xreplace(rep).subs({**STATE, Gamma_HI: 1e-9, sigma_HI: 3e-18, eps_HI: 4.8e-12,
                                          c_tilde: 0.0, T_bg: 2.73, dt: 1e11, F_IR: 1.0}))


def gizmo_dxHp_dt(Tv, xHp, xe, nHcgs, Gamma):
    """d x_H+/dt of find_abundances_and_rates' H balance (cooling.cc 900-904) at its converged n_e: the backward-Euler
    update x_H0 = (HI + f a)/(1 + f (a + ge + Gamma/n_e)), f = dt n_e, read as a rate"""
    aHp = float(gasphase_recombination_rates["H+"].subs(T, Tv))
    geH0 = float(collisional_ionization_rates["H"].subs(T, Tv))
    necgs = xe * nHcgs
    return Gamma * (1 - xHp) + necgs * geH0 * (1 - xHp) - necgs * aHp * xHp


@pytest.mark.parametrize("xHp", [1e-4, 0.3, 0.99])
def test_H_balance_runs_at_gizmo_speed(models, xHp):
    """The H+ row per n_Htot is GIZMO's dx_H+/dt: rates at nHcgs = (0.76/X) n_Htot, abundances per nHcgs"""
    rt, _ = models
    xe = xHp + 2e-4  # metal electrons on top
    xs = {"H": 1 - xHp, "H+": xHp, "e-": xe}
    rows = sum(process(rt, name).network["H+"].rhs for name in
               ("Photoionization of H by the RT band", "Gas-phase recombination of H+", "Collisional Ionization of H"))
    nHcgs = 0.76 / STATE[X_H] * STATE[n_Htot]
    assert per_H(rows, xs) / STATE[n_Htot] == pytest.approx(gizmo_dxHp_dt(STATE[T], xHp, xe, nHcgs, 1e-9), rel=1e-10)


def test_photoheating_is_heat_ion_from_RHD(models):
    """Heat_Ion_from_RHD (cooling.cc 1229-1263) times nHcgs^2: rt_ion_G_HI sigma c n_gamma nH0 / nHcgs per nHcgs^2"""
    rt, _ = models
    xH0 = 0.4
    heat = per_H(process(rt, "Photoionization of H by the RT band").heat, {"H": xH0})
    nHcgs = 0.76 / STATE[X_H] * STATE[n_Htot]
    assert heat == pytest.approx(4.8e-12 * 1e-9 * xH0 * nHcgs, rel=1e-12)


def test_frozen_law_at_zero_c_tilde_and_photon_budget():
    """S(0) = 1 exactly (GIZMO's law); GIZMO's Pade form of S = (1 - e^-x)/x is within 0.24% of it (worst near x = 5);
    with the time-averaged law the photons a cell takes over the step, (c_tilde/c) Gamma_eff n_H0 dt = n_gamma tau S(tau)
    ~ n_gamma (1 - e^-tau), exceed the band's n_gamma by at most the Pade error, 0.08%"""
    assert float(slab_average(sp.S.Zero)) == 1.0
    for x in np.logspace(-4, 4, 41):
        assert float(slab_average(sp.Float(x))) == pytest.approx(-np.expm1(-x) / x, rel=2.4e-3)
    for tau in np.logspace(-3, 6, 37):
        taken_per_gamma = tau * float(slab_average(sp.Float(tau)))
        assert taken_per_gamma <= 1.0008 and taken_per_gamma == pytest.approx(-np.expm1(-tau), rel=2.4e-3)


def test_slab_optical_depth_at_gizmo_density(models):
    """tau = c_tilde sigma n_H0 dt at nHcgs, GIZMO's kappa rho = sigma HI 0.76 rho/m_p (rt_kappa)"""
    rt, _ = models
    rate = process(rt, "Photoionization of H by the RT band").network["H+"].rhs
    xH0, ct, d = 0.5, 3e6, 1e11
    nHcgs = 0.76 / STATE[X_H] * STATE[n_Htot]
    tau = ct * 3e-18 * nHcgs * xH0 * d
    got = float(rate.xreplace({n_("H"): n_Htot * xH0}).subs({**STATE, Gamma_HI: 1e-9, sigma_HI: 3e-18, c_tilde: ct,
                                                               dt: d}))
    # the species row is per nHcgs: Gamma_eff x_H0 n_Htot
    assert got == pytest.approx(1e-9 * float(slab_average(sp.Float(tau))) * xH0 * STATE[n_Htot], rel=1e-12)


def test_emission_corrections_take_T_bg(models):
    """get_background_radiation_temperature_for_emission_corrections replaces the CMB in the bath factor of metal-line,
    molecular and fine-structure cooling; Compton keeps the CMB"""
    rt, legacy = models
    for name in ("Metal line cooling", "GIZMO H2 + HD cooling", "GIZMO C+, [CI] and CO cooling"):
        h_rt, h_lg = process(rt, name).heat.subs(F_IR, 1), process(legacy, name).heat
        assert z not in h_rt.free_symbols and T_bg in h_rt.free_symbols
        assert sp.simplify(h_rt.subs(T_bg, 2.73 * (1 + z)) - h_lg) == 0
    compton = process(rt, "Inverse Compton cooling (CMB)").heat
    assert T_bg not in compton.free_symbols and compton.subs(F_IR, 1) == process(legacy, "Inverse Compton cooling (CMB)").heat


def test_species_rows_per_nHcgs_heat_unchanged(models):
    """Rule GIZMO_ABUNDANCES divides species rows (not heat) by nHcgs/n_Htot; the GIZMO H2 network is exempt"""
    rt, legacy = models
    lam = 0.76 / X_H
    rec_rt, rec_lg = process(rt, "Gas-phase recombination of H+"), process(legacy, "Gas-phase recombination of H+")
    assert sp.simplify(rec_rt.network["H+"].rhs - rec_lg.network["H+"].rhs / lam) == 0
    assert sp.simplify(rec_rt.heat.subs(F_IR, 1) - rec_lg.heat) == 0
    h2_rt, h2_lg = process(rt, "GIZMO H2 network"), process(legacy, "GIZMO H2 network")
    assert sp.simplify(h2_rt.network["H_2"].rhs - h2_lg.network["H_2"].rhs) == 0
    assert "GIZMO H2 network" in GIZMO_ABUNDANCES.exempt
    assert GIZMO_EMISSION_BACKGROUND.exempt == {"Inverse Compton cooling (CMB)"}


def test_declarations(models):
    rt, legacy = models
    assert rt.solve_vars == legacy.solve_vars
    assert rt.time_dependent == ("T", "H+", "H_2")
    assert {"Gamma_HI", "sigma_HI", "eps_HI", "c_tilde", "T_bg"} <= {p.name for p in rt.parameters}
    assert set(rt.processes) - set(legacy.processes) == {"Photoionization of H by the RT band"}
    assert rt.rules[:2] == legacy.rules


@pytest.mark.slow
def test_check_with_an_ionizing_field():
    import jaco
    from jaco.model_check import CheckGrid
    grid = CheckGrid(np.logspace(np.log10(3.0), 9, 6), np.logspace(-3, 9, 3), n_random=2)
    field = dict(Gamma_HI=1e-8, sigma_HI=3e-18, eps_HI=4.8e-12, c_tilde=3e6, T_bg=20.0, f_IR_selfabs=0.7, f_recNUV=0.5)
    report = jaco.check(make_model(), grid, params=field, name="starforge_legacy_RT")
    assert report.ok, report.summary()


# --- B2: the outputs GIZMO's cooling-radiation return reads (CoolingRate 1317-1348) ---

NUV_TERMS = ["Metal line cooling", "Nebular forbidden-line cooling", "H-e- Line Cooling", "He+-e- Line Cooling"]
IR_TERMS = ["GIZMO C+, [CI] and CO cooling", "GIZMO H2 + HD cooling", "Inverse Compton cooling (CMB)"]
RECOMBINATION = [f"Gas-phase recombination of {i}" for i in ("H+", "He+", "He++")]
FREE_FREE = [f"Free-free emission from {i}" for i in ("H+", "He+", "He++")]


def _state(Tv):
    """A partly ionized, partly molecular state with metals, at Tv"""
    from ..symbols import G_0, NH, grad_v, Z_dust
    vals = {T: Tv, n_Htot: 50.0, X_H: 0.7, z: 0.0, sp.Symbol("y"): 0.0994, G_0: 2.0, NH: 1e21, grad_v: 1e-13,
            Z_dust: 1.0, sp.Symbol("f_metal"): 1.0, sp.Symbol("f_neb"): 1.0, sp.Symbol("Td"): 20.0, T_bg: 12.0,
            sp.Symbol("f_IR_selfabs"): 0.8, sp.Symbol("f_recNUV"): 0.6, sp.Symbol("x_C,tot"): 2.4e-4,
            sp.Symbol("f_d"): 1.0, sp.Symbol("C_2"): 1.0}
    xs = {"H+": 0.3, "He+": 0.02, "He++": 1e-3, "H_2": 0.1, "e-": 0.33}
    xs["H"] = 1 - xs["H+"] - 2 * xs["H_2"]
    for el, x in (("C", 2.4e-4), ("N", 6.8e-5), ("O", 4.9e-4), ("Ne", 8.5e-5), ("Mg", 3.2e-5), ("Si", 3.2e-5),
                  ("S", 1.3e-5), ("Ca", 2.2e-6), ("Fe", 2.5e-5)):
        xs[el] = x
    xs["C+"], xs["CO"] = 0.0, 0.0
    vals.update({x_(s): v for s, v in xs.items()})
    vals.update({n_(s): v * vals[n_Htot] for s, v in xs.items()})
    return vals


def _value(expr, vals):
    return float(expr.xreplace(vals).subs(vals))


@pytest.mark.parametrize("Tv", [80.0, 8.0e3, 3.0e5])
def test_band_outputs_are_gizmos_routing(models, Tv):
    """L_NUV = -fcorr (metal + nebular + H/He+ excitation + f_recNUV recombination + free-free above 1e5 K) and
    L_IR_gas = -fcorr (molecular/fine-structure + Compton + free-free below 1e5 K), from the legacy terms with T_bg in
    the bath factor; photoelectric heating, collisional ionization, cosmic rays and the dust coupling are in neither"""
    rt, legacy = models
    vals = _state(Tv)
    out = {o.name: _value(o.expr, vals) for o in rt.network.outputs if o.name != "photoionization_rate"}
    heat = {p.name: (p.heat if "Compton" in p.name else p.heat.subs(z, T_bg / 2.73 - 1))  # bath factor at T_bg
            for p in legacy.subprocesses}

    def h(names, w=1):
        return sum(_value(heat[n], vals) for n in names) * w

    fcorr, frec, hot = 0.8, 0.6, 1.0 if Tv >= 1e5 else 0.0
    expect_nuv = -fcorr * (h(NUV_TERMS) + frec * h(RECOMBINATION) + hot * h(FREE_FREE))
    expect_ir = -fcorr * (h(IR_TERMS) + (1 - hot) * h(FREE_FREE))
    assert out["L_NUV"] == pytest.approx(expect_nuv, rel=1e-6)
    assert out["L_IR_gas"] == pytest.approx(expect_ir, rel=1e-6)
    dust_row = process(rt, "Gas-dust collisions").network["dust heat"].rhs  # the dust_heat output, not self-absorbed
    assert _value(dust_row, vals) == pytest.approx(-h(["Gas-dust collisions"]), rel=1e-9)


def test_ir_self_absorption_scales_all_heat_but_dust_and_pdv(models):
    """fcorr multiplies Heat and Lambda (photoheating included), not the gas-dust coupling nor the hydro work; the
    chemistry is unscaled"""
    rt, _ = models
    f = sp.Symbol("f_IR_selfabs")
    for p in rt.subprocesses:
        has_f = f in sp.sympify(p.heat).free_symbols
        assert has_f == (p.name not in ("Gas-dust collisions", "PdV work") and p.heat != 0), p.name
        for k, e in p.network.items():
            if k != "heat":
                assert f not in sp.sympify(e.rhs).free_symbols, (p.name, k)


def test_photoionization_rate_output(models):
    """photoionization_rate = Gamma_eff n_H0 at nHcgs: eps_HI times it is the photoheating"""
    rt, _ = models
    vals = {**_state(8e3), Gamma_HI: 2e-9, sigma_HI: 3e-18, eps_HI: 4.8e-12, c_tilde: 3e6, dt: 1e11}
    out = {o.name: o.expr for o in rt.network.outputs}
    heat = process(rt, "Photoionization of H by the RT band").heat
    assert 4.8e-12 * _value(out["photoionization_rate"], vals) == pytest.approx(_value(heat, vals) / 0.8, rel=1e-12)
