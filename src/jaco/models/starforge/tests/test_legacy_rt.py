"""starforge_legacy_RT and starforge_legacy_RT_EUV against direct transcriptions of GIZMO's legacy RT coupling
(cooling/cooling.cc, radiation/rt_utilities.cc, rt_dust_opacity.cc, eos/eos.cc at gizmo_jaco_dev 2df0c6dd)."""

import math

import numpy as np
import pytest
import sympy as sp
from scipy.optimize import brentq

from jaco.symbols import n_, x_, dt, sanitize_symbols
from jaco.model_check import _c_lambdify
from jaco.interpolation import tables_in
from jaco.processes.recombination import gasphase_recombination_rates
from jaco.processes.ionization import collisional_ionization_rates
from ..symbols import T, n_Htot, X_H, z
from .. import radiation as rt
from ..radiation import EUV, FUV, NUV, ONIR, IR, C_LIGHT, EV
from .. import dust_opacity as do
from ...starforge_legacy import make_model as make_legacy
from ...starforge_legacy_RT import make_model, NUV_ROUTE, IR_ROUTE, DERIVED
from ...starforge_legacy_RT_EUV import make_model as make_euv

S = sp.Symbol
MP, KB, X_GIZMO = 1.6726e-24, 1.38066e-16, 0.76


@pytest.fixture(scope="module")
def models():
    return make_model(), make_legacy(), make_euv()


def process(model, name):
    return next(p for p in model.subprocesses if p.name == name)


def state(Tv=8e3, **over):
    """Symbol -> value of a partly ionized, partly molecular, irradiated state with metals"""
    v = {T: Tv, n_Htot: 50.0, X_H: 0.7, z: 0.0, S("y"): 0.0994, S("Z_d"): 1.0, S("N_H"): 1e21, S("∇v"): 1e-13,
         S("f_metal"): 1.0, S("f_neb"): 1.0, S("x_C,tot"): 2.4e-4, S("Td"): 25.0, S("T_rad"): 30.0, S("T_CMB"): 2.73,
         S("rho"): 50 * 2.34e-24, S("Z_metals"): 0.014, S("gamma_eos"): 5.0 / 3.0, S("u_initial"): 3e11,
         S("rsol"): 1e-3, S("sigma_HI"): 3e-18, S("eps_HI"): 4.8e-12, S("hnu_EUV"): 21.0, S("Δx"): 3e18,
         S("G_LW_bg"): 0.0, S("gamma_12_UVB"): 0.0, S("eps_H0_UVB"): 0.0, dt: 1e11, S("C_2"): 1.0, S("C_3"): 1.0}
    xs = {"H+": 0.3, "He+": 0.02, "He++": 1e-3, "H_2": 0.1, "e-": 0.33,
          EUV: 1e-3, FUV: 0.02, NUV: 0.01, ONIR: 0.5, IR: 2.0}
    for el, x in (("C", 2.4e-4), ("N", 6.8e-5), ("O", 4.9e-4), ("Ne", 8.5e-5), ("Mg", 3.2e-5), ("Si", 3.2e-5),
                  ("S", 1.3e-5), ("Ca", 2.2e-6), ("Fe", 2.5e-5)):
        xs[el] = x
    xs["C+"], xs["CO"] = 0.0, 0.0
    for k, val in over.items():
        if k in xs:
            xs[k] = val
        else:
            v[T if k == "T" else n_Htot if k == "n_Htot" else S(k)] = val
    if "Td_initial" not in over:
        v[S("Td_initial")] = v[S("Td")]  # the solve starting from the state's dust temperature
    xs["H"] = 1 - xs["H+"] - 2 * xs["H_2"]
    v.update({x_(s): val for s, val in xs.items()})
    v.update({n_(s): val * v[n_Htot] for s, val in xs.items()})
    v[S("x_photon_IR_initial")] = xs[IR]
    return v


KICK = dict(rt.kick_absorption_intermediates())


def value(expr, vals, derived=None):
    """expr at vals, the model's derived expressions and the kick-factor intermediates substituted first, with the
    generated code's float semantics"""
    expr = sp.sympify(expr)
    if derived:
        expr = expr.xreplace({S(k): e for k, e in derived.items()})
    expr = expr.xreplace(KICK)
    syms = sorted(expr.free_symbols, key=str)
    f = _c_lambdify([sanitize_symbols(s) for s in syms], [sanitize_symbols(expr)], tables_in([expr]))
    return float(np.ravel(f(*[np.array([float(vals[s])]) for s in syms])[0])[0])


# --- GIZMO transcriptions ----------------------------------------------------------------------------------------

def gizmo_dust_survival(Td):
    s = 9 * (1 - Td / 1500.0)
    return max(0.5 * (1 + s / math.sqrt(1 + s * s)) * math.exp(-min(40.0, (Td / 1500.0) ** 2 / 9)), 1e-25)


def gizmo_planck_mean(Trad, Tdust):
    """dust_planck_mean_opacity: zones by T_dust, linear in log T_rad, constant beyond the table"""
    zone = next((i for i, b in enumerate((160, 275, 425, 680, 1500)) if Tdust < b), 4)
    logT = math.log10(Trad)
    if logT >= 4:
        return 10 ** do.LOG_KAPPA[zone][-1]
    if logT <= 0:
        return 10 ** do.LOG_KAPPA[zone][0]
    idx = int(14 * logT / 4)
    w1 = 1 - (logT - do.LOG_T_RAD[idx]) / (do.LOG_T_RAD[1] - do.LOG_T_RAD[0])
    return 10 ** (w1 * do.LOG_KAPPA[zone][idx] + (1 - w1) * do.LOG_KAPPA[zone][idx + 1])


def gizmo_kappa_ir(Td, Trad, flag_ea, flag_dg, c):
    """rt_kappa_adaptive_IR_band [cm^2/g]; c: Ne, HII, HI, fmol, zmetals, Zfac, rho, Tgas"""
    if flag_ea == 1:
        Trad = Td
    dtm = gizmo_dust_survival(Td)
    kappa = 0.0
    if flag_dg >= 0:
        k = gizmo_planck_mean(Trad, Td)
        if flag_ea in (1, -1):
            k *= 1 - 0.5 / (1 + 725.0**2 / (1 + Trad**2))
        kappa += k * c["Zfac"] * dtm
    if flag_dg <= 0:
        X, xe, rho, Tg = X_GIZMO, c["Ne"], c["rho"], c["Tgas"]
        f_neutral, f_free = max(0, 1 - xe), c["zmetals"] * max(0, 1 - 0.5 * dtm)
        k_e = 0.4 * X * xe / ((1 + 2.7e11 * rho / Tg**2) * (1 + (Trad / 4.5e8) ** 0.86))
        k_mol = 0.1 * (f_free + 3e-9) * f_neutral * c["fmol"]
        k_K = 4.0e25 * (1 + X) * (f_free * math.exp(-min(1.5e5 / Trad, 40)) + 0.001 * xe) * rho / (Trad**3 * math.sqrt(Tg))
        k_K += (1.5e20 * f_free * rho / Trad**2 * math.exp(-min((0.8e4 / Trad) ** 4, 40))
                * math.exp(-min((Trad / 0.7e6) ** 2, 40)))
        k_R = f_neutral * min(5e-19 * Trad**4, 0.2 * (1 + X))
        tg = (Tg / 1.3e4) ** 2
        x_Hm = 4e-10 * Tg * xe * c["HI"] / ((1 + c["HII"] * 300 + xe * 1000 * tg / (1 + tg) + 4e-17) * (1 + Tg / 3e4))
        k_bf = 4.2e7 * (8760 / Trad) ** 1.5 * math.exp(-min(8760 / Trad, 40))
        phi = min(Tg / 5040, 2)
        k_ff = (1.9e6 * (8760 / Trad) ** 2 * math.exp(-min(8760 / Trad, 40))
                * (0.6 - 2.5 * math.sqrt(phi) + 2.5 * phi + 2.7 * phi * math.sqrt(phi)))
        k_rad = k_mol + k_K + x_Hm * (k_bf + k_ff) + k_e + k_R
        if flag_ea in (1, -1):
            k_rad -= k_e
        kappa += k_rad
    return kappa


def cell_of(vals):
    """the cell quantities rt_kappa_adaptive_IR_band reads, at a state"""
    tg = 1 + 0.59 * (vals[S("gamma_eos")] - 1) * (MP / KB) * vals[S("u_initial")]
    return dict(Ne=vals[x_("e-")], HII=vals[x_("H+")], HI=1 - vals[x_("H+")], fmol=2 * vals[x_("H_2")],
                zmetals=vals[S("Z_metals")], Zfac=vals[S("Z_d")], rho=vals[S("rho")], Tgas=tg)


def gizmo_kappa_band(kappa0, floor, vals):
    """rt_kappa for the photoelectric, NUV and optical bands [cm^2/g]"""
    zf = vals[S("Z_d")] * gizmo_dust_survival(vals[S("Td")])
    return max(0.02 + 0.35 * vals[x_("e-")], kappa0 * (max(floor, zf) if floor else zf))


def gizmo_gas_dust_coeff(Tg, Td, Zsol):
    return 1.116e-32 * math.sqrt(Tg) * (1 - 0.8 * math.exp(-75.0 / Tg)) * Zsol * gizmo_dust_survival(Td)


def nHcgs(vals):
    return X_GIZMO / vals[X_H] * vals[n_Htot]


def u_band(vals, band):
    """band energy density [erg cm^-3]"""
    return vals[n_(band)] * EV * (vals[S("hnu_EUV")] if band == EUV else 1)


# --- photoionization ----------------------------------------------------------------------------------------------

def test_photoionization_is_gizmos(models):
    """Ionizations per volume c sigma n_gamma nHcgs x_H, GIZMO's Gamma HI nHcgs on the atomic H (find_abundances_and_rates;
    GIZMO's HI also counts the nuclei in H2, see rt.photoionization); the band loses one photon each at c_tilde (the
    kick: kappa rho = sigma HI nHcgs, rt_kappa); the optical band gains hnu_EUV each (the donation); the gas gains eps_HI
    each (Heat_Ion_from_RHD), times fcorr"""
    rt_model, _, euv = models
    v = state()
    ionizations = C_LIGHT * v[S("sigma_HI")] * v[n_(EUV)] * v[x_("H")] * nHcgs(v)
    pi = [process(rt_model, name) for name in rt.EUV_SINKS]
    rows = {k: sum(value(p.network[k].rhs, v) for p in pi if k in p.network) for k in (EUV, ONIR, "H+")}
    assert -rows[EUV] == pytest.approx(v[S("rsol")] * ionizations, rel=1e-12)
    assert rows[ONIR] == pytest.approx(v[S("rsol")] * v[S("hnu_EUV")] * ionizations, rel=1e-12)
    fcorr = value(S("f_IR_selfabs"), v, DERIVED)
    assert sum(value(p.heat, v, DERIVED) for p in pi) == pytest.approx(fcorr * 4.8e-12 * ionizations, rel=1e-12)
    # the H+ row per n_Htot is GIZMO's dx_HII/dt, abundances per nHcgs (rule GIZMO_ABUNDANCES)
    assert rows["H+"] / v[n_Htot] == pytest.approx(ionizations / nHcgs(v), rel=1e-12)
    pi_euv = [process(euv, name) for name in rt.EUV_SINKS]
    assert -sum(value(p.network[EUV].rhs, v) for p in pi_euv) == pytest.approx(v[S("rsol")] * ionizations, rel=1e-12)
    assert all(ONIR not in p.network for p in pi_euv)  # no optical band: GIZMO's kick loses the absorbed energy


@pytest.mark.parametrize("xHp", [1e-4, 0.3, 0.99])
def test_H_balance_runs_at_gizmo_speed(models, xHp):
    """The H+ row per n_Htot is GIZMO's d x_H+/dt (cooling.cc 900-904) without H2: rates at nHcgs, abundances per nHcgs"""
    rt_model, _, _ = models
    v = state(**{"H+": xHp, "H_2": 1e-20, "e-": xHp + 2e-4})
    rows = sum(value(process(rt_model, name).network["H+"].rhs, v) for name in
               ("Photoionization of H by the ionizing band", "Gas-phase recombination of H+", "Collisional Ionization of H"))
    Tv, nH, xe = v[T], nHcgs(v), v[x_("e-")]
    aHp = float(gasphase_recombination_rates["H+"].subs(T, Tv))
    geH0 = float(collisional_ionization_rates["H"].subs(T, Tv))
    Gamma = C_LIGHT * v[S("sigma_HI")] * v[n_(EUV)]
    gizmo = Gamma * (1 - xHp) + xe * nH * geH0 * (1 - xHp) - xe * nH * aHp * xHp
    assert rows / v[n_Htot] == pytest.approx(gizmo, rel=1e-9)


# --- dust ------------------------------------------------------------------------------------------------------

@pytest.mark.parametrize("band", [FUV, NUV, ONIR])
@pytest.mark.parametrize("over", [{}, {"Td": 1450.0}, {"e-": 0.9, "H+": 0.9, "Z_d": 1e-3}])
def test_dust_band_absorption_is_the_kicks(models, band, over):
    """The kick takes e0 (1 - exp(-a dt)), a = c_tilde f_abs kappa rho (rt_update_driftkick, f_abs = 1/2 from
    rt_absorb_frac_albedo, kappa from rt_kappa), from the band, and the dust gains it at the true c (dust_dE_cooling's
    c/c_tilde): the backward-Euler step with the band's row lands on the kick's band, and the dust's row is the
    band's loss at c. Where rt_kappa's neutral/electron floor exceeds the dust's opacity (the last case, for the NUV
    and optical bands), the factor's exponent takes the dust's opacity alone and the band ends near the kick's"""
    rt_model, _, _ = models
    p = process(rt_model, f"Dust absorption of {band}")
    kappa0, floor = rt.DUST_BAND_OPACITY[band]
    for rho_ in (2.3e-22, 2.3e-20, 2.3e-18, 2.3e-14):  # a dt from ~1e-4 to beyond the exponent's cap
        v = state(**over)
        v[S("rho")] = rho_
        d, rsol_ = v[dt], v[S("rsol")]
        kappa = gizmo_kappa_band(kappa0, floor, v)
        zf = v[S("Z_d")] * gizmo_dust_survival(v[S("Td")])
        kappa_dust = kappa0 * (max(floor, zf) if floor else zf)
        a = rsol_ * C_LIGHT * 0.5 * kappa * rho_
        e1 = u_band(v, band)
        x = min(a * d * kappa_dust / kappa, rt.KICK_EXPONENT_CAP)
        e0 = e1 * (1 + a * d * math.expm1(x) / (a * d * kappa_dust / kappa))  # the implicit step back from e1
        if kappa_dust >= kappa:
            assert e0 == pytest.approx(e1 * math.exp(min(a * d, rt.KICK_EXPONENT_CAP)), rel=1e-9)  # the kick's exponential
        loss = -value(p.network[band].rhs, v) * EV * d
        assert loss == pytest.approx(e0 - e1, rel=1e-9)
        assert value(p.network["dust heat"].rhs, v) * d == pytest.approx((e0 - e1) / rsol_, rel=1e-9)


@pytest.mark.parametrize("over", [{}, {"e-": 0.9, "H+": 0.9}, {"rho": 2.3e-18}])
def test_ir_band_gets_the_donations_twice(models, over):
    """rt_update_driftkick: each donor band adds its absorbed energy de_abs to the IR band and to E_abs_tot_toIR; the
    IR band, last, reads its energy (now with the donations) as e0 and adds E_abs_tot_toIR dt again in total_de_dt,
    keeping both (its own absorption it re-emits but for the gas share). Transcribed over a step that ends at the
    state's donor energies: the IR band gains twice what the donors lose, as the model's copy plus the dust's
    re-emission of what it absorbs from them"""
    rt_model, _, _ = models
    v = state(**over)
    d = v[dt]
    E_abs_tot_toIR, first_copy = 0.0, 0.0
    for b in (FUV, NUV, ONIR):  # donors first, an exponential absorption each at c_tilde
        a = v[S("rsol")] * C_LIGHT * 0.5 * gizmo_kappa_band(*rt.DUST_BAND_OPACITY[b], v) * v[S("rho")]
        e1 = u_band(v, b)
        de_abs = e1 * math.expm1(min(a * d, rt.KICK_EXPONENT_CAP))  # from the e0 that ends at e1
        E_abs_tot_toIR += de_abs / d
        first_copy += de_abs
    kick_gain = first_copy + E_abs_tot_toIR * d  # no IR absorption: its re-emission returns it but for the gas share
    copy = value(process(rt_model, "GIZMO's second copy of the donated dust absorption in photon_IR").network[IR].rhs, v)
    absorbed = sum(value(process(rt_model, f"Dust absorption of {b}").network["dust heat"].rhs, v) for b in (FUV, NUV, ONIR))
    model_gain = (copy * EV + v[S("rsol")] * absorbed) * d  # the dust balance re-emits what it absorbs
    assert model_gain == pytest.approx(kick_gain, rel=1e-9)


def _near_switch(Td):
    return any(abs(Td - b) < h for b, h in zip(do.ZONE_BOUNDARIES, do.ZONE_HALF_WIDTHS))


IR_GRID = [(Td, Trad, Tg) for Td in (5.0, 40.0, 120.0, 200.0, 350.0, 550.0, 900.0, 1400.0, 3000.0)
           for Trad in (3.0, 30.0, 300.0, 3e3, 3e4) for Tg in (20.0, 3e3, 2e4)]


def test_ir_opacities_are_gizmos():
    """The dust absorption (flags -1, 1), emission (1, 1) and gas absorption (-1, -1) opacities against
    rt_kappa_adaptive_IR_band, exact away from the smoothed composition switches"""
    v0 = state()
    fns = {"abs": do.ir_dust_opacity(do.T_dust, do.T_rad), "em": do.ir_dust_opacity(do.T_dust, do.T_dust),
           "gas": do.ir_gas_opacity(do.T_rad, do.T_dust)}
    for Td, Trad, Tg in IR_GRID:
        u = (Tg - 1) / (0.59 * (2.0 / 3.0) * MP / KB)  # the opacity's gas temperature estimate is then Tg
        v = {**v0, S("Td"): Td, S("T_rad"): Trad, S("u_initial"): u}
        c = cell_of(v)
        ref = {"abs": gizmo_kappa_ir(Td, Trad, -1, 1, c), "em": gizmo_kappa_ir(Td, Trad, 1, 1, c),
               "gas": gizmo_kappa_ir(Td, Trad, -1, -1, c)}
        for k, e in fns.items():
            if k != "gas" and _near_switch(Td):
                continue
            assert value(e, v) == pytest.approx(ref[k], rel=1e-9, abs=1e-300), (k, Td, Trad, Tg)


def test_zone_switches_are_bounded():
    """Within a switch's window the smoothed opacity lies between the two zones' values; outside the windows the dust
    emission per surviving dust, T^4 kappa(T, T) / f_dust, increases with T_dust (sublimation above ~1500 K removes the
    dust from every term of its balance alike)"""
    e = do.T_dust**4 * do.ir_dust_opacity(do.T_dust, do.T_dust).subs(do.Z_dust, 1) / do.dust_survival(do.T_dust)
    em = _c_lambdify([sanitize_symbols(do.T_dust)], [sanitize_symbols(e)], tables_in([e]))
    Td = np.array([t for t in np.linspace(3, 5000, 100001) if not _near_switch(t)])
    rising = np.diff(em(Td)[0]) > 0
    assert np.all(rising | (np.diff(Td) > 0.1))  # only across a window may it fall
    for b, h in zip(do.ZONE_BOUNDARIES, do.ZONE_HALF_WIDTHS):
        for t in (b - 0.5 * h, b, b + 0.5 * h):
            got = value(do.semenov_planck_mean(30.0, t), {})
            lo, hi = sorted((gizmo_planck_mean(30.0, b - 1e-6), gizmo_planck_mean(30.0, b + 1e-6)))
            assert lo * (1 - 1e-9) <= got <= hi * (1 + 1e-9)


def test_dust_emission_and_ir_rows(models):
    """Dust emission 4 sigma kappa_P(Td) rho Td^4 into the IR band at c_tilde (dust_dEdt); dust and gas absorption of the
    IR band at kappa(Td, T_rad), the gas heated at c_tilde/c of the physical rate as GIZMO's kick does"""
    rt_model, _, _ = models
    v = state(Td=60.0, T_rad=45.0)
    c = cell_of(v)
    E_ir = u_band(v, IR)
    emission = 4 * 5.67e-5 * gizmo_kappa_ir(60.0, 60.0, 1, 1, c) * v[S("rho")] * 60.0**4
    dust_abs = C_LIGHT * gizmo_kappa_ir(60.0, 45.0, -1, 1, c) * v[S("rho")] * E_ir
    gas_abs = C_LIGHT * gizmo_kappa_ir(60.0, 45.0, -1, -1, c) * v[S("rho")] * E_ir
    r = v[S("rsol")]
    p = process(rt_model, "Dust emission into photon_IR")
    assert value(p.network["dust heat"].rhs, v) == pytest.approx(-emission, rel=1e-12)
    assert value(p.network[IR].rhs, v) * EV == pytest.approx(r * emission, rel=1e-12)
    p = process(rt_model, "Dust absorption of photon_IR")
    assert value(p.network["dust heat"].rhs, v) == pytest.approx(dust_abs, rel=1e-12)
    assert value(p.network[IR].rhs, v) * EV == pytest.approx(-r * dust_abs, rel=1e-12)
    p = process(rt_model, "Gas absorption of photon_IR")
    d = v[dt]
    x = r * C_LIGHT * gizmo_kappa_ir(v[S("Td_initial")], 45.0, -1, 1, c) * v[S("rho")] * d  # the kick's dust absorption
    share = 2 * -math.expm1(-x / 2) / x
    assert value(p.heat, v) == pytest.approx(r * share * gas_abs, rel=1e-9)
    assert value(p.network["dust heat"].rhs, v) == pytest.approx(gas_abs, rel=1e-12)  # rt_eqm_dust_temp counts it all
    assert value(p.network[IR].rhs, v) * EV == pytest.approx(-r * (share + 1) * gas_abs, rel=1e-9)  # the share, and the
    # dust balance's re-emission of the gas absorption, which the kick does not make


@pytest.mark.parametrize("rho_", [2.3e-22, 2.3e-18, 2.3e-14])
def test_gas_share_is_two_half_kicks(models, rho_):
    """rt_update_driftkick, twice per step: each half-step kick absorbs e0 (1 - exp(-a dt/2)) of the IR band and gives
    the gas its opacity share fgas of that (DtInternalEnergy), re-emitting the rest; the model's gas heating over the
    step is that, to first order in fgas, from optically thin (the rate) to thick (2 fgas e0 per step)"""
    rt_model, _, _ = models
    v = state(Td=60.0, T_rad=45.0, rho=rho_)
    c = cell_of(v)
    d, r = v[dt], v[S("rsol")]
    k_dust = gizmo_kappa_ir(v[S("Td_initial")], 45.0, -1, 1, c)
    k_gas = gizmo_kappa_ir(60.0, 45.0, -1, -1, c)
    e = u_band(v, IR)
    gas = 0.0
    for _ in range(2):
        de_abs = e * -math.expm1(-r * C_LIGHT * (k_dust + k_gas) * rho_ * d / 2)
        gas += de_abs * k_gas / (k_dust + k_gas)
        e -= de_abs * k_gas / (k_dust + k_gas)
    model = value(process(rt_model, "Gas absorption of photon_IR").heat, v) * d
    assert model == pytest.approx(gas, rel=3 * k_gas / k_dust + 1e-6)


def test_gas_dust_collisions_are_the_dust_balances(models):
    """gas_dust_heating_coeff nHcgs^2 (T - Td) with the surviving dust at Td and no high-temperature truncation
    (rt_ir_lambdadust); the dust reservoir gets what the gas loses"""
    rt_model, _, _ = models
    for Tv in (30.0, 8e3, 5e5):
        v = state(Tv, Td=40.0)
        expected = gizmo_gas_dust_coeff(Tv, 40.0, 1.0) * nHcgs(v) ** 2 * (Tv - 40.0)
        p = process(rt_model, "Gas-dust collisions")
        assert value(p.network["dust heat"].rhs, v) == pytest.approx(expected, rel=1e-12)
        assert value(p.heat, v) == pytest.approx(-expected, rel=1e-12)


def _dust_row_function(model, vals):
    """Td -> the dust's energy balance (the sum of every process's dust heat row) at vals"""
    e = sum(p.network["dust heat"].rhs for p in model.subprocesses if "dust heat" in p.network)
    e = sp.sympify(e).xreplace({S(k): ex for k, ex in DERIVED.items()}).xreplace(KICK)
    e = e.xreplace({s: vals[s] for s in e.free_symbols if s in vals and s != do.T_dust})
    f = _c_lambdify([sanitize_symbols(do.T_dust)], [sanitize_symbols(e)], tables_in([e]))
    return lambda t: float(f(np.array([t]))[0][0])


def gizmo_eqm_dust_temp(vals):
    """rt_eqm_dust_temp's root of dust_dEdt (volumetric): Lambda_gd nHcgs^2 (T - Td) + absorption - 4 sigma kappa_P rho
    Td^4, the absorption (the dust bands at rt_kappa and f_abs, the IR band at the gas and dust absorption opacity,
    flags -1, 0) evaluated at the root itself, as the cooling-side balance (dust_dE_cooling) has it. A dust band's
    absorption is the kick's
    E_abs_tot_toIR (c/c_tilde): e0 (1 - exp(-a dt)) / dt, for the e0 that the kick takes to the band's energy here
    (the exponent with the dust's opacity alone where rt_kappa's floor exceeds it, as the model's intermediate has it)"""
    c = cell_of(vals)
    d, rsol_ = vals[dt], vals[S("rsol")]

    def kick_absorbed(v, b):  # rt_kappa from the cell before the update: the start-of-step dust temperature
        kappa0, floor = rt.DUST_BAND_OPACITY[b]
        v0 = {**v, S("Td"): vals[S("Td_initial")]}
        kappa = gizmo_kappa_band(kappa0, floor, v0)
        zf = v[S("Z_d")] * gizmo_dust_survival(v0[S("Td")])
        x = rsol_ * C_LIGHT * 0.5 * kappa0 * (max(floor, zf) if floor else zf) * v[S("rho")] * d
        a_dt = rsol_ * C_LIGHT * 0.5 * kappa * v[S("rho")] * d
        factor = math.expm1(min(x, rt.KICK_EXPONENT_CAP)) / x if x > 0 else 1.0
        return u_band(v, b) * a_dt * factor / d / rsol_

    def dEdt(Td):
        v = {**vals, S("Td"): Td}
        absorbed = sum(kick_absorbed(v, b) for b in (FUV, NUV, ONIR))
        absorbed += C_LIGHT * gizmo_kappa_ir(Td, vals[S("T_rad")], -1, 0, c) * vals[S("rho")] * u_band(vals, IR)
        emission = 4 * 5.67e-5 * gizmo_kappa_ir(Td, Td, 1, 1, c) * vals[S("rho")] * Td**4
        return gizmo_gas_dust_coeff(vals[T], Td, vals[S("Z_d")]) * nHcgs(vals) ** 2 * (vals[T] - Td) + absorbed - emission
    return _walk_root(dEdt, vals[S("Td")])


def _walk_root(f, guess, lo=2.73, hi=1e4):
    """rt_eqm_dust_temp's search: from the guess, walk by factors 0.9/1.1 (growing) in the direction f points to the
    first sign change, then the root there; or the bound it reaches"""
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


DUST_GRID = [dict(T=Tg, n_Htot=n, rho=n * 2.34e-24, Td=guess, **{FUV: f, ONIR: 25 * f, IR: ir})
             for Tg in (10.0, 100.0, 3e3) for n in (1e2, 1e5, 1e8) for f in (0.0, 1e-3, 10.0) for ir in (1e-3, 1.0, 300.0)
             for guess in (15.0, 600.0)]


def test_dust_temperature_root_is_rt_eqm_dust_temps(models):
    """The steady state of the dust heat row (the model's Td equation) against rt_eqm_dust_temp's fixed point over gas
    temperatures, densities, band energies and starting dust temperatures, both found by GIZMO's walk from the start
    (across a composition switch the balance has two roots, and the walk takes the one it reaches first): equal where
    the root is not within a switch's smoothing window, within the window where it is"""
    rt_model, _, _ = models
    for over in DUST_GRID:
        v = state(**over, T_rad=40.0)
        ref = gizmo_eqm_dust_temp(v)
        got = _walk_root(_dust_row_function(rt_model, v), v[S("Td")])
        near = [(b, h) for b, h in zip(do.ZONE_BOUNDARIES, do.ZONE_HALF_WIDTHS) if abs(ref - b) < h or abs(got - b) < h]
        if near:
            assert abs(got - ref) < 2 * near[0][1], (over, got, ref)
        else:
            assert got == pytest.approx(ref, rel=1e-8), (over, got, ref)


# --- cooling radiation into the bands, Compton, the derived inputs --------------------------------------------------

@pytest.mark.parametrize("Tv,T_rad", [(80.0, 30.0), (8.0e3, 30.0), (3.0e5, 30.0), (8.0e3, 2e4)])
def test_cooling_radiation_routing(models, Tv, T_rad):
    """CoolingRate's return (1317-1348): each routed term's radiated energy (minus its fcorr-scaled heat) into the NUV
    band (the IR band where T_rad > 1e4 K) or the IR band, with the recombination share f_recNUV and the free-free split
    at 1e5 K, at c_tilde"""
    rt_model, _, _ = models
    v = state(Tv, T_rad=T_rad)
    r = v[S("rsol")]
    frec = value(S("f_recNUV"), v, DERIVED)
    hot = 1.0 if Tv >= 1e5 else 0.0
    to_ir = 1.0 if T_rad > 1e4 else 0.0
    for name in set(NUV_ROUTE) | set(IR_ROUTE):
        p = process(rt_model, name)
        lost = -value(p.heat, v, DERIVED)
        if "Free-free" in name:
            w_nuv, w_ir = hot, 1 - hot
        else:
            w_nuv, w_ir = (frec if "recombination" in name else NUV_ROUTE.get(name, 0)), IR_ROUTE.get(name, 0)
        got_nuv = value(p.network[NUV].rhs, v, DERIVED) if NUV in p.network else 0.0
        got_ir = value(p.network[IR].rhs, v, DERIVED) if IR in p.network else 0.0
        assert got_nuv * EV == pytest.approx(r * lost * w_nuv * (1 - to_ir), rel=1e-9, abs=1e-300), name
        assert got_ir * EV == pytest.approx(r * lost * (w_ir + w_nuv * to_ir), rel=1e-9, abs=1e-300), name
    # photoelectric heating, collisional ionization and cosmic rays return nothing; nor does H2 photodissociation take
    for name in ("Photoelectric Heating", "Collisional Ionization of H", "GIZMO H2 network"):
        assert not set(process(rt_model, name).network) & set(rt.BANDS), name


def test_compton_off_bands_is_gizmos(models):
    """evaluate_Compton_heating_cooling_rate's band terms (2410-2440), times nHcgs^2: 2.16e-35 nHcgs u[eV cm^-3]
    (T - T_eff), times x_e where T_eff < 3e4 K"""
    rt_model, _, _ = models
    v = state(1e4, T_rad=50.0)
    teff = {EUV: 2340 * v[S("hnu_EUV")], FUV: 24400.0, NUV: 12000.0, ONIR: 2800.0, IR: 50.0}
    lam = sum(2.16e-35 / nHcgs(v) * u_band(v, b) / EV * (v[x_("e-")] if te < 3e4 else 1.0) * (v[T] - te)
              for b, te in teff.items())
    p = process(rt_model, "Inverse Compton cooling (RT bands)")
    fcorr = value(S("f_IR_selfabs"), v, DERIVED)
    assert value(p.heat, v, DERIVED) == pytest.approx(-fcorr * lam * nHcgs(v) ** 2, rel=1e-9)


def gizmo_bb_frac(E1, E2, Teff):
    k = 8.617e-5
    x1, x2 = E1 / (k * Teff), E2 / (k * Teff)

    def f(x):
        return (131.4045728599595 * x**3) / (2560 + x * (960 + x * (232 + 39 * x))) if x < 3.40309 else \
            1 - 0.15398973382026504 * (6 + x * (6 + x * (3 + x))) * math.exp(-min(x, 40))
    df = f(x2) - f(x1)
    if df <= 0 and x1 > 4:
        return 0.15398973382026504 * (6 + x1 * (6 + x1 * (3 + x1))) * math.exp(-min(x1, 120)) if x1 < 120 else 2e-47
    return max(df, 0)


@pytest.mark.parametrize("over", [{}, {"T_rad": 3e4, FUV: 3e4}, {FUV: 0.0, IR: 0.0}, {"T": 50.0, "T_rad": 20.0},
                                  {"rho": 2.3e-20}, {"rho": 2.3e-18}])
def test_derived_inputs_are_gizmos(over):
    """get_FUV_G0 (RT_PHOTOELECTRIC under M1), update_explicit_molecular_fraction's G_LW, the emission corrections'
    T_bg, CoolingRate's fcorr, the recombination share f_recNUV and the dust survival, from the bands and the gas; the
    photoelectric band as the kick's exponential absorption has it over the step, (e0 - e1) / (a dt), for the e0 that
    the kick takes to the band's energy e1 here"""
    v = state(**over)
    a_dt = v[S("rsol")] * C_LIGHT * 0.5 * gizmo_kappa_band(*rt.DUST_BAND_OPACITY[FUV], v) * v[S("rho")] * v[dt]
    u_pe_end, u_ir = u_band(v, FUV), u_band(v, IR)
    u_pe = u_pe_end * math.expm1(min(a_dt, rt.KICK_EXPONENT_CAP)) / a_dt  # over the step, with the exponent's cap
    habing = 1.6e-3 / C_LIGHT
    assert value(DERIVED["G_0"], v) == pytest.approx(max(1e-56, min(u_pe_end / habing, 1e8)), rel=1e-12, abs=1e-56)
    assert value(DERIVED["G_0_step"], v) == pytest.approx(max(1e-56, min(u_pe / habing, 1e8)), rel=1e-12, abs=1e-56)
    g_lw = min(max(u_pe / habing, 1e-10), 1e10) + u_ir * gizmo_bb_frac(11.2, 500.0, v[S("T_rad")]) / habing
    assert value(DERIVED["G_LW"], v) == pytest.approx(g_lw, rel=1e-9)
    e_cmb, e_ir = 0.262, u_ir / EV
    assert value(DERIVED["T_bg"], v) == pytest.approx((e_ir * v[S("T_rad")] + e_cmb * 2.73) / (e_ir + e_cmb), rel=1e-12)
    tau = gizmo_kappa_ir(v[T], v[T], -1, -1, cell_of(v)) * 0.5 * v[S("rho")] * v[S("Δx")]
    assert value(DERIVED["f_IR_selfabs"], v) == pytest.approx(1 / (1 + tau * tau), rel=1e-9)
    nss = 0.0123 * 10 ** (0.173 * (math.log10(v[T]) - 4))
    q = nHcgs(v) / nss
    shield = 0.98 / (1 + q**1.64) ** 2.28 + 0.02 / (1 + q * (1 + 1e-4 * nHcgs(v) ** 4)) ** 0.84
    heat_rhd = 4.8e-12 * C_LIGHT * 3e-18 * v[n_(EUV)]
    assert value(DERIVED["f_recNUV"], v) == pytest.approx((1 - shield) * heat_rhd / (heat_rhd + 1e-30), rel=1e-9)
    assert value(DERIVED["f_d"], v) == pytest.approx(gizmo_dust_survival(v[S("Td")]), rel=1e-12)


def test_ir_radiation_temperature_output(models):
    """T_rad_new = E_final / (survivors/T_rad + rest x (dust share/Td + gas share/T)), survivors = (E_0 + direct
    donation) exp(-a dt), a the absorption rate, within the range of the three temperatures (GIZMO's photon-number
    weighting, rt_update_driftkick and rt_cooling_radiation_to_bands)"""
    rt_model, _, _ = models
    expr, sums = rt.ir_radiation_temperature(["a"], ["d"], ["g"], ["p"])
    v = state(300.0, Td=60.0, T_rad=40.0)
    n0, d = v[S("x_photon_IR_initial")] * v[n_Htot], v[dt]
    for k, D, G, P in ((1e-13, 3e-12, 1e-12, 0), (1e-13, -3e-12, 1e-12, 0), (1e-13, 3e-12, -1e-11, 0), (0, 0, 0, 0),
                       (1e-13, 3e-12, 1e-12, 2e-12), (1e-9, 3e-12, 1e-12, 2e-12), (1e-9, 3e-12, 0, 0)):
        n1 = (n0 + d * (max(D, 0) + max(G, 0) + P)) / (1 + k * d)  # the implicit step, A = -k n1
        vals = {**v, n_(IR): n1, S("A_IR"): -k * n1, S("D_IR"): D, S("G_IR"): G, S("P_IR"): P}
        kept = (n0 + d * P) * math.exp(-k * d)
        Dp, Gp = max(D, 0), max(G, 0)
        weights = kept / 40.0 + ((n1 - kept) * (Dp / 60.0 + Gp / 300.0) / (Dp + Gp) if Dp + Gp else 0)
        assert value(expr, vals) == pytest.approx(min(max(n1 / weights, 40.0), 300.0), rel=1e-12)
    n1 = (n0 + 1e-3 * d) / (1 + 1e-9 * d)
    vals = {**v, n_(IR): n1, S("A_IR"): -1e-9 * n1, S("D_IR"): 1e-3, S("G_IR"): 0, S("P_IR"): 0}
    assert value(expr, vals) == pytest.approx(60.0, rel=1e-3)  # thick: the band at the dust temperature
    out = next(o for o in rt_model.outputs if o.name == "T_rad_new")
    groups = {k: {p for (p, row), _ in terms} for k, terms in out.sums}
    assert groups["A_IR"] == {"Dust absorption of photon_IR", "Gas absorption of photon_IR"}
    assert groups["D_IR"] == {"Dust emission into photon_IR"}
    assert groups["P_IR"] == {"GIZMO's second copy of the donated dust absorption in photon_IR"}
    emitters = {p.name for p in rt_model.subprocesses if IR in p.network} - groups["A_IR"] - groups["D_IR"]
    assert groups["G_IR"] == emitters - groups["P_IR"]  # every other process that adds to the band is the gas's


def test_declarations(models):
    rt_model, legacy, euv = models
    assert rt_model.solve_vars == legacy.solve_vars + (EUV, FUV, NUV, ONIR, IR, "Td")
    assert rt_model.time_dependent == ("T", "H+", "H_2", EUV, FUV, NUV, ONIR, IR)
    assert [v.name for v in rt_model.variables] == ["Td"] and rt_model.variables[0].row == "dust heat"
    assert euv.solve_vars == legacy.solve_vars + (EUV,) and not euv.variables
    params = {p.name for p in rt_model.parameters}
    assert {"rsol", "sigma_HI", "eps_HI", "hnu_EUV", "T_rad", "rho"} <= params
    assert not {"G_0", "G_0_step", "G_LW", "Td", "f_d", "T_bg", "f_IR_selfabs", "f_recNUV"} & params  # expressions now
    pe = process(rt_model, "Photoelectric Heating").heat.free_symbols
    assert S("G_0_step") in pe and S("G_0") not in pe  # the heating over the step; the closures at its end
    assert {"Td", "G_0", "G_LW"} <= {p.name for p in euv.parameters}  # no IR band: GIZMO's non-RT dust temperature


@pytest.mark.slow
@pytest.mark.parametrize("name", ["starforge_legacy_RT", "starforge_legacy_RT_EUV"])
def test_check(name):
    import jaco
    from importlib import import_module
    from jaco.model_check import CheckGrid
    grid = CheckGrid(np.logspace(np.log10(3.0), 7, 6), np.logspace(-2, 8, 4), n_random=3)
    field = dict(rsol=1e-3, sigma_HI=3e-18, eps_HI=4.8e-12, hnu_EUV=20.0, T_rad=30.0, T_CMB=2.73, rho=2.3e-22,
                 Z_metals=0.014, gamma_eos=5 / 3., Z_d=1.0, **{"Δx": 3e18})
    report = jaco.check(import_module(f"jaco.models.{name}").make_model(), grid, params=field, name=name)
    assert report.ok, report.summary()
