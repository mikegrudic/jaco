"""
Tests the NL99 model chemistry and thermal terms against the equations of Nelson & Langer (1999) and (1997).
"""

import functools
import warnings

import numpy as np
import pytest
import sympy as sp
import jax
import jaco

from astropy import units
from jaco.declarations import Species
from jaco.model import Model
from jaco.processes import ThermalTerm
from jaco.symbols import T, n_

from .. import make_model, model as nelson_model, PAPER_PARAMETERS
from ..chemistry import reactions
from .. import thermal
from ..shielding import co_shielding, LOG_N_CO, LOG_N_H2, THETA


def approx(expected, rel=1e-6):
    """pytest.approx with a purely relative tolerance. Its default absolute tolerance of 1e-12 would make every
    comparison of the small numbers here (rates of 1e-17 to 1e-30 erg cm^-3 s^-1 and s^-1) pass whatever the value"""
    return pytest.approx(expected, rel=rel, abs=0.0)


P = PAPER_PARAMETERS


@pytest.fixture(scope="module", autouse=True)
def x64():
    """
    Running with float32 would break 'test_conservation' below due to rounding errors.
    This sets jax to use float64 for the tests in this file only.
    """
    enable_x64 = getattr(jax, "enable_x64", None)  # top level since jax 0.8
    if enable_x64 is None:
        from jax.experimental import enable_x64
    with enable_x64(True):
        yield


@pytest.fixture(scope="module")
def model():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return make_model()


# ------------------------------------------------------------------ declarations and chemistry

def test_check(model):
    assert jaco.check(model).ok


def test_species_of_the_paper(model):
    """NL99 §2.3 evolves He+, H3+, OHx, CHx, CO, C, C+, HCO+, O, M, M+ and e-; the cloud's H_2 is the background"""
    material = {s.name for s in model.species}
    assert material == {"H_2", "H3+", "He", "He+", "C", "C+", "CH", "O", "OH", "CO", "HCO+", "Mg", "Mg+", "e-"}
    assert set(model.solve_vars) == {"u", "T", "H3+", "He+", "CH", "OH", "CO", "C+", "HCO+", "Mg+"}


def table_2(Tv):
    """
    Rate coefficients of NL99 Table 2, transcribed here:
        key = reactants (the photon and cosmic ray left out)
    at G_0 = 1, A_V = 1, unshielded
    """
    table = {
        ("H_2",)       : 1.2e-17,
        ("He",)        : 6.8e-18,
        ("C", "H3+")   : 2e-9,
        ("H3+", "O")   : 8e-10,
        ("CO", "H3+")  : 1.7e-9,
        ("H_2", "He+") : 7e-15,
        ("CO", "He+")  : 1.6e-9,
        ("C+", "H_2")  : 4e-16,
        ("C+", "OH")   : 1e-9,
        ("CH", "O")    : 2e-10,
        ("C", "OH")    : 5.8e-12 * Tv**0.5,
        ("e-", "He+")  : 9e-11 * Tv**-0.64,
        ("e-", "H3+")  : 1.9e-6 * Tv**-0.54,
        ("C+", "e-")   : 1.4e-10 * Tv**-0.61,
        ("HCO+", "e-") : 3.3e-5 * Tv**-1.0,
        ("Mg+", "e-")  : 3.8e-10 * Tv**-0.65,
        ("H3+", "Mg")  : 2e-9,
        ("C",)         : 3e-10 * np.exp(-3.0),
        ("CH",)        : 1e-9 * np.exp(-1.5),
        ("CO",)        : 1e-10 * np.exp(-3.0),
        ("OH",)        : 5e-10 * np.exp(-1.7),
        ("Mg",)        : 2e-10 * np.exp(-1.9),
        ("HCO+",)      : 1.5e-10 * np.exp(-2.5),
    }
    return {tuple(sorted(k)): v for k, v in table.items()}


@pytest.mark.parametrize("Tv", [10.0, 60.0, 300.0])
def test_rate_coefficients_are_table_2(Tv):
    expected = table_2(Tv)
    seen = set()
    for r in reactions():
        key = tuple(sorted(r.reactants))
        assert key in expected, r.equation
        k = float(sp.sympify(r.rate_coefficient).subs({T: Tv, "G_0": 1, "A_V": 1, "N_CO": 1, "N_H2": 1}))
        assert k == approx(expected[key], rel=1e-12), r.equation
        seen.add(key)
    assert seen == set(expected)


# ------------------------------------------------------------------ thermal terms

def evaluate(term, **values):
    names = sorted(term.heat.free_symbols, key=str)
    return float(term.heat.subs({s: values[str(s)] for s in names}))


def densities(n, **extra):
    """n_Htot and the species densities of the molecular cloud with n = n(H_2) = n_Htot / 2"""
    return {"n_Htot": 2.0 * n, "n_H_2": n, **extra}


def test_cosmic_ray_heating():
    """NL97 eq. 11: 2.0e-17 s^-1 x 20 eV = 6.4e-28 erg per H_2"""
    assert evaluate(thermal.cosmic_ray_heating(), **densities(1e3)) == approx(6.4e-28 * 1e3, rel=3e-3)


def test_h2_formation_heating():
    """NL97 eq. 12"""
    cr = evaluate(thermal.cosmic_ray_heating(), **densities(1e3))
    assert evaluate(thermal.h2_formation_heating(), **densities(1e3)) == approx(cr / 3 * 2.0 / 4.48, rel=1e-12)


@pytest.mark.parametrize("n,Tv", [(1e2, 15.0), (1e4, 30.0), (1e6, 80.0)])
def test_molecular_line_cooling(n, Tv):
    """NL97 eq. 7: alpha (T/10)^beta n^delta, beta = 0.6 + 0.41 log n, delta = 2.3 - 0.18 log n"""
    beta, delta = 0.6 + 0.41 * np.log10(n), 2.3 - 0.18 * np.log10(n)
    expected = -5.0e-30 * (Tv / 10) ** beta * n**delta
    assert evaluate(thermal.molecular_line_cooling(), T=Tv, **densities(n)) == approx(expected, rel=1e-12)


def test_cplus_cooling():
    """NL97 eq. 8, n_H (1 - f/2) = n(H I) + n(H_2) = n for the molecular cloud, n_crit = 3e3"""
    n, Tv, xC = 1e3, 40.0, 5e-5
    boltz = np.exp(-92 / Tv)
    expected = -2.2e-23 * n * (xC * 2 * n) * boltz / (1 + n / 3e3 * (1 + 2 * boltz))
    assert evaluate(thermal.cplus_cooling(), T=Tv, **{"n_C+": xC * 2 * n}, **densities(n)) == approx(expected, rel=1e-12)


def test_oxygen_cooling():
    """NL97 eq. 9"""
    n, Tv, xO = 1e3, 40.0, 1e-4
    expected = -1e-26 * n * (xO * 2 * n) * (24 * np.exp(-228 / Tv) + 7 * np.exp(-326 / Tv)) * Tv**0.5
    heat = evaluate(thermal.oxygen_cooling(), T=Tv, n_O=xO * 2 * n, **densities(n))
    assert heat == approx(expected, rel=1e-12)


# The default dv_dr parameter: 1 km/s/pc in s^-1
DV_DR = 1e5 / 3.0856776e18

def c_i_reference(Tv, n_H2, n_C, dv=DV_DR):
    """A reference implementation of NL99 eq. 7 with Genzel (1991) Table 11. The population fractions f of excited states
    and escape probabilities beta are computed using 500 iterations of the governing eqns to ensure convergence, compared
    to the default 2 iterations used by the jaco implementation.
    """
    k_B, hc_k = 1.380649e-16, 1.438777
    lines = [(609.1354e-4, 7.9e-8, 1e4, 0, 1), (370.415e-4, 2.7e-7, 1.2e3, 1, 2)]  # (lambda [cm], A, n_crit(H_2), l, u)
    g = [1, 3, 5]
    E = [0.0, hc_k / lines[0][0], hc_k / lines[0][0] + hc_k / lines[1][0]]
    beta = [1.0, 1.0]

    for _ in range(500):
        R = np.zeros((3, 3))  # R[i, j]: rate i -> j
        for (lam, A, nc, l, u), b in zip(lines, beta):
            C_ul = n_H2 * A / nc
            R[u, l] += C_ul + A * b
            R[l, u] += C_ul * g[u] / g[l] * np.exp(-(E[u] - E[l]) / Tv)

        M = R.T - np.diag(R.sum(axis=1))
        M[-1, :] = 1
        f = np.linalg.solve(M, [0, 0, 1])

        new = []
        for (lam, A, nc, l, u) in lines:
            tau = A * lam**3 / (8 * np.pi * dv) * n_C * (f[l] * g[u] / g[l] - f[u]) + 1e-8
            new.append((1 - np.exp(-tau)) / tau)

        if np.allclose(new, beta, rtol=1e-13, atol=0):
            break
        beta = new

    return n_C * k_B * sum((E[u] - E[l]) * A * b * f[u] for (lam, A, nc, l, u), b in zip(lines, beta))


@functools.lru_cache(maxsize=None)
def _lambdified(passes):
    heat = thermal.carbon_cooling(passes=passes).heat
    names = sorted(heat.free_symbols, key=str)
    return [str(x) for x in names], sp.lambdify(names, heat, "numpy")


def c_i(Tv, n_H2, n_C, passes=thermal.C_I_PASSES):
    """
    Cooling of the C I term per cm^3.

    The lambdified expression is cached (as it is large).
    """
    names, f = _lambdified(passes)
    values = {"T": Tv, "n_H_2": n_H2, "n_C": n_C, "dv_dr": DV_DR}
    return -float(f(*[values[n] for n in names]))


def test_carbon_cooling_scales_as_collisions_at_low_density():
    """In the optically thin regime and below the critical density, every collisional excitation is radiated. So the
    carbon cooling term is approximately linear in both n(C) and n(H_2), which the test checks.

    The J = 2 level is reached in two collisional steps, so its (370 um) share grows as n(H_2)^2: ~0.3% per doubling at
    50 K, which sets the tolerance of the n(H_2) check"""
    base = c_i(50.0, 1.0, 1e-9)
    assert base > 0
    assert c_i(50.0, 1.0, 2e-9) == approx(2 * base, rel=1e-6)
    assert c_i(50.0, 2.0, 1e-9) == approx(2 * base, rel=1e-2)


def test_carbon_cooling_saturates_at_lte():
    """
    At high density, the levels are Boltzmann populated and the cooling no longer grows with n(H_2), which the test checks
    """
    assert c_i(50.0, 1e10, 1e-3) == approx(c_i(50.0, 1e8, 1e-3), rel=2e-2)


def test_carbon_cooling_matches_genzel_two_level_formula():
    """
    Genzel (1991) eq. 24 for the 609 um line alone, n_u / n_tot = (g_u/g_l) exp(-x) / (1 + (g_u/g_l) exp(-x) + n_crit/n),
    with n_crit = 1e4 for H_2 (Table 11); at 4 K the J = 2 level and the 370 um line are negligible (at 8 K the 370 um line
    adds 5%)
    """
    Tv, n = 4.0, 3e3
    boltz = 3 * np.exp(-1.438777 / 609.1354e-4 / Tv)
    expected = 1.380649e-16 * (1.438777 / 609.1354e-4) * 7.9e-8 * boltz / (1 + boltz + 1e4 / n)
    assert c_i(Tv, n, 1e-9) / 1e-9 == approx(expected, rel=5e-3)


@pytest.mark.parametrize("Tv", [8.0, 15.0, 30.0, 100.0])
@pytest.mark.parametrize("n_H2", [1e2, 1e3, 1e4, 1e5, 1e6])
@pytest.mark.parametrize("frac", [1e-8, 1e-6, 1e-5, 3e-5])
def test_carbon_cooling_agrees_with_the_converged_solution(Tv, n_H2, frac):
    """
    Here we compare the jaco carbon cooling term to a reference implementation that converges to a solution.

    The jaco implementation uses two steps of an iterative method to compute the correct term. For one iteration it
    is within 20% of the correct answer and 4% for two iterations (the worst cases measured are 16% and 3%). Two steps
    was chosen over three to avoid the computational cost.
    """
    tolerance = {1: 0.2, 2: 0.04}[thermal.C_I_PASSES]
    assert c_i(Tv, n_H2, frac * n_H2) == approx(c_i_reference(Tv, n_H2, frac * n_H2), rel=tolerance)


@pytest.mark.parametrize("Tv, expected", [(10.0, 2.0e-28), (20.0, 1.1e-27), (40.0, 2.5e-27)])
def test_carbon_cooling_against_the_C_I_peak_of_goldsmith_and_langer(Tv, expected):
    """
    GL78 Fig. 1: C I cooling per H_2 at the peak of its curve (n(H_2) = 4e3 cm^-3, X(C I) = 1e-5 per H_2), read from
    the plot.
    """
    n = 4e3
    assert c_i(Tv, n, 1e-5 * n) / n == approx(expected, rel=0.25)


def test_surface_cooling():
    """NL97 eq. 10"""
    n, Tv, tau = 10.0, 50.0, 0.7
    expected = -7.31e-25 * n * np.exp(-200 / Tv) * np.exp(-2 * tau)
    assert evaluate(thermal.surface_cooling(), T=Tv, tau_uv=tau, **densities(n)) == approx(expected, rel=1e-12)


def test_dust_temperature():
    """HTT91 eq. 5, T_0 = 12.2 G_0^0.2, floored at 10 K as in NL99"""
    def htt91(G0, AV):
        T0 = 12.2 * G0**0.2
        return (8.9e-11 * 3e15 * G0 * np.exp(-1.8 * AV) + 2.75**5
                + 3.4e-2 * (0.42 - np.log(3.45e-2 * 1e-3 * T0)) * 1e-3 * T0**6) ** 0.2

    for G0, AV in [(1.0, 0.0), (10.0, 0.3), (100.0, 1.0), (1e3, 0.0), (1.0, 3.0)]:
        got = float(thermal.dust_temperature().subs({"G_0": G0, "A_V": AV}))
        assert got == approx(max(10.0, htt91(G0, AV)), rel=1e-12)

    # Below we test the floor in the dust_temperature term.
    assert float(thermal.dust_temperature().subs({"G_0": 1.0, "A_V": 3.0})) == 10.0


def test_gas_dust_exchange():
    """NL97 eq. 13: the gas cools to cooler dust and is unchanged at T_d (10 K in shielded gas)"""
    term, n = thermal.gas_dust_exchange(), 1e5
    expected = -2.4e-33 * 30.0**0.5 * (30.0 - 10.0) * n**2
    assert evaluate(term, T=30.0, G_0=1.0, A_V=3.0, **densities(n)) == approx(expected, rel=1e-12)
    assert evaluate(term, T=10.0, G_0=1.0, A_V=3.0, **densities(n)) == 0
    assert evaluate(term, T=10.0, G_0=1e3, A_V=0.0, **densities(n)) > 0  # hot dust in a strong field warms the gas


@pytest.mark.parametrize("Tv,xe,G0,tau", [(60.0, 5e-5, 1.0, 0.0), (100.0, 1e-4, 10.0, 1.0), (30.0, 1e-6, 1.0, 0.5)])
def test_grain_charging_factor(Tv, xe, G0, tau):
    """HTT91 B7-B9: f = (1 - x)^2 / x, x = 1 - 0.28 / (gamma + 0.5), gamma = 2.9e-4 Y G_0 sqrt(T) exp(-tau) / n(e)"""
    n_e = xe * 1e3
    gamma = 2.9e-5 * G0 * np.sqrt(Tv) * np.exp(-tau) / n_e
    x = 1 - 0.28 / (gamma + 0.5)
    got = float(thermal.grain_charging_factor().subs({T: Tv, "n_e-": n_e, "G_0": G0, "tau_uv": tau}))
    assert got == approx((1 - x) ** 2 / x, rel=1e-12)


def test_grain_charging_limits():
    """Uncharged grains (weak field or many electrons) heat at (1 - 0.44)^2 / 0.44 = 0.71 of the unit-efficiency rate,
    and strongly charged grains (strong field, few electrons) hardly heat the gas"""
    f = lambda n_e, G0: float(thermal.grain_charging_factor().subs({T: 100.0, "n_e-": n_e, "G_0": G0, "tau_uv": 0.0}))
    assert f(1e3, 1e-3) == approx((1 - 0.44) ** 2 / 0.44, rel=1e-4)  # gamma -> 0: x = 1 - 0.28 / 0.5
    assert f(1e-6, 1e3) < 1e-3 * f(1e3, 1e-3)


def test_photoelectric_heating():
    """NL99 eq. 8 with HTT91's charging, per H nucleus: 4.86e-26 G_0 n_H f exp(-tau_uv)"""
    n_H, Tv, G0, tau, n_e = 2e3, 50.0, 2.0, 0.5, 0.05
    f = float(thermal.grain_charging_factor().subs({T: Tv, "n_e-": n_e, "G_0": G0, "tau_uv": tau}))
    expected = 4.86e-26 * G0 * n_H * f * np.exp(-tau)
    heat = evaluate(thermal.photoelectric_heating(), T=Tv, n_Htot=n_H, **{"n_e-": n_e}, G_0=G0, tau_uv=tau)
    assert heat == approx(expected, rel=1e-12)


# ------------------------------------------------------------------ CO shielding, van Dishoeck & Black (1988) Table 5

shield = sp.lambdify(sp.symbols("N_CO N_H2"), co_shielding(*sp.symbols("N_CO N_H2")), modules="numpy")


def test_shielding_reproduces_the_table():
    for i, ly in enumerate(LOG_N_H2):
        for j, lx in enumerate(LOG_N_CO):
            assert shield(10.0**lx, 10.0**ly) == approx(THETA[i, j], rel=1e-9)


@pytest.mark.parametrize("N_CO,N_H2,expected", [
    (1e15, 1e21, 0.1360),  # printed entry
    (1e17, 1e20, 8.769e-3),
    (1e19, 1e23, 3.572e-7),
])
def test_shielding_printed_values(N_CO, N_H2, expected):
    assert shield(N_CO, N_H2) == approx(expected, rel=1e-3)


def test_shielding_interpolation_and_limits():
    """Bilinear in ln Theta between nodes, unshielded for no columns, and held at the table edge beyond it"""
    mid = shield(10**13.5, 1.0)
    assert mid == approx(np.sqrt(THETA[0, 1] * THETA[0, 2]), rel=1e-9)
    assert shield(1.0, 1.0) == 1.0 and shield(0.0, 0.0) == 1.0
    assert shield(1e25, 1e25) == approx(THETA[-1, -1], rel=1e-9)
    N = np.logspace(0, 19, 60)
    assert (np.diff(shield(N, 1e21)) <= 1e-15).all()  # more CO shields more
    assert (np.diff(shield(1e16, np.logspace(0, 23, 60))) <= 1e-15).all()  # more H2 shields more


def test_co_photodissociation_carries_the_shielding():
    r = next(r for r in reactions() if r.equation == "CO -> C + O")
    subs = {T: 60.0, "G_0": 1.0, "A_V": 1.0, "N_CO": 1e16, "N_H2": 1e21}
    k = float(sp.sympify(r.rate_coefficient).subs(subs))
    assert k == approx(1e-10 * np.exp(-3.0) * THETA[3, 4], rel=1e-9)


# ------------------------------------------------------------------ initial abundances

def abundances(He_plus=1e-9, H3_plus=1e-9):
    """The initial conditions of the NL99 appendix per H nucleus: carbon as C+, oxygen atomic, M+, He neutral"""
    xC, xO, y, xM = P["x_C,tot"], P["x_O,tot"], P["y"], P["x_Mg,tot"]
    x = {"H_2": 0.5, "H3+": H3_plus, "He": y - He_plus, "He+": He_plus, "C": 0.0, "C+": xC, "CH": 0.0, "O": xO,
         "OH": 0.0, "CO": 0.0, "HCO+": 0.0, "Mg": 0.0, "Mg+": xM}
    x["e-"] = xC + He_plus + H3_plus + xM
    return x


# ------------------------------------------------------------------ the model through Model.solve

@pytest.fixture(scope="module")
def chemistry_model():
    """The reactions alone: T is an input, so the solve is the chemistry at fixed temperature"""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return Model(reactions(), solve_vars=list(nelson_model.SOLVED), time_dependent=list(nelson_model.SOLVED),
                     fixed={"H_2": nelson_model.H2_PER_H}, parameters=nelson_model.PARAMETERS,
                     species=[Species(s) for s in nelson_model.MATERIAL])


def chemistry_step(model, dt, temperature, density, extinction):
    """Take one batched backward-Euler step from the NL99 initial abundances."""
    initial = abundances()
    inputs = {"T": temperature, "n_Htot": density, "A_V": extinction, "G_0": 1.0, **P}
    inputs.update({f"x_{s}_initial": initial[s] for s in model.time_dependent})
    knowns = dict(zip(inputs, np.broadcast_arrays(*inputs.values())))
    shape = knowns["T"].shape
    guess = {s: np.full(shape, max(initial[s], 1e-15)) for s in model.time_dependent}
    return model.solve(knowns, guess, dt=dt * units.s, tol=1e-10)


def test_conservation(chemistry_model):
    """He, C, O, Mg and charge are conserved by the reductions, for the solved and the eliminated species"""
    x = chemistry_step(chemistry_model, 1e11, temperature=60.0, density=(1e2, 1e3, 1e4), extinction=0.0)
    assert np.all(x["CO"] > abundances()["CO"])
    np.testing.assert_allclose(x["He"] + x["He+"], P["y"], rtol=1e-9, atol=0)
    np.testing.assert_allclose(x["C"] + x["C+"] + x["CO"] + x["CH"] + x["HCO+"], P["x_C,tot"], rtol=1e-9, atol=0)
    np.testing.assert_allclose(x["O"] + x["OH"] + x["CO"] + x["HCO+"], P["x_O,tot"], rtol=1e-9, atol=0)
    np.testing.assert_allclose(x["Mg"] + x["Mg+"], P["x_Mg,tot"], rtol=1e-9, atol=0)
    ions = x["C+"] + x["He+"] + x["H3+"] + x["Mg+"] + x["HCO+"]
    np.testing.assert_allclose(x["e-"], ions, rtol=1e-9, atol=0)
    assert np.all(x["H_2"] == 0.5)
    assert all(np.all(np.asarray(v) >= 0) for v in x.values())


def test_co_forms_with_extinction(chemistry_model):
    """The photodissociation of CO and of the CH precursor weakens with A_V, so CO grows with the extinction"""
    x = chemistry_step(chemistry_model, 1e12, temperature=30.0, density=1e3, extinction=(0.0, 1.0, 3.0))
    assert np.all(np.diff(x["CO"]) > 0)


@pytest.fixture(scope="module")
def thermal_model():
    """The thermal terms alone: u and T unknown, the abundances inputs"""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return Model(processes=[*thermal.thermal_terms(), ThermalTerm(sp.Symbol("pdv_work"), name="PdV work")],
                     solve_vars=["u", "T"], time_dependent=["T"], parameters=nelson_model.PARAMETERS,
                     species=[Species(s) for s in nelson_model.MATERIAL])


def test_thermal_balance_has_a_root(thermal_model):
    """Heating wins at low T and cooling at high T, so a cloud relaxes to a temperature between, whether it starts
    at 3 K or at 300 K. A step of 1e19 s makes the implicit energy equation the equilibrium, whatever u_initial is."""
    x = abundances()
    knowns = {"n_Htot": 1e3, "pdv_work": 0.0, "G_0": 1.0, "A_V": 1.0, "tau_uv": 2.3, "dv_dr": 3.2408e-14, **P}
    knowns |= {f"x_{s}": x[s] for s in ("H_2", "H3+", "He+", "CH", "OH", "CO", "C+", "HCO+", "Mg+")}
    knowns = {k: np.array([v]) for k, v in knowns.items()}
    knowns["u_initial"] = np.array([1e-12])
    T_eq = [float(np.asarray(thermal_model.solve(knowns, {"u": np.array([1e-12]), "T": np.array([T0])},
                                                 dt=1e19 * units.s, tol=1e-10)["T"])[0]) for T0 in (3.0, 300.0)]
    assert 3.0 < T_eq[0] < 300.0
    assert T_eq[1] == approx(T_eq[0], rel=1e-4)


def test_coupled_timestep_balances_chemistry_and_energy(model):
    """The full Nelson model satisfies backward Euler and the EOS at the same final state."""
    network = model.network
    eos = network.eos
    initial = abundances()
    T0, n_H, dt = 30.0, 1e3, 1e11  # K, cm^-3, s
    initial_state = {T: T0, **{n_(s): n_H * x for s, x in initial.items()}}
    u0 = float(eos.internal_energy.xreplace(initial_state))
    params = {"n_Htot": n_H, "G_0": 1.0, "A_V": 1.0, "tau_uv": 2.3, "N_CO": 1e15, "N_H2": 1e21,
              "dv_dr": 3.2408e-14, "C_2": 1.0, "pdv_work": 0.0, **P}
    knowns = {**params, "u_initial": u0, **{f"x_{s}_initial": initial[s] for s in nelson_model.SOLVED}}
    guess = {"u": u0, "T": T0, **{s: max(initial[s], 1e-15) for s in nelson_model.SOLVED}}
    solution = model.solve({k: [v] for k, v in knowns.items()}, {k: [v] for k, v in guess.items()},
                           dt=dt * units.s, tol=1e-10)
    final = {k: float(np.asarray(v).item()) for k, v in solution.items()}
    assert 0.0 < final["T"] < T0
    assert final["CO"] > initial["CO"]

    # Evaluate the original rate equations at the returned composition, including the eliminated species.
    final_state = {sp.Symbol(k): v for k, v in params.items()}
    final_state.update({T: final["T"], **{n_(s): n_H * final[s] for s in initial}})
    for s in nelson_model.SOLVED:
        rate = float(network[s].rhs.xreplace(final_state))
        step_rate = n_H * (final[s] - initial[s]) / dt
        assert step_rate == approx(rate), s

    heat = float(network["heat"].rhs.xreplace(final_state))
    rho = float(eos.density.xreplace(final_state))
    assert rho * (final["u"] - u0) / dt == approx(heat)
    assert final["u"] == approx(float(eos.internal_energy.xreplace(final_state)))


def test_two_body_rates_carry_no_clumping(model):
    """The reactions carry jaco's default clumping factor, which the model fixes at 1: the solved system has no C_2
    for a host to leave unset"""
    C_2 = sp.Symbol("C_2")
    assert C_2 in model.network.symbols
    assert C_2 not in model.network.reduced({}, list(model.time_dependent)).symbols
