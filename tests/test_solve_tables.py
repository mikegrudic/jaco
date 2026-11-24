"""Model.solve with jaco's interpolants: the 1D piecewise tables and the 2D/3D uniform tables"""
import warnings

import jax
import jax.numpy as jnp
import numpy as np
import pytest
import sympy as sp
from astropy import units
from scipy.optimize import brentq

from jaco.declarations import Species
from jaco.interpolation import tables_in
from jaco.model import Model
from jaco.model_check import _c_lambdify
from jaco.numerics.interp import jax_lambdify
from jaco.processes import Reaction, ThermalTerm
from jaco.symbols import T, n_, n_Htot, piecewise_linear, piecewise_powerlaw, table_interp_2d, table_interp_3d


@pytest.fixture(scope="module", autouse=True)
def x64():
    with jax.enable_x64(True):
        yield


def quiet(f, *a, **k):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return f(*a, **k)


RNG = np.random.default_rng(1)
X, Y, Z = sp.symbols("x y z")
T2 = table_interp_2d("t2", RNG.random((5, 7)), [np.logspace(0, 4, 5), np.linspace(-1, 2, 7)], X, Y,
                     log_axes=[True, False])
T3 = table_interp_3d("t3", RNG.random((4, 5, 6)), [np.linspace(0, 3, 4), np.logspace(-2, 2, 5), np.linspace(0, 1, 6)],
                     Z, X, Y, log_axes=[False, True, False])
P1 = piecewise_powerlaw([1, 10, 1e3, 1e4], [2.0, 5.0, 1.0, 7.0], X, extrapolate=True, name="pl")
P2 = piecewise_linear([0.0, 1.0, 3.0], [1.0, -1.0, 4.0], Y, name="lin")
# points inside and outside every table
POINTS = np.array([RNG.uniform(0.2, 3e4, 300), RNG.uniform(-2, 3, 300), RNG.uniform(-1, 4, 300)])


def test_jax_interpolants_match_the_generated_code_semantics():
    """values and the symbolic derivative nodes, against jaco.check's emulation of the C helpers"""
    exprs = [T2, T3, P1, P2, T2.diff(X), T2.diff(Y), T3.diff(Z), T3.diff(X), T3.diff(Y), P1.diff(X), P2.diff(Y),
             sp.exp(-T2) * P2 + T3**2]
    got = np.array([np.asarray(v) for v in jax_lambdify([X, Y, Z], exprs, cse=True)(*map(jnp.asarray, POINTS))])
    want = np.array(_c_lambdify([X, Y, Z], exprs, tables_in(exprs))(*POINTS))
    np.testing.assert_allclose(got, want, rtol=1e-12, atol=0)


def test_autodiff_of_a_table_is_its_derivative_node():
    """jacfwd through a lookup gives what jaco_table.h returns, zero along a clamped axis included"""
    f = jax_lambdify([X, Y, Z], [T2, T3, P1], cse=True)
    J = np.asarray(jax.vmap(jax.jacfwd(lambda v: jnp.array(f(v[0], v[1], v[2]))))(jnp.asarray(POINTS.T)))
    nodes = [T2.diff(X), T2.diff(Y), T3.diff(X), T3.diff(Y), T3.diff(Z), P1.diff(X)]
    want = np.array(_c_lambdify([X, Y, Z], nodes, tables_in(nodes))(*POINTS))
    got = np.array([J[:, 0, 0], J[:, 0, 1], J[:, 1, 0], J[:, 1, 1], J[:, 1, 2], J[:, 2, 0]])
    np.testing.assert_allclose(got, want, rtol=1e-10, atol=1e-14)
    assert (want[0] == 0).any() and (want[4] == 0).any()  # the points do reach the clamped region


def test_thermal_equilibrium_with_a_1d_table():
    """the toy model's cooling curve is a piecewise power law"""
    from jaco.models.toy import make_model, cooling
    n = np.array([1e-2, 1.0, 1e2])
    sol = quiet(make_model().solve, {"n_Htot": n, "pdv_work": np.zeros(3)},
                {"u": np.full(3, 1e12), "T": np.full(3, 5e3), "H+": np.full(3, 1e-3)}, time_dependent=[], tol=1e-10)
    Tv, xp = np.asarray(sol["T"]), np.asarray(sol["H+"])
    assert np.all((Tv > 10) & (Tv < 2e4)) and np.all(np.diff(Tv) < 0)
    # the energy equation at the returned state, from the processes themselves
    net = quiet(lambda: make_model().network)
    heat = net["heat"].rhs
    syms = sorted(heat.free_symbols, key=str)
    vals = {"T": Tv, "n_Htot": n, "pdv_work": 0.0, "n_H": n * (1 - xp), "n_H+": n * xp, "n_e-": n * xp, "C_2": 1.0}
    net_heat, gross = _c_lambdify(syms, [heat, cooling.heat], tables_in([heat]))(
        *[np.broadcast_to(vals[str(s)], n.shape) for s in syms])
    assert np.all(np.abs(net_heat) < 1e-6 * np.abs(gross))


def test_backward_euler_step_with_a_1d_table():
    """one implicit step of the toy model: the returned state satisfies the energy equation it was solved from"""
    from jaco.models.toy import make_model
    model = quiet(make_model)
    net, u0, dt, n = model.network, 1.2e12, 1e12, 1.0
    sol = quiet(model.solve, {"n_Htot": [n], "pdv_work": [0.0], "u_initial": [u0]},
                {"u": [u0], "T": [8e3], "H+": [1e-2]}, time_dependent=["T"], dt=dt * units.s, tol=1e-12)
    Tv, xp, u = (float(np.asarray(sol[k])[0]) for k in ("T", "H+", "u"))
    state = {"T": Tv, "n_Htot": n, "pdv_work": 0.0, "n_H": n * (1 - xp), "n_H+": n * xp, "n_e-": n * xp, "C_2": 1.0}
    exprs = [net["heat"].rhs, net.eos.density]
    syms = sorted(set().union(*[e.free_symbols for e in exprs]), key=str)
    heat, rho = (float(v) for v in _c_lambdify(syms, exprs, tables_in(exprs))(*[state[str(s)] for s in syms]))
    assert u < u0 and rho * (u - u0) / dt == pytest.approx(heat, rel=1e-6, abs=0)


def test_thermal_equilibrium_with_a_2d_table():
    """a cooling function tabulated in (n, T): the batched solve lands on the root a bracketing solver finds"""
    n_ax, T_ax = np.logspace(-2, 4, 7), np.logspace(1, 5, 9)
    log_lambda = -27 + 1.5 * np.log10(T_ax)[None, :] - 0.3 * np.log10(n_ax)[:, None] ** 2 / 4
    lam = table_interp_2d("lambda_nT", log_lambda, [n_ax, T_ax], n_Htot, T, log_axes=[True, True])
    model = quiet(Model, [ThermalTerm(-(10**lam) * n_Htot**2, name="table cooling"),
                          ThermalTerm(1e-24 * n_Htot**2, name="heating")],
                  solve_vars=["u", "T"], species=[Species("H")])
    n = np.array([0.03, 1.0, 40.0, 3e3])
    sol = quiet(model.solve, {"n_Htot": n}, {"u": np.full(4, 1e12), "T": np.full(4, 300.0)}, tol=1e-12)
    table = _c_lambdify([n_Htot, T], [lam], tables_in([lam]))
    want = [brentq(lambda Tv: 1e-24 - 10 ** float(table(nv, Tv)[0]), 10.0, 1e5, xtol=1e-12, rtol=1e-14)
            for nv in n]
    np.testing.assert_allclose(np.asarray(sol["T"]), want, rtol=1e-8)


def test_fixed_abundance_from_a_table_is_returned():
    """an eliminated species whose abundance is an interpolant is evaluated when the solution is packaged"""
    f_H2 = piecewise_linear([1.0, 2.0, 4.0], [0.4, 0.1, 0.0], sp.log(T, 10), name="f_H2")
    model = quiet(Model, [Reaction("H -> H+ + e-", 1e-16), Reaction("H+ + e- -> H", 2.59e-13 * (T / 1e4) ** -0.7,
                                                                      clumping=1)],
                  solve_vars=["H+"], fixed={"H_2": f_H2}, species=[Species(s) for s in ("H", "H_2", "H+", "e-")])
    sol = quiet(model.solve, {"T": np.array([30.0, 1e3]), "n_Htot": np.array([1.0, 1.0])},
                {"H+": np.array([1e-3, 1e-3])}, time_dependent=[], tol=1e-10)
    np.testing.assert_allclose(np.asarray(sol["H_2"]), [0.4 - 0.3 * (np.log10(30.0) - 1), 0.05], rtol=1e-12)
    np.testing.assert_allclose(np.asarray(sol["H"]) + np.asarray(sol["H+"]) + 2 * np.asarray(sol["H_2"]), 1.0,
                               rtol=1e-12)
