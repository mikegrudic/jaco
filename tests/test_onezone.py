"""The explicit right-hand side and the one-zone modes on the toy model."""
import numpy as np

from jaco.models.toy import make_model
from jaco.ode import OdeSystem
from jaco.onezone import OneZone, t_freefall


def test_rates_are_the_rows():
    s = OdeSystem(make_model())
    assert s.species == ["H+"] and s.balances == []
    r = s.rates(1e4, {"H+": 0.5}, {"n_Htot": 1.0, "pdv_work": 0.0})
    assert np.isfinite([r["heat"], r["dTdt"], r["u"], r["P"], r["rho"], r["cv"]]).all()
    r2 = s.rates(1e4, {"H+": 0.5}, {"n_Htot": 1.0, "pdv_work": 1e-20})
    assert abs((r2["heat"] - r["heat"]) / 1e-20 - 1) < 1e-10  # PdV work enters the heat row one to one
    assert abs(r2["dTdt"] - r["dTdt"] - 1e-20 / r["rho"] / r["cv"]) < 1e-12 * abs(r["dTdt"])  # and dT/dt through the EOS


def test_free_fall_runs():
    tr = OneZone(make_model()).free_fall(n0=0.1, T0=300.0, n_stop=1e8, n_points=50)
    n, T, x = tr["n_Htot"], tr["T"], tr["H+"]
    assert np.all(np.diff(n) > 0) and n[0] == 0.1 and abs(n[-1] / 1e8 - 1) < 1e-9
    assert np.all(np.isfinite(T)) and np.all(T > 0) and np.all((x >= 0) & (x <= 1))
    assert abs(tr["pdv_work"][-1] - tr["P"][-1] / t_freefall(tr["rho"][-1])) < 1e-6 * tr["pdv_work"][-1]
    assert np.all(np.diff(tr["t"]) > 0)


def test_equilibrium_and_evolve_agree():
    zone = OneZone(make_model())
    eq = zone.equilibrium([10.0, 100.0], T0=1e4)
    assert eq["relaxed"].all() and np.all(eq["T"] > 0)
    # the time evolution at fixed density relaxes to the same state
    ev = zone.evolve(100.0, T0=1e4, t_end=eq["t"][-1] * 3, n_points=20)
    assert abs(ev["T"][-1] / eq["T"][-1] - 1) < 1e-3
    assert abs(ev["H+"][-1] - eq["H+"][-1]) < 1e-3 * max(eq["H+"][-1], 1e-6)
    assert np.allclose(ev["n_Htot"], 100.0, rtol=1e-12) and np.all(ev["pdv_work"] == 0)


def test_steady_species_are_solved_and_u_round_trips():
    s = OdeSystem(make_model())  # the toy's H+ is a steady-state variable: solved, not integrated
    assert s.integrated == [] and s.steady == ["H+"]
    p = {"n_Htot": 1.0, "pdv_work": 0.0}
    for T in (300.0, 1.5e4, 5e4):
        s.reset()
        r = s.rates(T, {"H+": 0.3}, p)  # the 0.3 is only a starting guess
        row = s._evaluate(T, r["x"], {}, {**s.defaults, **p})[s._i_species[0]]
        assert abs(row) < 1e-8 * abs(s._system(T, None, r["x"], {}, {**s.defaults, **p})[1][0, 0]) + 1e-300
        T2, x2 = s.temperature(r["u"], {}, p, T_guess=1e3)
        assert abs(T2 / T - 1) < 1e-8 and abs(x2["H+"] / r["x"]["H+"] - 1) < 1e-6
    assert s.rhs(s.state(1e4, {}, p), p).shape == (1,)


def test_gizmo_host_closures_are_finite():
    from jaco.onezone import GizmoHostClosures as H
    st = {"rho": 2e-22, "lambda_J": 1e19, "dlnrho_dt": 1e-14, "t_ff": 1e14, "T": 50.0, "x": {}, "params": {"Z_d": 1.0, "ISRF": 1.0}}
    vals = [H.N_H(st), H.G_0(st), H.A_dust(st), H.kappa_thick(st), H.grad_v(st)]
    assert np.all(np.isfinite(vals)) and np.all(np.array(vals) > 0)
    assert H.kappa_thick(st) >= 0.1 and H.kappa_thick({**st, "params": {"Z_d": 0.0}}) == 0.1
    assert H.A_dust({**st, "params": {"Z_d": 0.0}}) == 0.0
    assert H.G_0(st) < 1.7 * st["params"]["ISRF"]


def test_plain_host_closures_are_callable():
    from jaco.onezone import Closures
    st = {"rho": 2e-22, "P": 1e-12, "lambda_J": 1e19, "dlnrho_dt": 1e-14, "t_ff": 1e14, "n_Htot": 100.0}
    vals = {k: f(st) for k, f in Closures.DEFAULT.items()}
    assert vals["∇v"] == 1e-14 and vals["N_H"] == 100.0 * 1e19 / 2 and np.all(np.isfinite(list(vals.values())))


def test_free_fall_steps_converge_to_free_fall():
    """The implicit stepping (free_fall_steps: GIZMO's one-zone driver loop, one backward-Euler step of the model per
    Delta_t, OdeSystem(step=True).solve_step) against the integrated collapse (free_fall). The density follows free
    fall exactly, so the elapsed time is 2 t_ff(n0) (1 - sqrt(n0 / n)) at every step; the temperature approaches the
    integrated track at first order as the step shrinks. Compared from 1.5 n0, past the toy's initial heating from 300 K
    to ~1900 K at nearly fixed density (where T(n) is close to vertical), and up to 1e6, where the toy, which has no
    temperature floor, is still above the step solver's 1 K"""
    n0, T0, n_stop = 0.1, 300.0, 1e6
    ref = OneZone(make_model(), rtol=1e-9, atol=1e-14).free_fall(n0=n0, T0=T0, n_stop=n_stop, n_points=2000)
    zone = OneZone(make_model(), step=True)
    err = []
    for f_ff in (0.04, 0.01):
        tr = zone.free_fall_steps(n0=n0, T0=T0, n_stop=n_stop, f_ff=f_ff)
        n, T, x = tr["n_Htot"], tr["T"], tr["H+"]
        assert n[-1] >= n_stop and np.all(np.diff(n) > 0) and np.all(np.isfinite(T)) and np.all((x >= 0) & (x <= 1))
        t_ff0 = t_freefall(tr["rho"][0])  # the gizmo host's free-fall density is the gas density
        assert np.allclose(tr["t"], 2 * t_ff0 * (1 - np.sqrt(n0 / n)), rtol=1e-9, atol=0)
        m = (n >= 1.5 * n0) & (n <= n_stop)
        lnT_ref = np.interp(np.log(n[m]), np.log(ref["n_Htot"]), np.log(ref["T"]))
        err.append(np.max(np.abs(np.log(T[m]) - lnT_ref)))
    assert err[1] < 1e-4  # in ln T, at Delta_t <= 0.01 t_ff (5e-6 as written)
    assert err[1] < err[0] / 2.5  # first order: a quarter of the step, about a quarter of the error
