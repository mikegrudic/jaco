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
