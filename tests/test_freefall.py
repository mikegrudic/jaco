"""The explicit right-hand side and the free-fall collapse on the toy model."""
import numpy as np

from jaco.freefall import free_fall, t_freefall
from jaco.models.toy import make_model
from jaco.ode import OdeSystem


def test_rates_are_the_rows():
    s = OdeSystem(make_model())
    assert s.species == ["H+"] and s.balances == []
    r = s.rates(1e4, {"H+": 0.5}, {"n_Htot": 1.0, "pdv_work": 0.0})
    assert np.isfinite([r["heat"], r["dTdt"], r["u"], r["P"], r["rho"], r["cv"]]).all()
    # PdV work enters the heat row one to one
    r2 = s.rates(1e4, {"H+": 0.5}, {"n_Htot": 1.0, "pdv_work": 1e-20})
    assert abs((r2["heat"] - r["heat"]) / 1e-20 - 1) < 1e-10
    # dT/dt is the energy equation through the EOS at fixed composition
    assert abs(r2["dTdt"] - r["dTdt"] - 1e-20 / r["rho"] / r["cv"]) < 1e-12 * abs(r["dTdt"])


def test_free_fall_runs():
    tr = free_fall(make_model(), n0=0.1, T0=300.0, n_stop=1e8, n_points=50)
    n, T, x = tr["n_Htot"], tr["T"], tr["H+"]
    assert np.all(np.diff(n) > 0) and n[0] == 0.1 and abs(n[-1] / 1e8 - 1) < 1e-9
    assert np.all(np.isfinite(T)) and np.all(T > 0) and np.all((x >= 0) & (x <= 1))
    assert tr["pdv_work"][-1] > 0 and abs(tr["pdv_work"][-1] - tr["P"][-1] / t_freefall(tr["rho"][-1])) < 1e-6 * tr["pdv_work"][-1]
    assert np.all(np.diff(tr["t"]) > 0)
