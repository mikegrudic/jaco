"""STARFORGE's CO cooling against Whitworth & Jaffa (2018, A&A 611, A20; WJ18) evaluated in their own variables.

WJ18 fit Goldsmith & Langer (1978) for gas of X = 0.70, Y = 0.28 with all of its H in H2 (their Sec. 1): mean particle
mass 3.97e-24 g, mass per H2 molecule rho / n_H2 = 4.77e-24 g, and X_CO = n_CO / n_H2. Their Eq. 37 gives the cooling
per unit volume in rho, the isothermal sound speed a = (k T / m)^(1/2), X_CO and |div v|. Its two limits agree with
Eqs. 15-16 to the rounding of the printed coefficients (2%), but its beta prefactor (0.813) is not Eq. 17's under the
same normalisation (1.25), and Fig. 1 plots Eqs. 15-18: beta is taken from Eq. 17.

The gas here has WJ18's X, so rho = m_p n_H / X, and n_CO is the carbon outside C+ times the molecular fraction 2 x_H2.
The model's velocity-gradient norm stands in for |div v|.
"""

import numpy as np
import pytest
import sympy as sp
from jaco.symbols import n_, x_, boltzmann_cgs, protonmass_cgs
from ..starforge import make_model
from ..gas_phase import x_gas
from ..symbols import T, grad_v, n_Htot, z, G_0

X_WJ18 = 0.70  # H mass fraction of WJ18's gas
M_PARTICLE, M_H2 = 3.97e-24, 4.77e-24  # g: WJ18 Sec. 1
LO_COEFF, HI_COEFF = 4.63e8, 4.58e-45  # Eq. 37 (cgs)
X_C_TOT = 2.7e-4  # carbon per H nucleus, gas and dust
C2 = sp.Symbol("C_2")


def wj18_cooling(T_gas, n_H, n_CO, div_v):
    """Eq. 37 (LO and HI limits) with Eq. 17's beta, erg cm^-3 s^-1"""
    rho = protonmass_cgs * n_H / X_WJ18
    n_H2 = rho / M_H2
    a = np.sqrt(boltzmann_cgs * T_gas / M_PARTICLE)
    lo = LO_COEFF * (n_CO / n_H2) * rho**2 * a**3
    hi = HI_COEFF * a**8 * div_v
    beta = 1.23 * n_H2**0.0533 * T_gas**0.164
    return (lo ** (-1 / beta) + hi ** (-1 / beta)) ** -beta


@pytest.fixture(scope="module")
def model():
    return make_model()


def state(T_gas, n_H, x_H2, div_v):
    """Numerical values of every symbol the CO cooling and the model's carbon partition read"""
    return {T: T_gas, n_Htot: n_H, x_("H_2"): x_H2, n_("H_2"): x_H2 * n_H, sp.Symbol("x_C,tot"): X_C_TOT, G_0: 1.0,
            grad_v: div_v, C2: 1.0, z: 0.0}


def jaco_cooling(model, vals):
    """-heat of the model's CO cooling per unit volume, with the model's CO abundance"""
    heat = model.processes["CO Cooling"].heat.subs({n_("CO"): x_("CO") * n_Htot})
    return -float(heat.subs(x_("CO"), model.fixed["CO"]).subs(vals))


def model_n_CO(model, vals):
    """n_CO the WJ18 evaluation is given: the gas-phase carbon outside C+ times the molecular fraction 2 x_H2"""
    x_Cplus = float(model.fixed["C+"].subs(vals))
    return (float(x_gas("C").subs(vals)) - x_Cplus) * 2 * vals[x_("H_2")] * vals[n_Htot]


# (T, n_H, |div v| in s^-1): optically thin (LO), LVG (HI) and in between
REGIMES = {"LO": (30.0, 1e2, 1e-10), "HI": (10.0, 1e7, 1e-16), "mid": (20.0, 1e4, 1e-13)}


@pytest.mark.parametrize("x_H2", [0.5, 0.25, 0.05], ids=["molecular", "half", "atomic"])
@pytest.mark.parametrize("regime", list(REGIMES))
def test_CO_cooling_matches_WJ18(model, regime, x_H2):
    T_gas, n_H, div_v = REGIMES[regime]
    vals = state(T_gas, n_H, x_H2, div_v)
    cmb_bath = (T_gas - 2.73) / (T_gas + 2.73)  # GIZMO's factor on its molecular cooling, which the model keeps
    expected = wj18_cooling(T_gas, n_H, model_n_CO(model, vals), div_v) * cmb_bath
    assert jaco_cooling(model, vals) == pytest.approx(expected, rel=0.03, abs=0)


def test_regimes_are_the_limits_they_name(model):
    """The LO and HI states sit deep in their limits, so each tests its coefficient alone"""
    for regime in ("LO", "HI"):
        T_gas, n_H, div_v = REGIMES[regime]
        vals = state(T_gas, n_H, 0.5, div_v)
        n_CO = model_n_CO(model, vals)
        rho = protonmass_cgs * n_H / X_WJ18
        a = np.sqrt(boltzmann_cgs * T_gas / M_PARTICLE)
        limit = LO_COEFF * n_CO * M_H2 * rho * a**3 if regime == "LO" else HI_COEFF * a**8 * div_v
        assert wj18_cooling(T_gas, n_H, n_CO, div_v) / limit == pytest.approx(1, abs=0.01)
