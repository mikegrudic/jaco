"""GIZMO's CMB-bath and high-temperature truncation factors for the low-temperature cooling block."""

import numpy as np
import pytest
import sympy as sp
from ..symbols import T, z, cmb_bath_factor, lowtemp_truncation, dust_sputtering_truncation
from ..H2_cooling import H2_cooling, H2_cooling_rate
from ..gas_dust_collisions import gas_dust_heattransfer_coeff
from ..symbols import sqrt_T, Z_dust, f_dust


def gizmo_truncation_factor(Tv):
    logT = np.log10(Tv)
    if logT > 5.3:  # GIZMO skips the whole block; the jaco factor is < 1.2e-7 there
        return 0.0
    dx = (logT - 4.5) / 0.20
    return np.exp(-min(dx * dx, 40.0)) if logT > 4.5 else 1.0


def gizmo_dust_factor(Tv):
    factor = gizmo_truncation_factor(Tv)
    if Tv > 3.0e5:
        dx = (Tv - 3.0e5) / 2.0e5
        factor *= np.exp(-min(dx * dx, 40.0))
    return factor


trunc = sp.lambdify(T, lowtemp_truncation, modules="numpy")
dust = sp.lambdify(T, lowtemp_truncation * dust_sputtering_truncation, modules="numpy")
cmb = sp.lambdify((T, z), cmb_bath_factor, modules="numpy")


@pytest.mark.parametrize("Tv", [10.0, 1e3, 1e4, 10**4.5, 5e4, 1e5, 10**5.29, 1e6, 1e8])
def test_truncation_factors(Tv):
    assert trunc(Tv) == pytest.approx(gizmo_truncation_factor(Tv), rel=1e-12, abs=2e-7)
    assert dust(Tv) == pytest.approx(gizmo_dust_factor(Tv), rel=1e-12, abs=2e-7)


@pytest.mark.parametrize("Tv,zv", [(2.0, 0.0), (10.0, 0.0), (100.0, 0.0), (30.0, 6.0)])
def test_cmb_bath_factor(Tv, zv):
    T_cmb = 2.73 * (1 + zv)
    assert cmb(Tv, zv) == pytest.approx((Tv - T_cmb) / (Tv + T_cmb), rel=1e-12)


def test_factors_applied():
    assert sp.simplify(H2_cooling.heat + H2_cooling_rate() * cmb_bath_factor * lowtemp_truncation) == 0
    bare = 1.116e-32 * sqrt_T * (1.0 - 0.8 * sp.exp(-75.0 / T)) * Z_dust * f_dust
    assert sp.simplify(gas_dust_heattransfer_coeff() - bare * lowtemp_truncation * dust_sputtering_truncation) == 0
