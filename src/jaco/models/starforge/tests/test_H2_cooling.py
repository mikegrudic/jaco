"""H2 + HD cooling (STARFORGE) against the Glover & Abel 2008 thin rates with number-density colliders and the HM79 LTE
limit: n/n_crit is the thin rate per molecule (colliders included) over the LTE rate per molecule, and HD is excited by
collisions with H nuclei, as in GIZMO's CoolingRate."""

import numpy as np
import pytest
import sympy as sp
from jaco.symbols import n_
from ..H2_cooling import H2_cooling_rate
from ..symbols import T, n_Htot, x_

C2 = sp.Symbol("C_2")


def ga08_H2_HD(Tv, n, xH, xH2, xHe, xHp, xe):
    """Number-weighted GA08 thin rates, n/n_crit = thin rate per molecule / HM79 LTE rate per molecule, HD/H = 2 D/H x_H2
    excited by H nuclei, all with clumping 1"""
    T3, logT, q = Tv / 1e3, np.log10(Tv), np.log10(Tv) - 3
    poly = lambda c, v: sum(ci * v**i for i, ci in enumerate(c))
    lam = {"H": 10**poly([-103., 97.59, -48.05, 10.8, -0.9032], logT),
           "He": 10**poly([-23.6892, 2.18924, -0.815204, 0.290363, -0.165962, 0.191914], q),
           "H_2": 10**poly([-23.9621, 2.09434, -0.771514, 0.436934, -0.149132, -0.0336383], q),
           "H+": 10**poly([-21.7167, 1.38658, -0.379153, 0.114537, -0.232142, 0.0585389], q),
           "e-": 10**(poly([-22.1903, 1.5729, -0.213351, 0.961498, -0.910232, 0.137497], q) if logT > 2.30103 else poly([-34.2862, -48.5372, -77.1212, -51.3525, -15.1692, -0.981203], q))}
    thin = n * (xH * lam["H"] + xHe * lam["He"] + xH2 * lam["H_2"] + xHp * lam["H+"] + xe * lam["e-"])
    LTE = (6.7e-19*np.exp(-5.86/T3) + 1.6e-18*np.exp(-11.7/T3) + 3.e-24*np.exp(-0.51/T3) + 9.5e-22*T3**3.76*np.exp(-0.0022/T3**3)/(1+0.12*T3**2.1))
    HD = ((1.555e-25 + 1.272e-26*Tv**0.77)*np.exp(-128./Tv) + (2.406e-25 + 1.232e-26*Tv**0.92)*np.exp(-255./Tv)) * np.exp(-T3*T3/25.)
    xHD = 2 * 2.527e-5 * xH2
    r = thin / LTE
    return n * xH2 * thin / (1 + r) + n * xHD * n * HD / (1 + xHD / xH2 * r)


@pytest.mark.parametrize("Tv,n", [(100.0, 30.0), (4500.0, 105.0), (500.0, 1e5), (2000.0, 1e8)])
def test_starforge_H2_cooling(Tv, n):
    xH, xH2, xHe, xHp, xe = 0.3, 0.35, 0.094, 1e-4, 2e-4
    vals = {T: Tv, n_("H"): xH * n, n_("H_2"): xH2 * n, n_("He"): xHe * n, n_("H+"): xHp * n, n_("e-"): xe * n,
            n_("HD"): 2 * 2.527e-5 * xH2 * n, x_("HD"): 2 * 2.527e-5 * xH2, x_("H_2"): xH2, n_Htot: n, C2: 1.0}
    assert float(H2_cooling_rate().subs(vals)) == pytest.approx(ga08_H2_HD(Tv, n, xH, xH2, xHe, xHp, xe), rel=1e-10)


