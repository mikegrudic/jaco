"""Convergence of band-projected coupling with the number of bands, after He, Wibking & Krumholz (2024, §4.1.2): an
absorber chi ~ E^-2 and its Kirchhoff emission, on N bands log-spaced over 0.01-10 eV (fixed slope E u_E = const),
against a 512-band reference.

1. Thermal relaxation: gas at 3000 K with heat capacity a T0^3 and no radiation, in a closed cell. The gas temperature
   trajectory converges at second order in the band width. Emission into a band at its energy mean (chi_B,b :=
   chi_E,b, Kirchhoff per band, the default and the paper's choice) or at its exact Planck mean: the first reaches
   the exact equilibrium, the second has half the transient error but an equilibrium bias at coarse bands.
2. Absorption of a 6000 K blackbody through a column of optical depth 1 at 1 eV.

Run with -s to print the tables."""

import numpy as np
import pytest

from jaco.bands import Absorber, Band, BandSet, PowerLaw, Projector, ThermalEmission, integrate, planck
from jaco.bands.spectra import C_LIGHT, SIGMA_SB

KAPPA = PowerLaw(1.0, 1.0, -2.0)  # chi / chi(1 eV)
E_LO, E_HI = 0.01, 10.0  # eV, the resolved range
FLOOR, CEILING = 1e-12, 100.0  # outer bands close the cell: below 1e-12 eV a fraction ~1e-11 of the emission escapes
T0 = 3000.0
A_RAD = 4 * SIGMA_SB / C_LIGHT
C_V = A_RAD * T0**3  # erg cm^-3 K^-1
T_GRID = tuple(np.logspace(np.log10(1500.0), np.log10(3200.0), 341))
TAUS = np.logspace(-2, 3, 26)  # c chi(1 eV) t
REFERENCE_BANDS = 512


def band_set(n):
    """n bands log-spaced over [E_LO, E_HI], plus one below and one above"""
    e = np.concatenate([[FLOOR], np.geomspace(E_LO, E_HI, n + 1), [CEILING]])
    return BandSet([Band(f"b{i}", lo, hi) for i, (lo, hi) in enumerate(zip(e[:-1], e[1:]))])


def relax(n, kirchhoff="band", steps_per_decade=250):
    """Gas temperature at TAUS. Bands: dU_b/dtau = s_b(T) - chi_E,b U_b with emission s_b = 4 pi/c chi_B,b B_b(T);
    gas: C_v dT/dtau = -sum_b dU_b/dtau - escape(T). Two-stage L-stable SDIRK; each stage needs only a scalar Newton
    solve in T, the bands coupling through the gas alone"""
    proj = Projector(band_set(n), [Absorber("absorption", KAPPA)], [ThermalEmission("emission", KAPPA, "heat")],
                     T_grid=T_GRID, kirchhoff=kirchhoff)
    names = proj.bands.names
    chi = np.array([proj.absorption("absorption", b).chi_E for b in names])
    e = proj.emission("emission")
    chi_B = np.array([e.planck_mean[b].values for b in names])
    source = 4 * np.pi / C_LIGHT * chi_B * np.array([e.band_planck[b].values for b in names])
    escape = 4 * np.pi / C_LIGHT * e.escape.values * e.kappa_planck.values * SIGMA_SB * np.asarray(T_GRID) ** 4 / np.pi
    logT = np.log10(T_GRID)
    L, Le = np.log10(np.maximum(source, 1e-300)), np.log10(np.maximum(escape, 1e-300))

    def sources(T):
        """band sources, escape and their T derivatives, log-log interpolated"""
        x = np.log10(T)
        i = min(max(np.searchsorted(logT, x) - 1, 0), len(logT) - 2)
        w, h = (x - logT[i]) / (logT[i + 1] - logT[i]), logT[i + 1] - logT[i]
        s = 10 ** (L[:, i] + w * (L[:, i + 1] - L[:, i]))
        q = 10 ** (Le[i] + w * (Le[i + 1] - Le[i]))
        return s, s * (L[:, i + 1] - L[:, i]) / h / T, q, q * (Le[i + 1] - Le[i]) / h / T

    def stage(U_base, T_base, hg, T):
        """Solve U = U_base + hg f_U(U, T), T = T_base + hg f_T(U, T)"""
        d = 1 + hg * chi
        for _ in range(50):
            s, ds, q, dq = sources(T)
            g = C_V * (T - T_base) + hg * (np.sum((s - chi * U_base) / d) + q)
            step = g / (C_V + hg * (np.sum(ds / d) + dq))
            T -= step
            if abs(step) < 1e-13 * T:
                break
        return (U_base + hg * sources(T)[0]) / d, T

    gamma = 1 - np.sqrt(0.5)
    times = np.unique(np.concatenate([[0.0], np.logspace(-6, np.log10(TAUS[-1]), int(steps_per_decade * 9) + 1), TAUS]))
    times = times[np.concatenate([[True], np.diff(times) > 1e-9 * times[1:]])]
    U, T, out = np.zeros(len(chi)), T0, []
    for t0, t1 in zip(times[:-1], times[1:]):
        h = t1 - t0
        U1, T1 = stage(U, T, h * gamma, T)
        fU, fT = (U1 - U) / (h * gamma), (T1 - T) / (h * gamma)
        U, T = stage(U + h * (1 - gamma) * fU, T + h * (1 - gamma) * fT, h * gamma, T1)
        if len(out) < len(TAUS) and np.isclose(t1, TAUS[len(out)], rtol=1e-8, atol=0):
            out.append(T)
    assert len(out) == len(TAUS)
    return np.array(out)


@pytest.fixture(scope="module")
def reference():
    return relax(REFERENCE_BANDS)


def test_reference_is_converged(reference):
    T_eq = reference[-1]
    x = T_eq / T0
    assert x + x**4 == pytest.approx(1, rel=1e-5)  # C_v T + a T^4 conserved, radiation Planckian at the end
    scale = T0 - T_eq
    assert np.max(np.abs(relax(REFERENCE_BANDS // 2) - reference)) / scale < 1e-4
    assert np.max(np.abs(relax(REFERENCE_BANDS, steps_per_decade=500) - reference)) / scale < 1e-5


def test_relaxation_converges(reference):
    scale = T0 - reference[-1]
    rows = {}
    for kirchhoff in ("band", "planck"):
        for n in (4, 8, 16):
            err = np.abs(relax(n, kirchhoff) - reference) / scale
            rows[kirchhoff, n] = (err.max(), TAUS[err.argmax()], err[-1])
    print("\nrelaxation, |T - T_ref| / (T0 - T_eq): max over the trajectory (at tau), and at equilibrium")
    for (k, n), (m, at, final) in rows.items():
        print(f"  emission at chi_B={'Planck mean' if k == 'planck' else 'chi_E,b':12s} {n:3d} bands: "
              f"max {m:.2e} (tau {at:.2g})  equilibrium {final:.1e}")
    for k in ("planck", "band"):
        m = [rows[k, n][0] for n in (4, 8, 16)]
        assert m[0] > m[1] > m[2] and m[1] / m[2] > 3.0  # second order from 8 bands
    assert rows["planck", 16][0] < 0.01 and rows["band", 16][0] < 0.02
    assert rows["planck", 4][2] > 0.01 and rows["planck", 16][2] < 1e-5  # equilibrium bias of the Planck-mean emission
    assert all(rows["band", n][2] < 1e-4 for n in (4, 8, 16))  # Kirchhoff per band: exact equilibrium


def test_absorbed_fraction_converges():
    """A 6000 K blackbody through a column of chi(1 eV) Sigma = 1: the absorbed fraction with each band's energy-mean
    opacity (fixed slope) against the exact integral"""
    T_star, column = 6000.0, 1.0
    u = lambda E: planck(E, T_star)  # noqa: E731
    exact = integrate(lambda E: u(E) * -np.expm1(-KAPPA(E) * column), FLOOR, CEILING) / integrate(u, FLOOR, CEILING)
    print(f"\nabsorbed fraction of a {T_star:g} K blackbody, exact {exact:.5f}")
    errors = []
    for n in (4, 8, 16, 32):
        proj = Projector(band_set(n), [Absorber("absorption", KAPPA)])
        U = np.array([integrate(u, b.E_lo, b.E_hi) for b in proj.bands])
        chi = np.array([proj.absorption("absorption", b).chi_E for b in proj.bands.names])
        f = np.sum(U * -np.expm1(-chi * column)) / U.sum()
        errors.append(abs(f / exact - 1))
        print(f"  {n:3d} bands: {f:.5f}  relative error {f / exact - 1:+.2e}")
    assert errors[0] > errors[1] > errors[2] > errors[3] and errors[2] / errors[3] > 3.0 and errors[3] < 0.01
