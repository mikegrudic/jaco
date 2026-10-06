"""Spectral functions of photon energy E [eV]: power laws, blackbodies, photoionization cross-sections, and the
quadrature the band projections use."""

from dataclasses import dataclass
from typing import Callable

import numpy as np

# CODATA 2018
EV = 1.602176634e-12  # erg
K_B_EV = 8.617333262e-5  # eV / K
H_EV = 4.135667696e-15  # eV s
C_LIGHT = 2.99792458e10  # cm / s
SIGMA_SB = 5.670374419e-5  # erg cm^-2 s^-1 K^-4
MEGABARN = 1e-18  # cm^2


@dataclass(frozen=True)
class PowerLaw:
    """amplitude (E / E_ref)^index for E >= E_min, zero below. Projections of power laws against power-law spectra
    have closed forms."""

    amplitude: float = 1.0
    E_ref: float = 1.0
    index: float = 0.0
    E_min: float = 0.0

    def __call__(self, E, T=None):
        E = np.asarray(E, dtype=float)
        return np.where(E >= self.E_min, self.amplitude * (E / self.E_ref) ** self.index, 0.0)


@dataclass(frozen=True)
class Spectral:
    """A function f(E) (or f(E, T) with temperature_dependent) for E >= E_min, zero below; vectorized in E (and T)"""

    func: Callable
    E_min: float = 0.0
    temperature_dependent: bool = False
    name: str = ""

    def __call__(self, E, T=None):
        E = np.asarray(E, dtype=float)
        value = self.func(E, T) if self.temperature_dependent else self.func(E)
        return np.where(E >= self.E_min, value, 0.0)


def as_spectral(f):
    """A PowerLaw or Spectral from a number (a constant), a callable of E, or either class"""
    if isinstance(f, (PowerLaw, Spectral)):
        return f
    if np.isscalar(f):
        return PowerLaw(float(f), 1.0, 0.0)
    if callable(f):
        return Spectral(f)
    raise TypeError(f"not a spectral function: {f!r}")


def support(f):
    """Lowest energy at which f can be non-zero"""
    return getattr(f, "E_min", 0.0)


def planck(E, T):
    """B_E(T) [erg s^-1 cm^-2 sr^-1 eV^-1], the Planck function per unit photon energy: B_nu / h"""
    E, T = np.asarray(E, dtype=float), np.asarray(T, dtype=float)
    x = E / (K_B_EV * T)
    occupation = np.exp(-x) / -np.expm1(-x)  # 1 / (e^x - 1) without overflow
    return 2.0 * (E * EV) ** 3 / (H_EV * EV * C_LIGHT) ** 2 * occupation / H_EV


@dataclass(frozen=True)
class Blackbody:
    """The shape of a blackbody's energy spectrum at temperature T: u_E proportional to B_E(T)"""

    T: float

    def __call__(self, E, T=None):
        return planck(E, self.T)


def verner96(E_th, E0, sigma0, ya, P, yw=0.0, y0=0.0, y1=0.0):
    """Photoionization cross-section [cm^2] of Verner et al. (1996, ApJ 465, 487), Eq. 1, from its Table 1 parameters
    (sigma0 in Mb, energies in eV)"""
    def sigma(E):
        x = E / E0 - y0
        y = np.sqrt(x * x + y1 * y1)
        F = ((x - 1) ** 2 + yw * yw) * y ** (0.5 * P - 5.5) * (1 + np.sqrt(y / ya)) ** (-P)
        return sigma0 * MEGABARN * F
    return Spectral(sigma, E_min=E_th, name=f"Verner96 (E_th={E_th} eV)")


SIGMA_HI = verner96(13.6, 4.298e-1, 5.475e4, 3.288e1, 2.963)
SIGMA_HEI = verner96(24.59, 1.361e1, 9.492e2, 1.469, 3.188, 2.039, 4.434e-1, 2.136)
SIGMA_HEII = verner96(54.42, 1.720, 1.369e4, 3.288e1, 2.963)


# --- quadrature ----------------------------------------------------------------------------------------------------

_GL_NODES, _GL_WEIGHTS = np.polynomial.legendre.leggauss(8)
PANEL_WIDTH = 0.1  # in ln E


def quadrature_nodes(E_lo, E_hi, panel_width=PANEL_WIDTH):
    """(E, weights) of Gauss-Legendre quadrature in ln E on [E_lo, E_hi]: sum(w * f(E)) = integral of f dE"""
    if not 0 < E_lo < E_hi:
        return np.empty(0), np.empty(0)
    L0, L1 = np.log(E_lo), np.log(E_hi)
    n = max(1, int(np.ceil((L1 - L0) / panel_width)))
    edges = np.linspace(L0, L1, n + 1)
    half = 0.5 * np.diff(edges)
    mid = 0.5 * (edges[1:] + edges[:-1])
    lnE = (mid[:, None] + half[:, None] * _GL_NODES[None, :]).ravel()
    w = (half[:, None] * _GL_WEIGHTS[None, :]).ravel()
    E = np.exp(lnE)
    return E, w * E


def integrate(f, E_lo, E_hi, T=None):
    """integral of f(E[, T]) dE over [E_lo, E_hi]; with T an array, one integral per temperature"""
    E, w = quadrature_nodes(E_lo, E_hi)
    if T is None:
        return float(np.sum(w * f(E))) if E.size else 0.0
    T = np.asarray(T, dtype=float)
    if not E.size:
        return np.zeros_like(T)
    return np.sum(w[:, None] * f(E[:, None], T[None, :]), axis=0)


def power_integral(p, a, b):
    """integral of E^p dE over [a, b], stable at p = -1"""
    if not 0 < a < b:
        return 0.0
    L = np.log(b / a)
    q = p + 1.0
    return a**q * (L if abs(q * L) < 1e-12 else np.expm1(q * L) / q)
