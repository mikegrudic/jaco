"""Projection of spectral declarations onto a band set: band-mean opacities, energies per absorbed photon and emission
fractions, as constants or tables in a temperature, computed once at code-generation time.

For a band [a, c] with energy spectrum u_E and photon spectrum u_E / E, and an absorber sigma(E):

- sigma_N = int sigma u/E dE / int u/E dE    (photon mean: absorptions per absorber = c sigma_N n_photons)
- chi_E   = int sigma u dE / int u dE        (energy mean: absorbed power per absorber = c chi_E u_band)
- chi_F   = int (sigma + sigma_s) u dE / int u dE   (flux mean: the fixed-slope flux spectrum has u's shape)
- E_abs   = int sigma u dE / int sigma u/E dE      (mean energy of an absorbed photon; E_abs - E_th = <E - E_th>_b)

Where sigma and u are both power laws these are the closed forms of He, Wibking & Krumholz (2024), Eq. 30; otherwise
Gauss-Legendre quadrature in ln E.
"""

from dataclasses import dataclass, field, replace
from functools import lru_cache

import numpy as np
import sympy as sp

from .band import BandSet
from .interaction import Line, Continuum, ThermalEmission
from .spectra import PowerLaw, K_B_EV, planck, integrate, power_integral, support

DEFAULT_T_GRID = tuple(np.logspace(0.0, 8.0, 321))  # K, 0.025 dex
OPACITY_MODELS = ("edges", "exact")
TINY = 1e-300


class TTable:
    """A quantity tabulated in a temperature, interpolated linearly in log-log and clamped at the ends"""

    def __init__(self, T, values, name=""):
        self.T = np.asarray(T, dtype=float)
        self.values = np.asarray(values, dtype=float)
        self.name = name
        if self.T.shape != self.values.shape:
            raise ValueError("TTable: T and values differ in shape")

    def __call__(self, T):
        logv = np.log10(np.maximum(self.values, TINY))
        return 10 ** np.interp(np.log10(T), np.log10(self.T), logv)

    def __repr__(self):
        return f"TTable({self.name or '?'}: {self.values[0]:.4g} .. {self.values[-1]:.4g})"

    def expr(self, T, name=None):
        """sympy expression in the temperature symbol T"""
        from ..symbols import piecewise_linear
        if not np.any(self.values > 0):
            return sp.S.Zero
        logv = np.log10(np.maximum(self.values, TINY))
        return 10 ** piecewise_linear(np.log10(self.T), logv, sp.log(T, 10), name=name or self.name or None)


def symbolic(value, T=None, name=None):
    """sympy form of a projected quantity: a Float, a TTable's interpolant in T, or an override as given"""
    if isinstance(value, TTable):
        if T is None:
            raise ValueError("a tabulated quantity needs its temperature symbol")
        return value.expr(T, name)
    return sp.sympify(value)


# --- core integrals ------------------------------------------------------------------------------------------------

def _scaled_planck(E_ref):
    """B_E(T) e^(E_ref / kT): the Planck shape with the exponential at E_ref factored out, for band means that would
    otherwise underflow; the exponent is capped at 700 so that a band many hundred kT wide still has a finite mean"""
    def f(E, T):
        x = E / (K_B_EV * T)
        shift = np.minimum(x - E_ref / (K_B_EV * T), 700.0)
        return 2.0 * E**3 * np.exp(-shift) / -np.expm1(-np.minimum(x, 700.0))
    return f


def _scattering_only(absorber, chi_F, like):
    """Coefficients of an absorber that only scatters in a band"""
    zero = 0.0 * like
    return dict(sigma_N=zero, chi_E=zero, chi_F=chi_F, E_abs=zero + absorber.E_th, heat=zero, remainder=zero)


def _absorption_integrals(absorber, a, c, u, T=None):
    """dict of the absorption means of absorber on [a, c] under the energy spectrum u(E[, T]) by quadrature; None
    where it neither absorbs nor scatters in the band"""
    lo = max(a, absorber.E_min)
    s = absorber.sigma
    U = integrate(lambda E, *t: u(E, *t), a, c, T)
    scat = 0.0 * U
    if absorber.scattering is not None:
        scat = integrate(lambda E, *t: absorber.scattering(E) * u(E, *t),
                         max(a, support(absorber.scattering)), c, T)
    S_N = integrate(lambda E, *t: s(E) * u(E, *t) / E, lo, c, T) if lo < c else 0.0 * U
    if not np.any(S_N > 0):
        return _scattering_only(absorber, scat / U, U) if np.any(scat > 0) else None
    N = integrate(lambda E, *t: u(E, *t) / E, a, c, T)
    S_E = integrate(lambda E, *t: s(E) * u(E, *t), lo, c, T)
    E_abs = S_E / S_N
    excess = E_abs - absorber.E_th
    if absorber.heat_yield is None:
        heat = excess
    else:
        Y, Eth = absorber.heat_yield, absorber.E_th
        heat = integrate(lambda E, *t: s(E) * Y(E) * (E - Eth) * u(E, *t) / E, lo, c, T) / S_N
    return dict(sigma_N=S_N / N, chi_E=S_E / U, chi_F=(S_E + scat) / U, E_abs=E_abs, heat=heat,
                remainder=excess - heat)


def _is_power_law(f):
    return isinstance(f, PowerLaw)


def _power_law_moment(f, p, a, c):
    """integral of f(E) E^p dE over [a, c] for a PowerLaw f"""
    lo = max(a, f.E_min)
    return f.amplitude * f.E_ref ** -f.index * power_integral(p + f.index, lo, c) if lo < c else 0.0


def _absorption_closed_form(absorber, a, c, alpha):
    """_absorption_integrals for a power-law sigma (and constant yield, power-law scattering) under u = E^alpha"""
    lo = max(a, absorber.E_min)
    N, U = power_integral(alpha - 1, a, c), power_integral(alpha, a, c)
    scat = 0.0 if absorber.scattering is None else _power_law_moment(absorber.scattering, alpha, a, c)
    S_N = _power_law_moment(replace(absorber.cross_section, E_min=lo), alpha - 1, a, c) if lo < c else 0.0
    if not S_N > 0:
        return _scattering_only(absorber, scat / U, U) if scat > 0 else None
    S_E = _power_law_moment(replace(absorber.cross_section, E_min=lo), alpha, a, c)
    E_abs = S_E / S_N
    excess = E_abs - absorber.E_th
    heat = excess * (1.0 if absorber.heat_yield is None else absorber.heat_yield.amplitude)
    return dict(sigma_N=S_N / N, chi_E=S_E / U, chi_F=(S_E + scat) / U, E_abs=E_abs, heat=heat,
                remainder=excess - heat)


def _closed_form_applies(absorber):
    y = absorber.heat_yield
    return (_is_power_law(absorber.cross_section) and (y is None or (_is_power_law(y) and y.index == 0))
            and (absorber.scattering is None or _is_power_law(absorber.scattering)))


def edge_power_law(f, lo, hi):
    """The power law through f's values at lo and hi (zero below lo), or None where either is not positive"""
    f_lo, f_hi = float(f(lo)), float(f(hi))
    if not (f_lo > 0 and f_hi > 0):
        return None
    return PowerLaw(f_lo, lo, np.log(f_hi / f_lo) / np.log(hi / lo), E_min=lo)


def _edge_absorber(absorber, a, c):
    """absorber with sigma and scattering replaced by their power laws through the band's edges (each where it is
    positive at both)"""
    changes = {}
    lo = max(a, absorber.E_min)
    if lo < c:
        sigma = edge_power_law(absorber.sigma, lo, c)
        if sigma is not None:
            changes["cross_section"] = sigma
    scat = absorber.scattering
    if scat is not None and max(a, support(scat)) < c:
        scat = edge_power_law(scat, max(a, support(scat)), c)
        if scat is not None:
            changes["scattering"] = scat
    return replace(absorber, **changes) if changes else absorber


@lru_cache(maxsize=None)
def ppl_absorption(absorber, E_lo, E_hi, slope, opacity="edges"):
    """Absorption means on a fixed-slope band, u_E = E^slope (dict of floats, or None)"""
    if opacity not in OPACITY_MODELS:
        raise ValueError(f"opacity model must be one of {OPACITY_MODELS}")
    if opacity == "edges":
        absorber = _edge_absorber(absorber, E_lo, E_hi)
    if _closed_form_applies(absorber):
        return _absorption_closed_form(absorber, E_lo, E_hi, slope)
    return _absorption_integrals(absorber, E_lo, E_hi, PowerLaw(1.0, 1.0, slope))


@lru_cache(maxsize=None)
def tracked_absorption(absorber, E_lo, E_hi, T_grid):
    """Absorption means on a band holding a dilute blackbody at T_rad: dict of arrays over T_grid, or None"""
    T = np.asarray(T_grid)
    return _absorption_integrals(absorber, E_lo, E_hi, _scaled_planck(E_lo), T)


@lru_cache(maxsize=None)
def spectrum_absorption(absorber, E_lo, E_hi, spectrum):
    """Absorption means under a fixed reference spectrum u_E = spectrum(E), exact cross-section"""
    return _absorption_integrals(absorber, E_lo, E_hi, lambda E, *t: spectrum(E))


def mean_photon_energy(E_lo, E_hi, spectrum):
    """int u dE / int u/E dE [eV] of the energy spectrum u = spectrum(E)"""
    return integrate(spectrum, E_lo, E_hi) / integrate(lambda E: spectrum(E) / E, E_lo, E_hi)


def ppl_mean_photon_energy(E_lo, E_hi, slope):
    return power_integral(slope, E_lo, E_hi) / power_integral(slope - 1, E_lo, E_hi)


# --- results -------------------------------------------------------------------------------------------------------

ABSORPTION_FIELDS = ("sigma_N", "chi_E", "chi_F", "E_abs", "heat", "remainder")


@dataclass(frozen=True)
class AbsorptionCoefficients:
    """An absorber's coefficients in one band; floats for a "ppl" band, TTables in T_rad for a "tracked" one.

    sigma_N: photon-mean cross-section (absorptions per absorber per unit time: c sigma_N n_photons)
    chi_E: energy-mean (absorbed power per absorber: c chi_E u_band)
    chi_F: flux-mean extinction, scattering included (momentum)
    E_abs: mean energy of an absorbed photon [eV]
    heat: energy per absorbed photon to heat_to [eV]; remainder: to remainder_to [eV]; E_th to chemistry
    """

    process: str
    band: str
    sigma_N: object
    chi_E: object
    chi_F: object
    E_abs: object
    heat: object
    remainder: object
    E_th: float
    heat_to: str
    remainder_to: str = None
    overridden: frozenset = frozenset()

    @property
    def excess(self):
        """<E - E_th>_b per absorbed photon [eV]"""
        if isinstance(self.E_abs, TTable):
            return TTable(self.E_abs.T, self.E_abs.values - self.E_th, f"excess_{self.process}_{self.band}")
        return self.E_abs - self.E_th

    def expr(self, name, T=None):
        """sympy form of one coefficient (T: the band's radiation temperature symbol, for a tracked band)"""
        return symbolic(getattr(self, name), T, f"{name}_{self.process}_{self.band}")


EMISSION_FIELDS = ("fraction", "photons_per_eV", "planck_mean")


@dataclass(frozen=True)
class EmissionFractions:
    """Where an emitter's energy goes; floats, or TTables in the emitter's temperature.

    fraction: band -> share of the emitted energy
    photons_per_eV: band -> photons emitted into the band per eV emitted in total
    escape: share in no band (1 - sum of fraction, computed independently)
    For ThermalEmission also: planck_mean (band -> int_b kappa B / int_b B), band_planck (band -> int_b B_E dE
    [erg s^-1 cm^-2 sr^-1]) and kappa_planck (int kappa B / int B over all energies).
    """

    process: str
    source: str
    fraction: dict
    photons_per_eV: dict
    escape: object
    planck_mean: dict = None
    band_planck: dict = None
    kappa_planck: object = None
    overridden: frozenset = frozenset()


# --- emission ------------------------------------------------------------------------------------------------------

def _thermal_range(bands, T):
    E_lo = min([1e-5 * K_B_EV * T[0]] + [0.1 * b.E_lo for b in bands])
    E_hi = max([80.0 * K_B_EV * T[-1]] + [b.E_hi for b in bands])
    return E_lo, E_hi


@lru_cache(maxsize=None)
def thermal_emission(emitter, bands, T_grid):
    """EmissionFractions of a ThermalEmission over a BandSet, tables over T_grid"""
    T = np.asarray(T_grid)
    j = lambda E, t: emitter.kappa(E) * planck(E, t)  # noqa: E731
    E_lo, E_hi = _thermal_range(bands, T)
    total = integrate(j, E_lo, E_hi, T)
    B_total = integrate(planck, E_lo, E_hi, T)
    frac, photons, pmean, bplanck = {}, {}, {}, {}
    for b in bands:
        jb = integrate(j, b.E_lo, b.E_hi, T)
        frac[b.name] = TTable(T, jb / total, f"f_{emitter.name}_{b.name}")
        photons[b.name] = TTable(T, integrate(lambda E, t: j(E, t) / E, b.E_lo, b.E_hi, T) / total,
                                 f"photons_{emitter.name}_{b.name}")
        scaled = _scaled_planck(b.E_lo)
        pm = (integrate(lambda E, t: emitter.kappa(E) * scaled(E, t), b.E_lo, b.E_hi, T)
              / integrate(scaled, b.E_lo, b.E_hi, T))
        pmean[b.name] = TTable(T, pm, f"kappaB_{emitter.name}_{b.name}")
        bplanck[b.name] = TTable(T, integrate(planck, b.E_lo, b.E_hi, T), f"B_{b.name}")
    escape = sum((integrate(j, lo, hi, T) for lo, hi in bands.gaps(E_lo, E_hi)), np.zeros_like(T)) / total
    return EmissionFractions(emitter.name, emitter.source, frac, photons, TTable(T, escape, f"escape_{emitter.name}"),
                             pmean, bplanck, TTable(T, total / B_total, f"kappaP_{emitter.name}"))


@lru_cache(maxsize=None)
def continuum_emission(emitter, bands, T_grid):
    """EmissionFractions of a Continuum over a BandSet; tables over T_grid if it depends on temperature"""
    lo, hi = emitter.E_range
    if emitter.temperature_dependent:
        T = np.asarray(T_grid)
        j = emitter.j
        wrap = lambda v, n: TTable(T, v, f"{n}_{emitter.name}")  # noqa: E731
    else:
        T = None
        j = lambda E, *t: emitter.j(E)  # noqa: E731
        wrap = lambda v, n: float(v)  # noqa: E731
    total = integrate(j, lo, hi, T)
    frac, photons = {}, {}
    for b in bands:
        a, c = max(lo, b.E_lo), min(hi, b.E_hi)
        frac[b.name] = wrap(integrate(j, a, c, T) / total, f"f_{b.name}")
        photons[b.name] = wrap(integrate(lambda E, *t: j(E, *t) / E, a, c, T) / total, f"photons_{b.name}")
    escape = sum((integrate(j, a, c, T) for a, c in bands.gaps(lo, hi)), 0.0 if T is None else np.zeros_like(T))
    return EmissionFractions(emitter.name, emitter.source, frac, photons, wrap(escape / total, "escape"))


def line_emission(emitter, bands):
    b = bands.band_at(emitter.E0)
    frac = {x.name: float(x is b) for x in bands}
    photons = {x.name: (1.0 / emitter.E0 if x is b else 0.0) for x in bands}
    return EmissionFractions(emitter.name, emitter.source, frac, photons, 0.0 if b else 1.0)


# --- the projector -------------------------------------------------------------------------------------------------

@dataclass
class Projector:
    """The spectral declarations of a model projected onto a band set.

    Parameters
    ----------
    bands: BandSet
    absorbers: sequence of Absorber
    emitters: sequence of Line, Continuum or ThermalEmission
    overrides: dict, optional
        (process, band) -> {coefficient: value}: values (numbers or sympy expressions) that replace the projected ones,
        e.g. a legacy model's constants. Absorption coefficients are named as in AbsorptionCoefficients; emission ones
        are "fraction", "photons_per_eV" and "planck_mean", and an overridden fraction makes escape 1 - sum(fractions).
        An override may couple an absorber to a band it does not overlap (the other coefficients are then 0).
    opacity: str
        "edges": in a "ppl" band an opacity is the power law through its values at the band's edges (where it is
        positive at both); "exact": the declared function, by quadrature. Tracked bands always use the declared one.
    T_grid: sequence of float
        Temperatures [K] of the tables (tracked bands' T_rad, emitters' temperature).
    """

    bands: BandSet
    absorbers: tuple = ()
    emitters: tuple = ()
    overrides: dict = field(default_factory=dict)
    opacity: str = "edges"
    T_grid: tuple = DEFAULT_T_GRID

    def __post_init__(self):
        if not isinstance(self.bands, BandSet):
            self.bands = BandSet(self.bands)
        self.absorbers = {a.name: a for a in self.absorbers}
        self.emitters = {e.name: e for e in self.emitters}
        if set(self.absorbers) & set(self.emitters):
            raise ValueError(f"processes both absorb and emit: {sorted(set(self.absorbers) & set(self.emitters))}")
        if self.opacity not in OPACITY_MODELS:
            raise ValueError(f"opacity model must be one of {OPACITY_MODELS}")
        self.T_grid = tuple(float(t) for t in self.T_grid)
        self.overrides = {k: dict(v) for k, v in (self.overrides or {}).items()}
        for (process, band), values in self.overrides.items():
            self.bands[band]  # noqa: B018  (raises for an unknown band)
            if process in self.absorbers:
                allowed = ABSORPTION_FIELDS
            elif process in self.emitters:
                allowed = EMISSION_FIELDS
            else:
                raise KeyError(f"override for unknown process {process!r}")
            unknown = set(values) - set(allowed)
            if unknown:
                raise KeyError(f"override ({process}, {band}): unknown coefficients {sorted(unknown)}; "
                               f"allowed {allowed}")
        self._cache = {}

    def mean_photon_energy(self, band):
        """<h nu>_b [eV], number-weighted: a float, or a TTable in T_rad for a tracked band"""
        b = self.bands[band]
        if b.shape == "ppl":
            return ppl_mean_photon_energy(b.E_lo, b.E_hi, b.slope)
        T = np.asarray(self.T_grid)
        f = _scaled_planck(b.E_lo)
        return TTable(T, integrate(f, b.E_lo, b.E_hi, T) / integrate(lambda E, t: f(E, t) / E, b.E_lo, b.E_hi, T),
                      f"hnu_{b.name}")

    def absorption(self, process, band):
        """AbsorptionCoefficients of absorber process in band, or None where it neither absorbs nor scatters there"""
        key = ("abs", process, band)
        if key not in self._cache:
            self._cache[key] = self._absorption(process, band)
        return self._cache[key]

    def _absorption(self, process, band):
        a, b = self.absorbers[process], self.bands[band]
        if b.shape == "ppl":
            raw = ppl_absorption(a, b.E_lo, b.E_hi, b.slope, self.opacity)
        else:
            raw = tracked_absorption(a, b.E_lo, b.E_hi, self.T_grid)
            if raw is not None:
                raw = {k: TTable(self.T_grid, v, f"{k}_{process}_{band}") for k, v in raw.items()}
        over = self.overrides.get((process, band), {})
        if raw is None:
            if not over:
                return None
            raw = dict.fromkeys(ABSORPTION_FIELDS, 0.0)
        raw = {**raw, **over}  # never modify the cached projection
        return AbsorptionCoefficients(process, band, **raw, E_th=a.E_th, heat_to=a.heat_to,
                                      remainder_to=a.remainder_to, overridden=frozenset(over))

    def absorptions(self):
        """(process, band) -> AbsorptionCoefficients for every pair that absorbs"""
        out = {}
        for p in self.absorbers:
            for b in self.bands.names:
                c = self.absorption(p, b)
                if c is not None:
                    out[(p, b)] = c
        return out

    def emission(self, process):
        """EmissionFractions of emitter process"""
        key = ("emit", process)
        if key not in self._cache:
            self._cache[key] = self._emission(process)
        return self._cache[key]

    def _emission(self, process):
        e = self.emitters[process]
        if isinstance(e, Line):
            res = line_emission(e, self.bands)
        elif isinstance(e, Continuum):
            res = continuum_emission(e, self.bands, self.T_grid)
        elif isinstance(e, ThermalEmission):
            res = thermal_emission(e, self.bands, self.T_grid)
        else:
            raise TypeError(f"unknown emitter {e!r}")
        over = {b: v for (p, b), v in self.overrides.items() if p == process}
        if not over:
            return res
        changes, names = {}, set()
        for fld in EMISSION_FIELDS:
            values = dict(getattr(res, fld) or {})
            for b, v in over.items():
                if fld in v:
                    values[b] = v[fld]
                    names.add((b, fld))
            changes[fld] = values if getattr(res, fld) is not None or values else None
        if any(f == "fraction" for _, f in names):
            frac = changes["fraction"]
            if any(isinstance(v, TTable) for v in frac.values()):
                raise ValueError(f"{process}: an overridden fraction needs every band's fraction as a constant")
            changes["escape"] = 1 - sum(frac.values())
        return replace(res, **changes, overridden=frozenset(names))
