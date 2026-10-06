"""Fitting a fixed-slope band's spectral index to a reference spectrum by matching chosen band means."""

from dataclasses import dataclass, replace

import numpy as np
from scipy.optimize import brentq, minimize_scalar

from .band import Band
from .projection import (ABSORPTION_FIELDS, mean_photon_energy, ppl_absorption, ppl_mean_photon_energy,
                         spectrum_absorption)

MEAN_PHOTON_ENERGY = "mean_photon_energy"


def _value(target, band, absorbers, opacity, spectrum=None):
    """A target's value on band: under its power law, or under spectrum (with the exact cross-sections) if given.

    target: "mean_photon_energy", or (absorber name, coefficient) with coefficient one of ABSORPTION_FIELDS or
    "excess" (E_abs - E_th)
    """
    if target == MEAN_PHOTON_ENERGY:
        if spectrum is None:
            return ppl_mean_photon_energy(band.E_lo, band.E_hi, band.slope)
        return mean_photon_energy(band.E_lo, band.E_hi, spectrum)
    name, coeff = target
    a = absorbers[name]
    if spectrum is None:
        raw = ppl_absorption(a, band.E_lo, band.E_hi, band.slope, opacity)
    else:
        raw = spectrum_absorption(a, band.E_lo, band.E_hi, spectrum)
    if raw is None:
        raise ValueError(f"{name} absorbs nothing in band {band.name}")
    return raw["E_abs"] - a.E_th if coeff == "excess" else raw[coeff]


@dataclass(frozen=True)
class SlopeFit:
    """band: the band with the fitted slope; rows: target -> (reference, fitted, fitted / reference - 1), the fitted
    targets first, then the reported ones"""

    band: Band
    fitted: tuple
    rows: dict

    @property
    def slope(self):
        return self.band.slope

    def table(self):
        lines = [f"{self.band.name} [{self.band.E_lo:g}, {self.band.E_hi:g}] eV: slope {self.slope:.4f}"]
        for t, (ref, mod, res) in self.rows.items():
            tag = "fitted" if t in self.fitted else "reported"
            lines.append(f"  {str(t):45s} reference {ref:.5g}  band {mod:.5g}  residual {res:+.2%}  ({tag})")
        return "\n".join(lines)


def fit_slope(band, reference, targets, absorbers=(), opacity="edges", weights=None, report=(),
              bounds=(-40.0, 10.0)):
    """The slope of band's power-law spectrum that best reproduces the targets of the reference spectrum.

    With one target the slope matches it exactly (where one exists in bounds); with several it minimizes
    sum(weight * ln(fitted / reference)^2).

    Parameters
    ----------
    band: Band ("ppl")
    reference: callable of E
        The reference energy spectrum u_E (e.g. spectra.Blackbody(4e4)), any normalization.
    targets: sequence
        "mean_photon_energy" or (absorber name, coefficient), coefficient in ABSORPTION_FIELDS or "excess".
    absorbers: sequence of Absorber
    opacity: str
        The band's opacity model (Projector.opacity); the reference always uses the exact cross-sections.
    weights: sequence of float, optional
    report: sequence
        Further targets evaluated at the fitted slope, not fitted.
    """
    if band.shape != "ppl":
        raise ValueError("only a ppl band has a slope to fit")
    absorbers = {a.name: a for a in absorbers}
    for t in list(targets) + list(report):
        if t != MEAN_PHOTON_ENERGY and t[1] not in ABSORPTION_FIELDS + ("excess",):
            raise ValueError(f"unknown target {t!r}")
    targets = tuple(t if isinstance(t, str) else tuple(t) for t in targets)
    weights = np.ones(len(targets)) if weights is None else np.asarray(weights, dtype=float)
    ref = {t: _value(t, band, absorbers, opacity, reference) for t in targets}

    def log_residuals(slope):
        b = band.with_slope(slope)
        return np.array([np.log(_value(t, b, absorbers, opacity) / ref[t]) for t in targets])

    if len(targets) == 1:
        lo, hi = bounds
        if np.sign(log_residuals(lo)[0]) == np.sign(log_residuals(hi)[0]):
            raise ValueError(f"no slope in {bounds} matches {targets[0]}")
        slope = brentq(lambda s: log_residuals(s)[0], lo, hi, xtol=1e-12)
    else:
        slope = minimize_scalar(lambda s: np.sum(weights * log_residuals(s) ** 2), bounds=bounds, method="bounded",
                                options={"xatol": 1e-10}).x
    fitted_band = replace(band, slope=float(slope))
    rows = {}
    for t in targets + tuple(t if isinstance(t, str) else tuple(t) for t in report):
        r = ref[t] if t in ref else _value(t, band, absorbers, opacity, reference)
        m = _value(t, fitted_band, absorbers, opacity)
        rows[t] = (r, m, m / r - 1)
    return SlopeFit(fitted_band, targets, rows)
