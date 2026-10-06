"""Dust opacities as functions of photon energy.

draine_mw31(): the Milky Way R_V = 3.1 carbonaceous-silicate model of Weingartner & Draine (2001, case A), renormalized
as in Draine (2003, ARA&A 41, 241), with Li & Draine (2001) PAH and Draine (2003, ApJ 598, 1017) graphite/silicate
optical properties: Draine's table kext_albedo_WD_MW_3.1_60_D03.all (DRAINE_MW31_URL, computed 2009 Oct 4), stored
unmodified as package data in jaco/bands/data/.
"""

import re
from dataclasses import dataclass
from functools import cached_property
from pathlib import Path

import numpy as np

from .interaction import Absorber
from .spectra import Spectral

HC_EV_MICRON = 1.239841984  # h c [eV micron]
DATA_DIR = Path(__file__).with_name("data")
DRAINE_MW31_FILE = "kext_albedo_WD_MW_3.1_60_D03.all"
DRAINE_MW31_URL = "https://www.astro.princeton.edu/~draine/dust/extcurvs/" + DRAINE_MW31_FILE


@dataclass(eq=False)
class DustModel:
    """A dust model's opacities per gram of dust, tabulated in photon energy.

    E: photon energies [eV], increasing; kappa_abs, kappa_sca: absorption and scattering cross-sections per gram of
    dust [cm^2 g^-1] at E; dust_to_gas: the model's dust-to-gas mass ratio. Between the tabulated energies the
    opacities are interpolated in log-log; beyond them they are the end values.
    """

    name: str
    E: np.ndarray
    kappa_abs: np.ndarray
    kappa_sca: np.ndarray
    dust_to_gas: float
    source: str = ""

    def __post_init__(self):
        order = np.argsort(self.E)
        self.E, self.kappa_abs, self.kappa_sca = (np.asarray(x, dtype=float)[order]
                                                  for x in (self.E, self.kappa_abs, self.kappa_sca))
        if np.any(np.diff(self.E) <= 0) or np.any(self.kappa_abs <= 0) or np.any(self.kappa_sca < 0):
            raise ValueError(f"{self.name}: energies must be distinct and opacities positive")

    def _loglog(self, values):
        lnE, lnv = np.log(self.E), np.log(np.maximum(values, 1e-300))

        def f(E):
            return np.exp(np.interp(np.log(E), lnE, lnv))
        return f

    @cached_property
    def absorption(self):
        """kappa_abs(E) per gram of dust [cm^2 g^-1]"""
        return Spectral(self._loglog(self.kappa_abs), name=f"{self.name} absorption")

    @cached_property
    def scattering(self):
        """kappa_sca(E) per gram of dust [cm^2 g^-1]"""
        return Spectral(self._loglog(self.kappa_sca), name=f"{self.name} scattering")

    def absorber(self, name="dust", per="gas", scale=1, scattering_scale=None, heat_to="dust heat"):
        """An Absorber of the dust: per gram of gas at the model's dust-to-gas ratio (per="gas"), or per gram of dust
        (per="dust"); scale the state-dependent factor (e.g. Z_d f_survival(T_d))"""
        if per not in ("gas", "dust"):
            raise ValueError("per must be 'gas' or 'dust'")
        if per == "dust":
            return Absorber(name, self.absorption, heat_to=heat_to, scattering=self.scattering, scale=scale,
                            scattering_scale=scattering_scale)
        d2g, abs_, sca = self.dust_to_gas, self.absorption, self.scattering
        return Absorber(name, Spectral(lambda E: d2g * abs_(E), name=f"{self.name} absorption per gas"),
                        heat_to=heat_to, scattering=Spectral(lambda E: d2g * sca(E), name=f"{self.name} scattering"),
                        scale=scale, scattering_scale=scattering_scale)


def band_means(dust, bands, per="gas", **projector_options):
    """band name -> (chi_E, chi_S, chi_F) of a DustModel's absorber in each band (floats, or TTables in T_rad)"""
    from .projection import Projector
    proj = Projector(bands, [dust.absorber(per=per)], **projector_options)
    out = {}
    for b in bands.names:
        c = proj.absorption("dust", b)
        out[b] = (c.chi_E, c.chi_S, c.chi_F)
    return out


def read_draine(path, name=None):
    """DustModel from one of Draine's kext_albedo_*.all tables: lambda [micron], albedo, <cos>, C_ext/H [cm^2/H],
    K_abs [cm^2/g dust], <cos^2>, and the header's M_dust per H nucleon and M_gas/M_dust.

    kappa_abs is the K_abs column; kappa_sca = albedo C_ext/H / (M_dust/H)."""
    text = Path(path).read_text()
    m_dust = re.search(r"^\s*([0-9.Ee+-]+)\s*=\s*M_dust per H nucleon", text, re.M)
    gas_to_dust = re.search(r"^\s*([0-9.Ee+-]+)\s*=\s*M_gas/M_dust", text, re.M)
    if not (m_dust and gas_to_dust):
        raise ValueError(f"{path}: no M_dust per H nucleon / M_gas/M_dust in the header")
    rows = []
    for line in text.splitlines():
        tokens = line.split()
        try:
            rows.append([float(t) for t in tokens[:6]])
        except ValueError:
            continue
        if len(rows[-1]) < 6:
            rows.pop()
    data = np.array(rows)
    if data.ndim != 2 or len(data) < 2:
        raise ValueError(f"{path}: no data rows")
    lam, albedo, _, cext_H, k_abs, _ = data.T
    return DustModel(name or Path(path).name, HC_EV_MICRON / lam, k_abs, albedo * cext_H / float(m_dust.group(1)),
                     1.0 / float(gas_to_dust.group(1)), source=str(path))


def draine_mw31():
    """Draine's Milky Way R_V = 3.1 dust (module docstring)"""
    path = DATA_DIR / DRAINE_MW31_FILE
    if not path.is_file():
        raise FileNotFoundError(f"{path} is missing; fetch it unmodified with: curl -o {path} {DRAINE_MW31_URL}")
    return read_draine(path, "Draine MW R_V=3.1 (WD01 case A, D03)")
