"""Dust opacities as functions of photon energy.

draine_mw31(): the Milky Way R_V = 3.1 carbonaceous-silicate model of Weingartner & Draine (2001, case A), renormalized
as in Draine (2003, ARA&A 41, 241), with Li & Draine (2001) PAH and Draine (2003, ApJ 598, 1017) graphite/silicate
optical properties: Draine's table kext_albedo_WD_MW_3.1_60_D03.all (DRAINE_MW31_URL, computed 2009 Oct 4), stored
unmodified as package data in jaco/bands/data/.

hensley_draine_astrodust(): the astrodust + PAH model of Hensley & Draine (2023, ApJ 948, 55), from the table
hensley_draine_2022_astrodust_opacity.dat (per-H absorption, scattering and extinction cross sections of the PAH and
astrodust populations, 0.1 micron to 3 cm) as distributed with meshoid, stored unmodified in jaco/bands/data/. It stops
at 12.4 eV: above it the model is extrapolated flat, so it is for the non-ionizing bands.
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


ASTRODUST_FILE = "hensley_draine_2022_astrodust_opacity.dat"
ASTRODUST_GAS_PER_H = 1.4  # gas mass per H nucleus in m_H, the convention of Draine's tables (M_gas/M_dust there)


def read_astrodust(path, name=None):
    """DustModel from the Hensley & Draine astrodust table: lambda [micron], then C_abs, C_sca and C_ext per H
    [cm^2/H] of the PAH and astrodust populations, C_ext total and the albedo; the header gives Md/MH, the dust mass
    per H mass. kappa = (C_PAH + C_Ad) / (Md/MH m_H); dust_to_gas = (Md/MH) / ASTRODUST_GAS_PER_H."""
    text = Path(path).read_text()
    m = re.search(r"Md/MH\s*=\s*([0-9.Ee+-]+)", text)
    if not m:
        raise ValueError(f"{path}: no Md/MH in the header")
    md_mh = float(m.group(1))
    data = np.loadtxt(path)
    if data.ndim != 2 or data.shape[1] < 5:
        raise ValueError(f"{path}: expected lambda and four cross-section columns")
    lam, cabs_pah, cabs_ad, csca_pah, csca_ad = data[:, :5].T
    m_dust_per_H = md_mh * 1.6735575e-24  # g
    return DustModel(name or Path(path).name, HC_EV_MICRON / lam, (cabs_pah + cabs_ad) / m_dust_per_H,
                     (csca_pah + csca_ad) / m_dust_per_H, md_mh / ASTRODUST_GAS_PER_H, source=str(path))


def hensley_draine_astrodust():
    """Hensley & Draine's astrodust + PAH dust (module docstring)"""
    path = DATA_DIR / ASTRODUST_FILE
    if not path.is_file():
        raise FileNotFoundError(f"{path} is missing: it is meshoid's radiation/opacity_tables/{ASTRODUST_FILE}")
    return read_astrodust(path, "Hensley & Draine astrodust + PAH (2023)")


def draine_mw31():
    """Draine's Milky Way R_V = 3.1 dust (module docstring)"""
    path = DATA_DIR / DRAINE_MW31_FILE
    if not path.is_file():
        raise FileNotFoundError(f"{path} is missing; fetch it unmodified with: curl -o {path} {DRAINE_MW31_URL}")
    return read_draine(path, "Draine MW R_V=3.1 (WD01 case A, D03)")
