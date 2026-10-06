"""Band specs of models.

GIZMO_STARFORGE: GIZMO's STARFORGE bands as they are (one ionizing band, photoelectric, NUV, optical/NIR, IR at the
radiation temperature T_rad). GIZMO's own coefficients for them are overrides in the model (starforge radiation).

STARFORGE_RT: GIZMO's STARFORGE bands with the ionizing band split at the He I threshold. Each ionizing sub-band's
slope is fitted to the stellar spectrum GIZMO assumes for its ionizing band (rt_get_sigma: a 4e4 K blackbody), matching
the H photoionization heat per ionization G_HI and the mean photon energy nu_eff in the sub-band jointly."""

from dataclasses import dataclass

import numpy as np

from .band import Band, BandSet
from .fit import fit_slope, MEAN_PHOTON_ENERGY
from .interaction import Absorber
from .projection import mean_photon_energy, ppl_absorption, ppl_mean_photon_energy, spectrum_absorption
from .spectra import Blackbody, SIGMA_HI, integrate

GIZMO_T_EFF = 4.0e4  # K, rt_get_sigma's T_eff under GALSF
H_PHOTOIONIZATION = Absorber("H photoionization", SIGMA_HI, absorber="H", E_th=13.6)
G_HI = ("H photoionization", "excess")
SIGMA_HI_MEAN = ("H photoionization", "sigma_N")
IONIZING_EDGES = (13.6, 24.59, 500.0)  # eV: H and He I thresholds, GIZMO's upper limit


def ionizing_fits(T_eff=GIZMO_T_EFF, opacity="exact"):
    """SlopeFit per ionizing sub-band (EUV_H, EUV_He): G_HI and nu_eff of a T_eff blackbody matched jointly, sigma_HI
    reported"""
    bands = [Band(n, lo, hi, unit="photons") for n, lo, hi in
             zip(("EUV_H", "EUV_He"), IONIZING_EDGES[:-1], IONIZING_EDGES[1:])]
    return [fit_slope(b, Blackbody(T_eff), [G_HI, MEAN_PHOTON_ENERGY], [H_PHOTOIONIZATION], opacity,
                      report=[SIGMA_HI_MEAN]) for b in bands]


def starforge_rt(T_eff=GIZMO_T_EFF, opacity="exact"):
    """BandSet EUV_H, EUV_He (photons), FUV, NUV, ONIR (energy, E u_E = const), IR (energy, tracked at T_rad)"""
    return BandSet([*(f.band for f in ionizing_fits(T_eff, opacity)),
                    Band("FUV", 8.0, 13.6), Band("NUV", 3.444, 8.0), Band("ONIR", 0.4133, 3.444),
                    Band("IR", 0.001, 0.4133, shape="tracked", temperature="T_rad")])


@dataclass(frozen=True)
class IonizingBandCheck:
    """One photons band against a reference spectrum, per absorption of the absorber: the band's mean photon energy
    hnu (what the band loses per photon in energy) and the mean energy of an absorbed photon E_abs (what the products
    receive: E_th to chemistry, the rest to heat), under the band's spectrum and under the reference."""

    band: str
    hnu: float
    E_abs: float
    hnu_reference: float
    E_abs_reference: float

    @property
    def mismatch(self):
        """(hnu - E_abs) / hnu under the band's spectrum"""
        return 1 - self.E_abs / self.hnu

    @property
    def mismatch_reference(self):
        return 1 - self.E_abs_reference / self.hnu_reference


def hardening_mismatch(band, absorber=H_PHOTOIONIZATION, reference=Blackbody(GIZMO_T_EFF), opacity="exact"):
    """IonizingBandCheck of a ppl band"""
    raw = ppl_absorption(absorber, band.E_lo, band.E_hi, band.slope, opacity)
    ref = spectrum_absorption(absorber, band.E_lo, band.E_hi, reference)
    return IonizingBandCheck(band.name, ppl_mean_photon_energy(band.E_lo, band.E_hi, band.slope), raw["E_abs"],
                             mean_photon_energy(band.E_lo, band.E_hi, reference), ref["E_abs"])


def combined_ionizing_residuals(bands, reference=Blackbody(GIZMO_T_EFF), absorber=H_PHOTOIONIZATION, opacity="exact"):
    """(sigma_HI, G_HI, nu_eff) over the union of bands, each band holding the reference's photons in it, relative to
    the reference over the union (value / reference - 1)"""
    N = np.array([integrate(lambda E: reference(E) / E, b.E_lo, b.E_hi) for b in bands])
    raw = [ppl_absorption(absorber, b.E_lo, b.E_hi, b.slope, opacity) for b in bands]
    sigma = np.array([r["sigma_N"] for r in raw])
    excess = np.array([r["E_abs"] - absorber.E_th for r in raw])
    hnu = np.array([ppl_mean_photon_energy(b.E_lo, b.E_hi, b.slope) for b in bands])
    lo, hi = min(b.E_lo for b in bands), max(b.E_hi for b in bands)
    ref = spectrum_absorption(absorber, lo, hi, reference)
    return (np.sum(N * sigma) / N.sum() / ref["sigma_N"] - 1,
            np.sum(N * sigma * excess) / np.sum(N * sigma) / (ref["E_abs"] - absorber.E_th) - 1,
            np.sum(N * hnu) / N.sum() / mean_photon_energy(lo, hi, reference) - 1)


STARFORGE_RT = starforge_rt()

GIZMO_STARFORGE = BandSet([
    Band("EUV", 13.6, 500.0, unit="photons", doc="ionizing photons (13.6-500 eV) per H nucleus"),
    Band("FUV", 8.0, 13.6, doc="photoelectric band (8-13.6 eV) energy per H nucleus [eV]"),
    Band("NUV", 3.444, 8.0, doc="NUV band (3.444-8 eV) energy per H nucleus [eV]"),
    Band("ONIR", 0.4133, 3.444, doc="optical/NIR band (0.4133-3.444 eV) energy per H nucleus [eV]"),
    Band("IR", 0.001, 0.4133, shape="tracked", temperature="T_rad",
         doc="IR band (0.001-0.4133 eV) energy per H nucleus [eV]"),
])
