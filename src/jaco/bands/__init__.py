"""Radiation bands as a spec: processes declare the photon energies they absorb or emit, and a band set turns those
declarations into per-band coefficients. See docs/band_spec.md."""

from .band import Band, BandSet
from .spectra import (PowerLaw, Spectral, Blackbody, planck, integrate, SIGMA_HI, SIGMA_HEI, SIGMA_HEII, verner96,
                      EV, K_B_EV, C_LIGHT, SIGMA_SB)
from .interaction import Absorber, Line, Continuum, RoutedEmission, ThermalEmission
from .projection import Projector, AbsorptionCoefficients, EmissionFractions, TTable, symbolic, edge_power_law
from .fit import fit_slope, SlopeFit, MEAN_PHOTON_ENERGY
from .structure import band_structure, BandStructure
from .specs import STARFORGE_RT, GIZMO_STARFORGE, starforge_rt, ionizing_fits, hardening_mismatch, combined_ionizing_residuals
