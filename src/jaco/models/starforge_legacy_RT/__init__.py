"""GIZMO's legacy cooling module with its matter-radiation coupling (RADTRANSFER + RT_CHEM_PHOTOION), as a jaco model.

starforge_legacy plus the terms the legacy module adds when the ionizing band is evolved by the RT solver. Line numbers
refer to cooling/cooling.cc and radiation/rt_chem.cc at gizmo_jaco_dev bd3500d7.

- Photoionization of H by the ionizing band (find_abundances_and_rates, 846-904): Gamma = c sigma_HI n_gamma per
  neutral H at the true c, with heat eps_HI per ionization (Heat_Ion_from_RHD, 1225-1266). eps_HI is rt_ion_G_HI, the
  cross-section-weighted photoelectron energy of the band, not h nu_eff - 13.6 eV (rt_get_sigma). He is not
  photoionized (rt_ion_G_HeI = rt_ion_G_HeII = 0 without RT_CHEM_PHOTOION_HE). The rate acts on atomic H; GIZMO's
  neutral H includes H2, which the solver's H budget cannot give up to H+ directly.
- Rate law: Gamma_eff = Gamma S(tau), S(x) = (1 - e^-x)/x in GIZMO's Pade form (slab_averaging_function), with
  tau = c_tilde sigma_HI n_H0 Delta_t: the photon-limited average over the step, under which the photons a cell takes,
  n_H0 Gamma_eff Delta_t (c_tilde/c), never exceed n_gamma. c_tilde = 0 gives GIZMO's law, Gamma frozen over the step,
  which is right where the RT solver's kick has already attenuated the band (cooling.cc 875 and 1244 comment S out).
- H+ is time-dependent, as GIZMO's linearized backward-Euler H balance (903) is. GIZMO integrates abundances per
  nHcgs = 0.76 rho/m_p, at which its rates are evaluated (rule GIZMO_DENSITY); rule GIZMO_ABUNDANCES divides the species
  rows by nHcgs/n_Htot so that the backward-Euler term, which jaco writes per n_Htot, runs at GIZMO's speed.
- The emission corrections (the CMB-bath factor of metal-line, molecular and fine-structure cooling) take the
  background radiation temperature T_bg (get_background_radiation_temperature_for_emission_corrections, 2467-2482): the
  energy-weighted IR + CMB temperature under RT_INFRARED. Compton cooling keeps the CMB alone; GIZMO adds the RT bands
  to it, a term negligible against photoheating and the dust coupling.
"""

import sympy as sp
from jaco.model import Rule
from jaco.declarations import Parameter
from jaco.process import Process
from jaco.processes import Reaction
from jaco.symbols import n_, dt
from ..starforge.starforge import GIZMO_FAMILY, SHARED_SPECIES  # noqa: F401  (GIZMO_FAMILY: the host-code family)
from ..starforge.symbols import z, nH_gizmo_cooling, n_Htot
from ..starforge_legacy import make_model as make_legacy

T_CMB_Z0 = 2.73  # GIZMO's CMB temperature at z = 0, as jaco's starforge symbols write it: 2.73 (1 + z)

Gamma_HI = sp.Symbol("Gamma_HI")
sigma_HI = sp.Symbol("sigma_HI")
eps_HI = sp.Symbol("eps_HI")
c_tilde = sp.Symbol("c_tilde")
T_bg = sp.Symbol("T_bg")

RT_PARAMETERS = [
    Parameter("Gamma_HI", "s^-1", 0.0, "photoionization rate per neutral H of the ionizing band at the true c, "
                                       "c sigma_HI n_gamma (cooling.cc gJH0ne n_e)"),
    Parameter("sigma_HI", "cm^2", 0.0, "band-averaged H photoionization cross-section (rt_ion_sigma_HI)"),
    Parameter("eps_HI", "erg", 0.0, "heat per photoionization (rt_ion_G_HI)"),
    Parameter("c_tilde", "cm s^-1", 0.0, "speed of light of the band's absorption in the rate law; 0: Gamma frozen "
                                         "over the step, as where the RT solver absorbs the band"),
    Parameter("T_bg", "K", T_CMB_Z0, "background radiation temperature of the emission corrections"),
]


def slab_average(x):
    """(1 - e^-x)/x: GIZMO's slab_averaging_function, a Pade fit accurate to ~0.1% that is 1 at x = 0 and 1/x at x >> 1"""
    return ((1 + x * (0.21772719088733913 + x * (0.047076512011644776 + x * 0.005068307557496351)))
            / (1 + x * (0.71772719088733920 + x * (0.239273440788647680 + x * (0.046750496137263675
                                                                             + x * 0.005068307557496351)))))


def photoionization_optical_depth():
    """c_tilde sigma_HI n_H0 Delta_t, the band's absorption optical depth over the step"""
    return c_tilde * sigma_HI * n_("H") * dt


photoionization = Reaction("H -> H+ + e-", rate=Gamma_HI * slab_average(photoionization_optical_depth()) * n_("H"),
                           heat_per_reaction=eps_HI, clumping=1, name="Photoionization of H by the RT band",
                           bibliography=["GIZMO cooling.cc find_abundances_and_rates and Heat_Ion_from_RHD",
                                         "GIZMO rt_chem.cc rt_get_sigma"])


def _material_rows_scaled(species_rows, factor):
    """Rule: every row of a species in species_rows divided by factor, the heat and energy rows kept"""
    def rule(p):
        rows = {k: e.rhs for k, e in p.network.items()}
        if not any(k in species_rows and rows[k] != 0 for k in rows):
            return p
        return Process(p.name, p.bibliography, {k: (v / factor if k in species_rows else v) for k, v in rows.items()})
    return rule


def _cmb_to_background(expr):
    """expr with the CMB temperature 2.73 (1 + z) replaced by T_bg wherever it appears as a sum"""
    def fix(a):
        terms = set(a.args)
        for s in (1, -1):
            pair = {s * sp.Float(T_CMB_Z0) * z, s * sp.Float(T_CMB_Z0)}
            if pair <= terms:
                return sp.Add(*(terms - pair)) + s * T_bg
        return a
    return expr.replace(lambda a: a.is_Add, fix)


GIZMO_ABUNDANCES = Rule("GIZMO integrates abundances per nHcgs",
                        _material_rows_scaled({s.name for s in SHARED_SPECIES if s.kind in ("material", "trace")},
                                              nH_gizmo_cooling / n_Htot),
                        exempt={"GIZMO H2 network"})
GIZMO_EMISSION_BACKGROUND = Rule("GIZMO's background temperature in the emission corrections",
                                 lambda p: p.transformed(_cmb_to_background),
                                 exempt={"Inverse Compton cooling (CMB)"})


def make_model():
    """The STARFORGE_LEGACY_RT model"""
    base = make_legacy()
    return base.evolve(
        list(base.processes.values()) + [photoionization],
        time_dependent=("T", "H+", "H_2"),
        rules=[*base.rules, GIZMO_ABUNDANCES, GIZMO_EMISSION_BACKGROUND],
        parameters=[*base.parameters, *RT_PARAMETERS],
    )
