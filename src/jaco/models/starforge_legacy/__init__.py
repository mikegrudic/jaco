"""GIZMO's legacy cooling module as a jaco model, built from the starforge process library.

It reproduces cooling/cooling.cc, eos/eos.cc and cooling/simple_chemistry.cc in the COOL_LOW_TEMPERATURES +
SIMPLE_STEADYSTATE_CHEMISTRY + COOL_MOLECFRAC_NONEQM configuration term for term. Line numbers refer to cooling.cc at
gizmo_jaco_dev e5865f97. What separates it from STARFORGE (jaco.models.starforge):

- Clumping (rule): GIZMO clumps only the H2 formation/dissociation terms update_explicit_molecular_fraction multiplies by
  clumping_factor (2025-2075), which the GIZMO H2 network carries in its rates; every declared clumping factor is dropped.
- Density in the rates (rule): GIZMO's nHcgs = 0.76 rho/m_p in CoolingRate (1046) and simple_chemistry.cc, and rho/m_p
  in update_explicit_molecular_fraction (1985, written into the GIZMO H2 network), which are not the H density n_Htot.
- H2 network: update_explicit_molecular_fraction (1969-2148) transcribed: half-speed rates for its mass-fraction
  variable, n_crit weights from fixed mass fractions, its own H- estimate, no D channels; no heat of H2 formation or
  dissociation.
- Free electrons: find_abundances_and_rates' metal budget (937-948), whose heavy-ion term is its cosmic-ray ionization
  of the neutrals; no explicit cosmic-ray ionization of H.
- C+, [CI] 609 um and CO cooling: Lambda_Metals_Neutral (1158-1176), [CI] weighted by the C+ fraction, HM79 CO with an
  LVG cap; carbon counted as atoms in the EOS and the metal-line tables.
- H2/HD cooling: colliders weighted by the fixed mass fractions X_H, Y_He, HD/H2 = min(0.00126, 4e-5 x_H0 / x_H2)
  (1177-1192).

Everything else is shared: KWH ionization balance and its cooling, gas-dust coupling, cosmic-ray heating, photoelectric
heating, CMB Compton cooling, tabulated metal lines (f_metal), nebular lines (f_neb), the CMB-bath and high-temperature
truncation factors and the PdV term.
"""

import sympy as sp
from jaco.model import Model, Rule
from ..starforge.starforge import (SOLVE_VARS, TIME_DEPENDENT, GIZMO_FAMILY, DERIVED, SHARED_SPECIES, shared_processes,
                                   pdv_work)
from ..starforge.gizmo_lowtemp import gizmo_carbon_cooling, gizmo_H2_cooling
from ..starforge.h2_chemistry.gizmo_network import GizmoH2Network
from ..starforge.metal_electrons import metal_electrons
from ..starforge.symbols import n_Htot, nH_gizmo_cooling, PARAMETERS


def scaled_densities(factor):
    """Expression rule multiplying every number density (n_X, n_Htot) by factor"""
    def rule(e):
        return e.xreplace({s: factor * s for s in e.free_symbols if str(s).startswith("n_")})
    return rule


GIZMO_CLUMPING = Rule("GIZMO clumps only its H2 network", lambda p: p.unclumped())
GIZMO_DENSITY = Rule("GIZMO's nHcgs in the rates", lambda p: p.transformed(scaled_densities(nH_gizmo_cooling / n_Htot)),
                     exempt={"GIZMO H2 network"})


def make_model():
    """The STARFORGE_LEGACY model"""
    intermediates, fixed_electrons = metal_electrons(nH_gizmo_cooling)
    return Model(
        shared_processes() + [gizmo_H2_cooling, gizmo_carbon_cooling, GizmoH2Network(), pdv_work],
        solve_vars=SOLVE_VARS,
        time_dependent=TIME_DEPENDENT,
        fixed={"C+": sp.S.Zero, "CO": sp.S.Zero},  # carbon counted as atoms in the EOS and metal-line tables
        derived=DERIVED,
        intermediates=intermediates,
        fixed_electrons=fixed_electrons,
        rules=[GIZMO_CLUMPING, GIZMO_DENSITY],
        parameters=PARAMETERS,
        species=SHARED_SPECIES,
    )
