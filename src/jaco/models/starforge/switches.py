"""The switches that separate the two models built from the starforge process library.

``STARFORGE`` is the physically preferred model. ``STARFORGE_LEGACY`` reproduces GIZMO's legacy cooling module
(cooling/cooling.cc, eos/eos.cc and cooling/simple_chemistry.cc in the COOL_LOW_TEMPERATURES +
SIMPLE_STEADYSTATE_CHEMISTRY + COOL_MOLECFRAC_NONEQM configuration) term for term. Line numbers refer to cooling.cc
at gizmo_jaco_dev e5865f97. Everything not listed here is shared: KWH ionization balance and its cooling, gas-dust
coupling, cosmic-ray heating, photoelectric heating, CMB Compton cooling, tabulated metal lines (f_metal), nebular
lines (f_neb), the CMB-bath and high-temperature truncation factors and the PdV term.
"""

from dataclasses import dataclass


@dataclass(frozen=True)
class Switches:
    # Sub-grid clumping: "all" puts C_2 on every two-body rate and C_3 on three-body ones; "h2_chemistry" keeps it only
    # on the H2 formation/dissociation terms update_explicit_molecular_fraction multiplies by clumping_factor (2025-2075)
    clumping: str
    # Density in the rates: "true" uses n_Htot; "gizmo" uses GIZMO's nHcgs = 0.76 rho/m_p in CoolingRate (1046) and
    # simple_chemistry.cc, and rho/m_p in update_explicit_molecular_fraction (1985), which are not the H density
    rate_density: str
    # H2 network: "per_molecule" is the reaction network with H- in equilibrium coupled to the ions; "gizmo" transcribes
    # update_explicit_molecular_fraction (1969-2148): half-speed rates for its mass-fraction variable, n_crit weights
    # from fixed mass fractions, its own H- estimate, no D channels
    h2_network: str
    # Heat of H2 formation (Hollenbach & McKee 1979 / Omukai 2000 critical-density partition) and of collisional
    # dissociation; GIZMO's cooling module carries neither
    h2_chemical_heat: bool
    # Free electrons beyond the solved H and He ions, and where cosmic-ray ionization of H goes: "solved" ionizes H to
    # H+ explicitly, with radiative, grain-assisted and charge-transfer (to Mg) sinks, and balances C+, Mg+ and
    # molecular ions; "gizmo" is find_abundances_and_rates' metal budget (937-948), whose heavy-ion term is its CR
    # ionization of the neutrals
    electrons: str
    # C+, [CI] 609 um and CO cooling: "network" uses the model's carbon species ([CI] on neutral C, Whitworth & Jaffa
    # 2018 CO); "gizmo" is Lambda_Metals_Neutral (1158-1176): [CI] weighted by the C+ fraction, HM79 CO with an LVG cap
    carbon_cooling: str
    # H2/HD cooling: "network" weights colliders by number and locks all D in HD; "gizmo" weights them by the fixed
    # mass fractions X_H, Y_He and takes HD/H2 = min(0.00126, 4e-5 x_H0 / x_H2) (1177-1192)
    h2_cooling: str


STARFORGE = Switches(
    clumping="all",
    rate_density="true",
    h2_network="per_molecule",
    h2_chemical_heat=True,
    electrons="solved",
    carbon_cooling="network",
    h2_cooling="network",
)

STARFORGE_LEGACY = Switches(
    clumping="h2_chemistry",
    rate_density="gizmo",
    h2_network="gizmo",
    h2_chemical_heat=False,
    electrons="gizmo",
    carbon_cooling="gizmo",
    h2_cooling="gizmo",
)
