"""Menon et al. 2025's (arXiv:2610.17577) thermochemistry as their one-zone collapse test runs it (KROME network
react_popsicle_semenov_photo_cr_full with the driver tests/popsicle_semenov_photo_cr_full_collapseTest), for the
free-fall comparison with their Fig. 2. Not a GIZMO model: its photo- and cosmic-ray inputs (chi_0, Gamma_ion0,
zeta_CR, lambda_J) are parameters GIZMO's host does not fill; closures.MenonClosures supplies them.

Taken from starforge: Lyman-alpha and He+ line cooling, collisional ionization, He recombination, free-free, Compton,
Kim+23 nebular cooling (on), the [CII]/[CI]/[OI] lines (with KROME's DESPOTIC optical-depth factors), CO cooling, the
C+/CO partition, grain-assisted H+ recombination, the free electrons on metals and molecular ions (with KROME's
unattenuated cosmic-ray rate), H- collisional detachment, the HM79 grain formation rate. KROME's (processes.py):
- photochemistry: H photoionization by the driver's ionizing field 2.19e-12 chi exp(-1.14e4 A_V) s^-1 (3.4 eV heat),
  A_V = Z N_H / 1.87e21, so unattenuated at Z = 0; LW dissociation 5.6e-11 chi exp(-2.805 A_V) with WG11 H2 (b = 1
  km/s) and H Lyman-series shielding (0.4 + 16 f eV heat); H- photodetachment 7.1e-7 chi exp(-0.5 A_V);
- cosmic rays: zeta = 3e-17 s^-1 unattenuated; H at (0.5 + 3.731 x_H2) zeta, He at 0.55 zeta, Bialy & Sternberg
  heat; H2 dissociation 0.7 zeta;
- the H-/H2+/H3+/H2 rates (de Jong, Kreckel, Stenrup, Poulaert, Savin, Hunter+23 collisional dissociation by H and
  H2, Forrey three-body with H and H2/8), KROME's chemical heat (3.53, 4.48 and 0.2 + 4.2 f eV per H2 formed, f the
  HM79 de-excitation fraction; 4.48 eV per collisional dissociation), H+ recombination (Abel+97, Cen cooling);
- WD01 photoelectric heating on the attenuated chi_FUV;
- H2 cooling: GA08 + Glover 2015 with Omukai's escape probability; HD: Lipovka 2005, HD from D fractionation
  equilibrium (D/H = 3e-5); gas continuum (Lenzuni);
- dust: Semenov (rho, Td) Planck mean with escape over the Jeans length, gas-dust coupling alpha_gd, FUV and CMB
  absorption.
Not carried (they need the O/C network): OH and H2O cooling, explicit C/O chemistry. No clumping, no tabulated metal
lines, no optically-thick cap.
"""

import sympy as sp
from jaco.model import Model, Rule
from jaco.process import Process
from jaco.declarations import Parameter, Species
from jaco.processes import CollisionalIonization, GasPhaseRecombination, FreeFreeEmission, ThermalTerm
from ..starforge.starforge import SHARED_SPECIES, carbon_abundances
from ..starforge.line_cooling import LineCoolingSimple
from ..starforge.CO_cooling import CO_cooling
from ..starforge.compton import compton_cooling
from ..starforge.nebular_cooling import nebular_cooling
from ..starforge.ionization_balance import ionization_processes
from ..starforge.h2_chemistry.collisional_detachment import Hminus_collisional_detachment
from ..starforge.dust_balance import TD_VARIABLE
from ..starforge.symbols import ISRF, cosmicray_attenuation_fac, X_H_SOLAR
from . import processes as mp

GIZMO_FAMILY = "starforge_menon"
SOLVE_VARS = ("u", "T", "H+", "He+", "He++", "H_2", "H_3+", "Td")
TIME_DEPENDENT = ("T", "H_2", "H+", "He+", "He++", "H_3+")

# starforge's cosmic-ray rate of H is 1.6e-17 sqrt(ISRF) times a column attenuation; KROME's is zeta/2, unattenuated
CR_SUBSTITUTION = {cosmicray_attenuation_fac: sp.S.One, sp.sqrt(ISRF): mp.zeta / (2 * 1.6e-17)}

PARAMETERS = [
    Parameter("chi_0", "Draine", 1.0, "FUV and ionizing field strength"),
    Parameter("Gamma_ion0", "s^-1", 2.19e-12, "unattenuated H photoionization rate of the driver's field"),
    Parameter("zeta_CR", "s^-1", 3e-17, "cosmic-ray ionization rate of H2 (KROME user_crate)"),
    Parameter("lambda_J", "cm", 3e19, "Jeans length sqrt(pi k T / (G rho mu m_p))"),
    Parameter("N_H", "cm^-2", 1e20, "shielding column n_Htot lambda_J / 2"),
    Parameter("rho", "g cm^-3", 2.3e-22, "gas density"),
    Parameter("x_D,tot", "", 3e-5, "deuterium nuclei per H nucleus"),
    Parameter("alpha_gd", "erg cm^3 s^-1 K^-3/2", 3.2e-34, "gas-dust energy exchange coefficient"),
    Parameter("Δx", "cm", 3e19, "length for the line-width estimate of the reused coolants"),
    Parameter("∇v", "s^-1", 1e-14, "velocity gradient, for the reused coolants"),
    Parameter("X", "", X_H_SOLAR, "H mass fraction"),
    Parameter("Z_d", "", 1.0, "dust-to-gas ratio relative to solar"),
    Parameter("f_d", "", 1.0, "dust survival fraction"),
    Parameter("z", "", 0.0, "redshift"),
    Parameter("G_0", "Habing", 1.69, "FUV field (derived from chi_0 and A_V)"),
    Parameter("f_neb", "", 1.0, "nebular cooling switch (derived: on)"),
]

SPECIES = SHARED_SPECIES + [
    Species("H-"), Species("H_2+"), Species("H_3+"),
    Species("HD", "trace", "HD per H nucleus; outside the EOS and the H budget"),
    Species("photon_assoc,H", "radiation", "photons of the radiative association of H-"),
]


def _cr_substituted(p):
    return Process(p.name, p.bibliography, {k: e.rhs.xreplace(CR_SUBSTITUTION) for k, e in p.network.items()})


def make_model():
    ion_procs, intermediates, electrons = ionization_processes()
    grain_recombination = [p for p in ion_procs if p.name == "Grain-assisted recombination of H+"]
    reused = (
        [*LineCoolingSimple("H"), *LineCoolingSimple("He+")]
        + [CollisionalIonization(s) for s in ("H", "He", "He+")]
        + [GasPhaseRecombination(i) for i in ("He+", "He++")]
        + [FreeFreeEmission(i) for i in ("H+", "He+", "He++")]
        + [compton_cooling, nebular_cooling, CO_cooling]
        + grain_recombination
        + [Hminus_collisional_detachment("H"), Hminus_collisional_detachment("e-")]
    )
    own = (mp.photochemistry() + mp.fine_structure_with_escape() + mp.cosmic_rays() + [mp.photoelectric_heating()]
           + mp.h2_network()
           + [mp.h2_cooling(), mp.hd_cooling(), mp.continuum_cooling()] + mp.dust_processes())
    model = Model(
        reused + own + [ThermalTerm(sp.Symbol("pdv_work"), name="PdV work")],
        solve_vars=SOLVE_VARS,
        time_dependent=TIME_DEPENDENT,
        steady_state=["H-"],
        fixed={"H_2+": sp.S.Zero, "HD": mp.hd_abundance(), **carbon_abundances()},
        derived={"C_2": sp.S.One, "C_3": sp.S.One, "G_0": mp.G_0_habing, "f_neb": sp.S.One},
        intermediates=[(s, e.xreplace(CR_SUBSTITUTION)) for s, e in intermediates],
        fixed_electrons=electrons.xreplace(CR_SUBSTITUTION),
        rules=[Rule("KROME cosmic-ray rate in the reused processes", _cr_substituted)],
        parameters=PARAMETERS,
        species=SPECIES,
        variables=[TD_VARIABLE],
    )
    return model
