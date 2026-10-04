"""GIZMO's legacy cooling module and its matter-radiation coupling under full-physics STARFORGE RT (M1 RADTRANSFER with
the ionizing, photoelectric, NUV, optical/NIR and infrared bands: SINGLE_STAR_FB_RAD), as one jaco model whose
unknowns include the energy of every band and the dust temperature.

starforge_legacy plus the processes of jaco.models.starforge.radiation, which hold the GIZMO references. GIZMO keeps
radiation transport only: the host hands the model each cell's post-transport band energies as the step's initial
values (x_photon_*_initial) and writes the solved ones back, and every exchange between the bands, the gas and the
dust is one of these processes, solved with the gas energy and the chemistry in one backward-Euler Newton system:

- photoionization of H by the ionizing band, each event taking one photon at c_tilde and donating its energy to the
  optical band (the kick's donation), heating the gas by eps_HI; the IR band's tail above 13.6 eV ionizes too;
- dust absorption of the photoelectric, NUV and optical bands and of the IR band, the dust's thermal emission into the
  IR band, gas-dust collisions, and gas absorption of the IR band (heating the gas at c_tilde/c, as GIZMO's kick does);
- the cooling radiation GIZMO returns to the bands (CoolingRate 1317-1348): metal and nebular lines, H and He+
  excitation, the recombination share f_recNUV and free-free emission above 1e5 K to the NUV band (to the IR band where
  T_rad > 1e4 K); molecular, fine-structure and atomic carbon lines, Compton and free-free emission below 1e5 K to the
  IR band;
- the dust temperature Td solves the dust's energy balance (zero heat capacity: absorption + gas-dust heating = emission),
  in place of rt_eqm_dust_temp and rt_ir_lambdadust;
- G_0, the H2 network's G_LW, the background temperature of the emission corrections T_bg, the IR self-absorption
  factor of the heating and cooling rates and the recombination share f_recNUV are expressions of the bands, the dust
  and the gas. The photoelectric heating and G_LW, rates over the step, see the photoelectric band averaged over it
  (as the dust absorbs it); the closures of the state at the step's end (the C+ fraction, the metals' free electrons)
  see the band there, so that the EOS does not depend on the step.

Deviations from GIZMO, each forced by doing the coupling in one implicit step:

- GIZMO absorbs each band in the kick with an exponential at frozen opacity and its dust temperature from that, then
  cools at fixed band energies; here absorption, emission and the dust temperature are terms of the same step, at the
  opacities of the solved state (Td, x_e, x_H+, H2). The dust's absorption of the photoelectric, NUV and optical bands
  keeps the kick's exponential at the start-of-step dust temperature (radiation.kick_absorption_factor,
  dust_band_rate), so those bands end the step at e^(-a Delta_t) of their initial energy, with the exponent capped at
  10 rather than GIZMO's 50, but what is added to them within the step (the NUV cooling return, the optical band's
  donation) is absorbed as if present from its start; the ionizing and IR bands are backward Euler;
- GIZMO's cooling return is limited by the gas energy change (de_u_touse); the processes give the bands exactly what
  the gas emits;
- the photon flux is left to the RT kick (whose relaxation with Rad_Kappa is the absorption's damping of it) but for the
  M1 limit; GIZMO's cooling return also scales the flux with the band's energy;
- the dust opacity table's composition switches are smoothed (jaco.models.starforge.dust_opacity);
- T_rad, the IR band's radiation temperature, is an input; the output T_rad_new is GIZMO's update of it from what the
  band kept and gained over the step, for the host to store. Solving for it would take the band's photon number as a
  second IR species (T_rad = its energy over its number times the mean photon energy's constant), with every IR
  process given a number row at the temperature it emits at.

Reproduced as GIZMO does them, though they do not conserve energy: the gas absorption of the IR band heats the gas at
c_tilde/c of the physical rate, at most its share of the band's energy per half-step kick, and is counted in full in
the dust balance as well, though the band gets only the dust's share back (radiation.gas_ir_absorption); the kick
puts the dust-absorbed energy of the photoelectric, NUV and optical bands into the IR band twice
(radiation.legacy_ir_donation_copy, a separate process: Model.without removes it); photoheating takes eps_HI per
photoionization while the band loses hnu_EUV to the optical band; photoelectric heating and H2
photodissociation do not take from the bands.
"""

import sympy as sp
from jaco.model import Rule
from jaco.declarations import Parameter, Species, Variable, Output
from jaco.process import Process
from ..starforge.starforge import GIZMO_FAMILY, SHARED_SPECIES  # noqa: F401  (GIZMO_FAMILY: the host-code family)
from ..starforge.symbols import z, nH_gizmo_cooling, n_Htot, T
from ..starforge.dust_opacity import T_dust, dust_survival
from ..starforge import radiation as rt
from ..starforge.radiation import EUV, FUV, NUV, ONIR, IR, BANDS
from ..starforge_legacy import make_model as make_legacy, scaled_densities, GIZMO_CLUMPING

T_CMB_Z0 = 2.73  # GIZMO's CMB temperature at z = 0, as jaco's starforge symbols write it: 2.73 (1 + z)
T_bg, f_IR_selfabs, f_recNUV = sp.Symbol("T_bg"), sp.Symbol("f_IR_selfabs"), sp.Symbol("f_recNUV")
DUST_TEMPERATURE_FLOOR = 2.73  # GIZMO's: max(MinGasTemp, T_CMB) under GALSF; the RT tests run MinGasTemp = 2.73
MAX_DUST_TEMP = 1.0e4  # GIZMO's MAX_DUST_TEMP

BAND_SPECIES = {
    EUV: Species(EUV, "radiation", "ionizing photons (13.6-500 eV) per H nucleus", floor=0.0),
    FUV: Species(FUV, "radiation", "photoelectric band (8-13.6 eV) energy per H nucleus [eV]", floor=0.0),
    NUV: Species(NUV, "radiation", "NUV band (3.444-8 eV) energy per H nucleus [eV]", floor=0.0),
    ONIR: Species(ONIR, "radiation", "optical/NIR band (0.4133-3.444 eV) energy per H nucleus [eV]", floor=0.0),
    IR: Species(IR, "radiation", "IR band (0.001-0.4133 eV) energy per H nucleus [eV]", floor=0.0),
}
RT_PARAMETERS = [
    Parameter("rsol", "", 1e-4, "reduced speed of light of the bands, c_tilde / c (RT_SPEEDOFLIGHT_REDUCTION)"),
    Parameter("sigma_HI", "cm^2", 0.0, "band-averaged H photoionization cross-section (rt_ion_sigma_HI)"),
    Parameter("eps_HI", "erg", 0.0, "heat per photoionization (rt_ion_G_HI)"),
    Parameter("hnu_EUV", "eV", 20.0, "mean energy of the ionizing photons (rt_nu_eff_eV)"),
]
IR_PARAMETERS = [
    Parameter("T_rad", "K", 10.0, "radiation temperature of the IR band (Radiation_Temperature)"),
    Parameter("Td_initial", "K", 20.0, "dust temperature at the start of the step (the dust-absorbed bands' opacity)"),
    Parameter("T_CMB", "K", T_CMB_Z0, "CMB temperature"),
    Parameter("rho", "g cm^-3", 2.3e-22, "gas density"),
    Parameter("Z_metals", "", 0.014, "metal mass fraction (Metallicity[0])"),
    Parameter("gamma_eos", "", 5.0 / 3.0, "the cell's adiabatic index at the start of the step"),
    Parameter("G_LW_bg", "Habing", 0.0, "UV background's contribution to the Lyman-Werner field"),
    Parameter("gamma_12_UVB", "", 0.0, "UV background's H photoionization rate / 1e-12 s^-1 (local_gammamultiplier gJH0)"),
    Parameter("eps_H0_UVB", "erg s^-1", 0.0, "UV background's photoheating per neutral H before shielding"),
]


def _material_rows_scaled(species_rows, factor):
    """Rule: every row of a species in species_rows divided by factor, the heat, energy and radiation rows kept"""
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


def _heat_scaled(factor):
    """Rule: the heat row multiplied by factor, the species and band rows kept"""
    def rule(p):
        rows = {k: e.rhs for k, e in p.network.items()}
        if rows.get("heat", 0) == 0:
            return p
        return Process(p.name, p.bibliography, {k: (v * factor if k == "heat" else v) for k, v in rows.items()})
    return rule


# GIZMO evaluates the rates at nHcgs = 0.76 rho/m_p (starforge_legacy's GIZMO_DENSITY); the bands are physical densities
GIZMO_DENSITY = Rule("GIZMO's nHcgs in the rates",
                     lambda p: p.transformed(scaled_densities(nH_gizmo_cooling / n_Htot, exclude=rt.PHOTON_DENSITIES)),
                     exempt={"GIZMO H2 network"})
GIZMO_ABUNDANCES = Rule("GIZMO integrates abundances per nHcgs",
                        _material_rows_scaled({s.name for s in SHARED_SPECIES if s.kind in ("material", "trace")},
                                              nH_gizmo_cooling / n_Htot),
                        exempt={"GIZMO H2 network"})
GIZMO_EMISSION_BACKGROUND = Rule("GIZMO's background temperature in the emission corrections",
                                 lambda p: p.transformed(_cmb_to_background),
                                 exempt={"Inverse Compton cooling (CMB)"})

# CoolingRate's routing of its cooling terms into the bands (1317-1348), by process: weight of the heat lost
ABOVE_1E5K = sp.Piecewise((1, T >= 1e5), (0, True))  # GIZMO's free-free routing (logT >= 5)
RECOMBINATION = [f"Gas-phase recombination of {i}" for i in ("H+", "He+", "He++")]
FREE_FREE = [f"Free-free emission from {i}" for i in ("H+", "He+", "He++")]
NUV_ROUTE = {"Metal line cooling": 1, "Nebular forbidden-line cooling": 1, "H-e- Line Cooling": 1,
             "He+-e- Line Cooling": 1, **{r: f_recNUV for r in RECOMBINATION}, **{f: ABOVE_1E5K for f in FREE_FREE}}
IR_ROUTE = {"GIZMO C+, [CI] and CO cooling": 1, "GIZMO H2 + HD cooling": 1, "Inverse Compton cooling (CMB)": 1,
            "Inverse Compton cooling (RT bands)": 1, **{f: 1 - ABOVE_1E5K for f in FREE_FREE}}
NUV_TO_IR = sp.Piecewise((1, sp.Symbol("T_rad") > 1e4), (0, True))  # where the IR band is hotter than 1e4 K


def _route_cooling(p):
    lost = -p.heat * rt.rsol / rt.EV  # radiated energy, at c_tilde, in the bands' eV
    w_nuv, w_ir = NUV_ROUTE.get(p.name, 0), IR_ROUTE.get(p.name, 0)
    rows = {}
    if w_nuv != 0:
        rows[NUV] = lost * w_nuv * (1 - NUV_TO_IR)
    if w_ir != 0 or w_nuv != 0:
        rows[IR] = lost * (w_ir + w_nuv * NUV_TO_IR)
    return p.with_rows(rows)


GIZMO_COOLING_RADIATION = Rule("GIZMO's cooling radiation into the NUV and IR bands", _route_cooling,
                               only=set(NUV_ROUTE) | set(IR_ROUTE))
# The photoelectric heating, a rate over the step, sees the photoelectric band over the step, as the H2 network's G_LW
# does; the closures of the state at the step's end (the C+ fraction, the metals' free electrons, grain charging in
# them) see the band there, so that u(T, x) does not depend on the step
G_0_step = sp.Symbol("G_0_step")
PHOTOELECTRIC_OVER_STEP = Rule("the photoelectric heating over the step",
                               lambda p: p.transformed(lambda e: e.xreplace({sp.Symbol("G_0"): G_0_step})),
                               only={"Photoelectric Heating"})
GIZMO_IR_SELF_ABSORPTION = Rule("GIZMO's IR self-absorption of the heating and cooling rates", _heat_scaled(f_IR_selfabs),
                                exempt={"Gas-dust collisions", "PdV work", "Gas absorption of photon_IR"})

RADIATION_PROCESSES = [
    rt.photoionization(donation=ONIR),
    rt.ir_tail_photoionization(),
    *[rt.dust_band_absorption(b) for b in (FUV, NUV, ONIR)],
    rt.legacy_ir_donation_copy(),
    rt.dust_ir_absorption(),
    rt.gas_ir_absorption(),
    rt.dust_ir_emission(),
    rt.compton_off_bands(BANDS),
]
IR_ABSORBERS = ["Dust absorption of photon_IR", "Gas absorption of photon_IR"]
DONATION_COPY = "GIZMO's second copy of the donated dust absorption in photon_IR"
DERIVED = {
    "f_d": dust_survival(T_dust),
    "G_0": rt.G0_of_band(),
    "G_0_step": rt.G0_of_band(over_step=True),
    "G_LW": rt.G_LW_of_bands(),
    "T_bg": rt.background_temperature(),
    "f_IR_selfabs": rt.ir_self_absorption(),
    "f_recNUV": rt.recombination_return_fraction(),
}


def _outputs():
    # GIZMO's direct donation (here the dust's re-emission of what it absorbs) is in the IR band before the kick's
    # update counts it at T_rad; its second copy is the update's dust emission: the copy takes the first's weight
    T_new, sums = rt.ir_radiation_temperature(IR_ABSORBERS, ["Dust emission into photon_IR"],
                                              sorted(set(NUV_ROUTE) | set(IR_ROUTE)), prior_sources=[DONATION_COPY])
    return [
        Output("T_rad_new", T_new, units="K", sums=sums,
               doc="IR radiation temperature after the step, GIZMO's photon-number weighting"),
        Output("photoionization_rate", heat_of={(p, EUV): -1 / rt.rsol for p in rt.EUV_SINKS},
               units="cm^-3 s^-1", doc="photoionizations by the ionizing band per unit volume (each takes one photon)"),
    ]


def make_model():
    """The STARFORGE_LEGACY_RT model"""
    base = make_legacy()
    processes = [rt.gas_dust_collisions() if i == "Gas-dust collisions" else p for i, p in base.processes.items()]
    dropped = set(DERIVED) | {"Td"}
    return base.evolve(
        processes + RADIATION_PROCESSES,
        solve_vars=base.solve_vars + (EUV, FUV, NUV, ONIR, IR, "Td"),
        time_dependent=("T", "H+", "H_2", EUV, FUV, NUV, ONIR, IR),
        rules=[GIZMO_CLUMPING, GIZMO_DENSITY, GIZMO_ABUNDANCES, GIZMO_EMISSION_BACKGROUND, GIZMO_IR_SELF_ABSORPTION,
               GIZMO_COOLING_RADIATION, PHOTOELECTRIC_OVER_STEP],
        derived={**base.derived, **DERIVED},
        intermediates=rt.kick_absorption_intermediates() + list(base.intermediates),
        parameters=[p for p in base.parameters if p.name not in dropped] + RT_PARAMETERS + IR_PARAMETERS,
        species=[*base.species, *BAND_SPECIES.values()],
        variables=[Variable("Td", "dust heat", floor=DUST_TEMPERATURE_FLOOR, ceiling=MAX_DUST_TEMP, units="K",
                            doc="dust temperature: the steady state of the dust's energy balance")],
        outputs=_outputs(),
    )
