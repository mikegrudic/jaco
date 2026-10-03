"""STARFORGE ISM thermochemistry: one process library, two models.

``make_model(STARFORGE)`` is the physically preferred model; ``make_model(STARFORGE_LEGACY)`` (the
``starforge_legacy`` package) reproduces GIZMO's legacy cooling module term for term. The switches that separate
them, with their GIZMO references, are in switches.py. Both solve u, T, x_H+, x_He+, x_He++ and x_H2 (T and H2 in
time, the ions in steady state) and share the KWH H/He ionization balance and its cooling, gas-dust coupling,
cosmic-ray and photoelectric heating, CMB Compton cooling, tabulated metal lines (switched by f_metal), Kim+23 nebular
lines (f_neb) and the PdV term.

STARFORGE in addition: the H2/H- reaction network with the heat of H2 formation and dissociation; cosmic-ray
ionization of H with radiative, grain-assisted and charge-transfer sinks, and the C+, Mg+ and molecular-ion balances
for the other free electrons (ionization_balance.py); C+ (on the Tielens C+/CO interpolation), [CI] 609 um on neutral
carbon and Whitworth & Jaffa CO cooling; H2/HD cooling with number-weighted colliders; clumping on every two-body rate.
"""

import sympy as sp
from jaco.processes import CollisionalIonization, GasPhaseRecombination, FreeFreeEmission, ThermalProcess
from jaco.process import Process
from jaco.equation import Equation
from jaco.symbols import x_
from .switches import Switches, STARFORGE, STARFORGE_LEGACY
from .line_cooling import LineCoolingSimple, CI_cooling
from .h2_chemistry import H2_chemistry
from .h2_chemistry.gizmo_network import GizmoH2Network
from .H2_cooling import H2_cooling
from .CO_cooling import CO_cooling
from .gizmo_lowtemp import gizmo_carbon_cooling, gizmo_H2_cooling
from .gas_dust_collisions import gas_dust_collisions
from .cosmic_ray_ionization import cosmic_ray_heating
from .photoelectric_heating import photoelectric_heating
from .metal_line_cooling import metal_line_cooling
from .metal_electrons import metal_electrons
from .ionization_balance import ionization_processes
from .nebular_cooling import nebular_cooling
from .compton import compton_cooling
from .symbols import T, n_Htot, G_0, clumping_factor, nH_gizmo_cooling

# Solve variables in index order (GIZMO's jaco.cc assumes u, T first, abundances after) and the
# subset that gets a backward-Euler term; the ions are solved in steady state.
SOLVE_VARS = ["u", "T", "H+", "He+", "He++", "H_2"]
TIME_DEPENDENT = ["T", "H_2"]
GIZMO_FAMILY = "starforge"  # host-code family macro JACO_FAMILY_STARFORGE, shared by both models
DEUTERIUM_PER_H = 2.527e-5  # Cooke, Pettini & Steidel 2018


def _transformed(process, rule):
    """Copy of process with rule (a function of an expression) applied to every rate and heat expression"""
    out = Process(name=process.name, bibliography=process.bibliography)
    for key, eq in process.network.items():
        if key == "heat":
            out.heat = rule(eq.rhs)
        else:
            out.network[key] = Equation(eq.lhs, rule(eq.rhs))
    return out


def scaled_densities(process, factor):
    """process with every number density in its rates (n_X, n_Htot) multiplied by factor"""
    def rule(e):
        e = sp.sympify(e)
        return e.xreplace({s: factor * s for s in e.free_symbols if str(s).startswith("n_")})
    return _transformed(process, rule)


def unclumped(process):
    """process with the clumping factors C_2 and C_3 set to 1"""
    return _transformed(process, lambda e: sp.sympify(e).xreplace({sp.Symbol("C_2"): 1, sp.Symbol("C_3"): 1}))


def make_model(switches: Switches = STARFORGE):
    """Build the model selected by switches (STARFORGE or STARFORGE_LEGACY)."""
    s = switches
    if s.h2_network == "gizmo" and s.h2_chemical_heat:
        raise ValueError("GIZMO's H2 network carries no chemical heat")
    C2 = sp.Symbol("C_2")

    processes = [
        *[LineCoolingSimple(sp_) for sp_ in ("H", "He+")],
        *[CollisionalIonization(sp_, clumping=C2) for sp_ in ("H", "He", "He+")],
        *[GasPhaseRecombination(i) for i in ("H+", "He+", "He++")],
        *[FreeFreeEmission(i) for i in ("H+", "He+", "He++")],
        gas_dust_collisions,
        cosmic_ray_heating,
        photoelectric_heating,
        compton_cooling,
        metal_line_cooling(z=0.0),
        nebular_cooling,
    ]
    processes += [H2_cooling] if s.h2_cooling == "network" else [gizmo_H2_cooling]
    processes += [LineCoolingSimple("C+"), CO_cooling, CI_cooling] if s.carbon_cooling == "network" else [gizmo_carbon_cooling]
    if s.electrons == "solved":
        ion_processes, intermediates, fixed_electrons = ionization_processes()
        processes += ion_processes
    else:
        n_H = nH_gizmo_cooling if s.rate_density == "gizmo" else n_Htot
        intermediates, fixed_electrons = metal_electrons(n_H)

    if s.rate_density == "gizmo":
        processes = [scaled_densities(p, nH_gizmo_cooling / n_Htot) for p in processes]
    if s.clumping == "h2_chemistry":
        processes = [unclumped(p) for p in processes]
    # the H2 network takes its densities and clumping as written (GIZMO's own, for the gizmo network)
    h2 = H2_chemistry(s.h2_chemical_heat) if s.h2_network == "per_molecule" else GizmoH2Network()
    model = sum(processes + [h2, ThermalProcess(sp.Symbol("pdv_work"), name="PdV work")])

    fixed = {}
    if s.h2_network == "per_molecule":
        fixed["H_2+"] = sp.S.Zero  # negligible in ISM conditions
        # H-: equilibrium of its rate equation, linear in n_H-: source + sink_coeff * n_H- = 0
        nHm = sp.Symbol("n_H-")
        source, sink_coeff = sp.S.Zero, sp.S.Zero
        for term in sp.Add.make_args(model.network["H-"].rhs):
            if nHm not in term.free_symbols:
                source += term
            else:
                sink_coeff += sp.simplify(term / nHm)
        fixed["H-"] = sp.Max(0, -source / sink_coeff / n_Htot)
    if s.h2_cooling == "network":
        fixed["HD"] = 2 * DEUTERIUM_PER_H * x_("H_2")  # all D in HD in proportion to the molecular fraction 2 x_H2
    if s.carbon_cooling == "network":
        # C+ and CO from GIZMO's cooling-curve interpolation, fco/(1 - fco) ~ (n / 340 G0)^2 / sqrt(T) (Tielens)
        x_C_tot = sp.Symbol("x_C,tot")
        f_Cp = 1 / (1 + (n_Htot / (340 * sp.Max(sp.Rational(1, 10), G_0)))**2 / sp.sqrt(sp.Max(T, 10)))
        fixed["C+"] = x_C_tot * f_Cp
        fixed["CO"] = sp.Max(1e-30, x_C_tot * (1 - f_Cp) * x_("H_2"))
    else:
        fixed.update({"C+": sp.S.Zero, "CO": sp.S.Zero})  # carbon counted as atoms in the EOS and metal-line tables
    model.network.fixed_species = fixed
    model.network.equilibrium_overrides = {k: v for k, v in fixed.items() if k in ("H_2+", "HD")}

    # Expressions of the solve variables, substituted before code generation so their derivatives enter the Jacobian
    model.network.derived_params = {"C_2": clumping_factor, "C_3": clumping_factor**3}  # <n^3>/<n>^3 = C_2^3, lognormal
    # Free electrons beyond the solved H and He ions
    model.network.intermediates, model.network.fixed_electrons = intermediates, fixed_electrons
    return model
