"""STARFORGE ISM thermochemistry: one process library, two models.

``make_model()`` here is the physically preferred model; ``jaco.models.starforge_legacy`` reproduces GIZMO's legacy
cooling module term for term, and its docstring lists what separates the two, with GIZMO references. Both solve u, T,
x_H+, x_He+, x_He++ and x_H2 (T and H2 in time, the ions in steady state) and share the KWH H/He ionization balance and
its cooling, gas-dust coupling, cosmic-ray and photoelectric heating, CMB Compton cooling, tabulated metal lines
(switched by f_metal), Kim+23 nebular lines (f_neb) and the PdV term (shared_processes()).

STARFORGE in addition: the H2/H- reaction network with the heat of H2 formation and dissociation, H- in steady state;
cosmic-ray ionization of H with radiative, grain-assisted and charge-transfer sinks, and the C+, Mg+ and molecular-ion
balances for the other free electrons (ionization_balance.py); C+ (on the Tielens C+/CO interpolation), [CI] 609 um on
neutral carbon and Whitworth & Jaffa CO cooling; H2/HD cooling with number-weighted colliders; clumping on every
two-body rate, estimated from the trace-free velocity gradient and tapered off above ~5000 K.
"""

import sympy as sp
from jaco.model import Model
from jaco.declarations import Parameter, Species
from jaco.processes import CollisionalIonization, GasPhaseRecombination, FreeFreeEmission, ThermalTerm
from jaco.symbols import x_
from .line_cooling import LineCoolingSimple, CI_cooling
from .h2_chemistry import h2_chemistry_processes
from .H2_cooling import H2_cooling
from .CO_cooling import CO_cooling
from .gas_dust_collisions import gas_dust_collisions
from .cosmic_ray_ionization import cosmic_ray_heating
from .photoelectric_heating import photoelectric_heating
from .metal_line_cooling import metal_line_cooling
from .ionization_balance import ionization_processes
from .nebular_cooling import nebular_cooling
from .compton import compton_cooling
from .symbols import T, n_Htot, G_0, clumping_factor, clumping_factor_starforge, PARAMETERS

# Solve variables in index order (GIZMO's jaco.cc assumes u, T first, abundances after) and the subset that gets a
# backward-Euler term; the ions are solved in steady state
SOLVE_VARS = ("u", "T", "H+", "He+", "He++", "H_2")
TIME_DEPENDENT = ("T", "H_2")
GIZMO_FAMILY = "starforge"  # host-code family macro JACO_FAMILY_STARFORGE, shared by both models
DEUTERIUM_PER_H = 2.527e-5  # Cooke, Pettini & Steidel 2018


def clumping(C_2):
    """C_2 and C_3 = C_2^3 (<n^3>/<n>^3 of a lognormal), as expressions of the solve variables, substituted before code
    generation so their derivatives enter the Jacobian"""
    return {"C_2": C_2, "C_3": C_2**3}


DERIVED = clumping(clumping_factor)  # GIZMO's estimator, on the full velocity gradient (starforge_legacy)
STARFORGE_DERIVED = clumping(clumping_factor_starforge)  # trace-free gradient, cold gas only
STARFORGE_PARAMETERS = PARAMETERS + [
    Parameter("∇v_tf", "s^-1", 1e-14, "Frobenius norm of the trace-free velocity gradient, for the clumping factor")]

pdv_work = ThermalTerm(sp.Symbol("pdv_work"), name="PdV work")

# Species of both models. C+ and CO are fixed; the other metals enter the EOS and the metal-line tables with their
# total abundances
SHARED_SPECIES = [Species(s) for s in ("H", "H+", "H_2", "He", "He+", "He++", "e-", "C", "C+", "CO", "O")]
SHARED_SPECIES += [Species(el, doc=f"{el} nuclei per H nucleus, all of them") for el in ("N", "Ne", "Mg", "Si", "S", "Ca", "Fe")]
SHARED_SPECIES.append(Species("dust heat", "energy", "energy the gas gives the dust"))
STARFORGE_SPECIES = SHARED_SPECIES + [
    Species("H-"), Species("H_2+"),
    Species("HD", "trace", "HD per H nucleus; outside the EOS and the H budget"),
    Species("photon_assoc,H", "radiation", "photons of the radiative association of H-"),
]


def shared_processes():
    """The terms both models carry, in the order they are summed"""
    return [
        *LineCoolingSimple("H"),
        *LineCoolingSimple("He+"),
        *[CollisionalIonization(s) for s in ("H", "He", "He+")],
        *[GasPhaseRecombination(i) for i in ("H+", "He+", "He++")],
        *[FreeFreeEmission(i) for i in ("H+", "He+", "He++")],
        gas_dust_collisions,
        cosmic_ray_heating,
        photoelectric_heating,
        compton_cooling,
        metal_line_cooling(z=0.0),
        nebular_cooling,
    ]


def carbon_abundances():
    """C+ and CO from GIZMO's cooling-curve interpolation, fco/(1 - fco) ~ (n / 340 G0)^2 / sqrt(T) (Tielens); the
    carbon outside C+ is in CO in proportion to the molecular fraction 2 x_H2"""
    x_C_tot = sp.Symbol("x_C,tot")
    f_Cp = 1 / (1 + (n_Htot / (340 * sp.Max(sp.Rational(1, 10), G_0)))**2 / sp.sqrt(sp.Max(T, 10)))
    return {"C+": x_C_tot * f_Cp, "CO": sp.Max(1e-30, x_C_tot * (1 - f_Cp) * 2 * x_("H_2"))}


def make_model():
    """The STARFORGE model"""
    ion_processes, intermediates, fixed_electrons = ionization_processes()
    return Model(
        shared_processes()
        + [H2_cooling, *LineCoolingSimple("C+"), CO_cooling, CI_cooling, *ion_processes]
        + h2_chemistry_processes(chemical_heat=True)
        + [pdv_work],
        solve_vars=SOLVE_VARS,
        time_dependent=TIME_DEPENDENT,
        steady_state=["H-"],
        fixed={
            "H_2+": sp.S.Zero,  # negligible in ISM conditions
            "HD": 2 * DEUTERIUM_PER_H * x_("H_2"),  # all D in HD in proportion to the molecular fraction 2 x_H2
            **carbon_abundances(),
        },
        derived=STARFORGE_DERIVED,
        intermediates=intermediates,
        fixed_electrons=fixed_electrons,  # free electrons beyond the solved H and He ions
        parameters=STARFORGE_PARAMETERS,
        species=STARFORGE_SPECIES,
    )
