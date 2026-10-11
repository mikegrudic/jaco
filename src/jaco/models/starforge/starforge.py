"""STARFORGE ISM thermochemistry: one process library, two models.

``make_model()`` here is the physically preferred model; ``jaco.models.starforge_legacy`` reproduces GIZMO's legacy
cooling module term for term, and its docstring lists what separates the two, with GIZMO references. Both solve u, T,
x_H+, x_He+, x_He++ and x_H2 (T and H2 in time, the ions in steady state) and share the KWH H/He ionization balance and
its cooling, gas-dust coupling, cosmic-ray and photoelectric heating, CMB Compton cooling, tabulated metal lines
(switched by f_metal), Kim+23 nebular lines (f_neb) and the PdV term (shared_processes()).

STARFORGE in addition: the H2/H- reaction network with the heat of H2 formation and dissociation, H- in steady state,
and the energy formation radiates as a cooling term under the optically-thick limit;
cosmic-ray ionization of H with radiative, grain-assisted and charge-transfer sinks, and the C+, Mg+ and molecular-ion
balances for the other free electrons (ionization_balance.py); the fine-structure lines of C+, C, O, Si+ and Fe+ from
levels in statistical equilibrium (fine_structure.py; C+ and CO on the Tielens interpolation) and Whitworth & Jaffa CO
cooling; H2/HD cooling with number-weighted colliders; clumping on every two-body rate, estimated from the trace-free
velocity gradient and tapered off above ~5000 K. The cooling and the electron budget see the same gas-phase
(undepleted) abundances (gas_phase.py). Its cooling has GIZMO's optically-thick limit (optically_thick.py), which
GIZMO's CoolingRate applies and the starforge_legacy transcription lacks.

Both solve the dust temperature with the gas from GIZMO's non-RT dust balance (dust_balance.py).
"""

import sympy as sp
from jaco.model import Model
from jaco.declarations import Parameter, Species
from jaco.processes import CollisionalIonization, GasPhaseRecombination, FreeFreeEmission, ThermalTerm
from jaco.symbols import x_
from .line_cooling import LineCoolingSimple
from .fine_structure import fine_structure_cooling
from .h2_chemistry import h2_chemistry_processes, h2_formation_radiation
from .H2_cooling import H2_cooling
from .CO_cooling import CO_cooling
from .gas_dust_collisions import gas_dust_collisions
from .cosmic_ray_ionization import cosmic_ray_heating
from .photoelectric_heating import photoelectric_heating
from .metal_line_cooling import metal_line_cooling
from .ionization_balance import ionization_processes
from .nebular_cooling import nebular_cooling
from .compton import compton_cooling
from .gas_phase import x_gas
from .optically_thick import with_optically_thick_cap
from .dust_balance import with_dust_balance
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
    """C+ and CO of the gas-phase carbon from GIZMO's cooling-curve interpolation, fco/(1 - fco) ~ (n / 340 G0)^2 /
    sqrt(T) (Tielens); the carbon outside C+ is in CO in proportion to the molecular fraction 2 x_H2"""
    f_Cp = 1 / (1 + (n_Htot / (340 * sp.Max(sp.Rational(1, 10), G_0)))**2 / sp.sqrt(sp.Max(T, 10)))
    return {"C+": x_gas("C") * f_Cp, "CO": sp.Max(1e-30, x_gas("C") * (1 - f_Cp) * 2 * x_("H_2"))}


def make_model(chemical_heat=True, optically_thick=True, three_body_detailed_balance=False):
    """The STARFORGE model. chemical_heat: the gas takes the heat of H2 formation and pays that of dissociation;
    optically_thick: GIZMO's optically-thick limit on the cooling (optically_thick.py), which a model without a
    transported IR band needs in any gas dense enough to trap its cooling radiation; three_body_detailed_balance:
    three-body H2 formation by detailed balance with the dissociation (h2_chemistry/detailed_balance.py)"""
    ion_processes, intermediates, fixed_electrons = ionization_processes()
    h2_chemistry = h2_chemistry_processes(chemical_heat=chemical_heat,
                                          three_body_detailed_balance=three_body_detailed_balance)
    model = Model(
        shared_processes()
        + [H2_cooling, *fine_structure_cooling(), CO_cooling, *ion_processes]
        + h2_chemistry
        + h2_formation_radiation(h2_chemistry)
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
    model = with_dust_balance(model)
    if optically_thick:  # chemical heat and the gas-dust exchange are not radiation; what formation radiates is
        model = with_optically_thick_cap(model, exempt=[p.name for p in h2_chemistry] + ["Gas-dust collisions"])
    return model
