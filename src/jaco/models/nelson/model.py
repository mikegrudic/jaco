"""
The Nelson & Langer (1999) model.

Chemical reactions and rates are in chemistry.py.
Thermal terms are in thermal.py.
"""

import sympy as sp
from jaco.declarations import Parameter, Species
from jaco.model import Model
from jaco.processes import ThermalTerm
from .chemistry import reactions
from .thermal import thermal_terms

# Parameters from NL99 appendix. NL99 totals are per H_2, but here they are per H nucleus (x2 fewer).
# These are used to reproduce results in NL99 during testing.
PAPER_PARAMETERS = {
    "y"        : 0.14,
    "x_C,tot"  : 5e-5,
    "x_O,tot"  : 1e-4,
    "x_Mg,tot" : 5e-8
}

SOLVED   = ("H3+", "He+", "CH", "OH", "CO", "C+", "HCO+", "Mg+")
MATERIAL = ("H_2", "H3+", "He", "He+", "C", "C+", "CH", "O", "OH", "CO", "HCO+", "Mg", "Mg+", "e-")
H2_PER_H = sp.Rational(1, 2)

DOCS = {
    "CH"  : "CHx of NL99: CH + CH2",
    "OH"  : "OHx of NL99: OH + H2O + O2",
    "Mg"  : "M of NL99: Mg, Fe, Ca, Na",
    "Mg+" : "M+ of NL99",
    "H_2" : "molecular hydrogen, all of the H"
}

PARAMETERS = (
    Parameter(   "G_0", "Habing", 1.0, "FUV field relative to the Habing (1968) field"),
    Parameter(   "A_V",    "mag", 0.0, "visual extinction of the photoreactions and the dust temperature"),
    Parameter("tau_uv",       "", 0.0, "UV optical depth of the photoelectric and surface terms"),
    Parameter(  "N_CO",  "cm^-2", 1.0, "CO column density to the surface, for the CO self-shielding (1 = unshielded)"),
    Parameter(  "N_H2",  "cm^-2", 1.0, "H2 column density to the surface, for the shielding of CO (1 = unshielded)"),
    Parameter( "dv_dr",    "s^-1", 3.2408e-14, "line-width (velocity) gradient of the C I escape probabilities, NL99: 1 km/s/pc"),
)


def make_model():
    """
    The NL99 network with the thermal terms of NL97 and NL99, and the gas internal energy and temperature solved.

    The cloud is taken to be molecular: H_2 is fixed at 1/2 per H nucleus.

    Pass PAPER_PARAMETERS for the elemental abundances of the paper.
    """
    pdv_work = ThermalTerm(sp.Symbol("pdv_work"), name="PdV work")
    return Model(
        processes      = [*reactions(), *thermal_terms(), pdv_work],
        solve_vars     = ["u", "T", *SOLVED],
        time_dependent = ["T", *SOLVED],
        fixed          = {"H_2": H2_PER_H},
        parameters     = PARAMETERS,
        species        = [Species(s, doc=DOCS.get(s, "")) for s in MATERIAL]
    )
