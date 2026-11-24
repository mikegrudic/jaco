"""Gas-phase fractions of the heavy elements (the rest is in dust), shared by STARFORGE's cooling and its free-electron
budget so that both see the same gas.

TODO: one depletion law for every element. Jenkins (2009, ApJ 700, 1299; 2009ApJ...700.1299J) Eq. 10,
[X_gas/H] = B_X + A_X (F* - z_X), with the Table 4 coefficients gives Si and Fe here; C, O and Mg could follow it too
(F* = 1: C 0.61, O 0.58, Mg 0.054). F* is a modelling choice: F* = 1 is the cold-neutral-medium reference (the
zeta Oph cloud). Jenkins' F*-<n(H)> fit (Sec. 10.2) is in sight-line mean densities: do not feed it local cell
densities.
"""

import sympy as sp
from jaco.symbols import x_
from .symbols import x_solar

F_STAR = 1  # Jenkins' depletion strength: the cold neutral medium
# Jenkins 2009 Table 4: (A_X, B_X, z_X)
JENKINS_2009_TABLE4 = {
    "C": (-0.101, -0.193, 0.803),
    "O": (-0.225, -0.145, 0.598),
    "Mg": (-0.997, -0.800, 0.531),
    "Si": (-1.136, -0.570, 0.305),
    "Fe": (-1.285, -1.513, 0.437),
}
# gas-phase C/H at solar metallicity in diffuse clouds (Sofia et al. 2004, ApJ 605, 272), as GIZMO's electron budget
# takes it (return_electron_fraction_from_Cplus)
SOFIA_2004_C_GAS_PER_H = 1.6e-4


def jenkins_gas_fraction(element, F_star=F_STAR):
    """10^[X_gas/H] of Jenkins 2009 Eq. 10"""
    A, B, z = JENKINS_2009_TABLE4[element]
    return 10 ** (B + A * (F_star - z))


GAS_PHASE_FRACTION = {
    "C": SOFIA_2004_C_GAS_PER_H / x_solar("C"),
    "O": 1,  # undepleted
    "Mg": 1,  # undepleted, as GIZMO's heavy-ion cap assumes
    "Si": jenkins_gas_fraction("Si"),
    "Fe": jenkins_gas_fraction("Fe"),
}

# the model's total abundance per H nucleus: x_X,tot where network species share the element, else the species itself
TOTAL_ABUNDANCE = {"C": sp.Symbol("x_C,tot"), "O": sp.Symbol("x_O,tot"), "Mg": x_("Mg"), "Si": x_("Si"), "Fe": x_("Fe")}


def x_gas(element):
    """Gas-phase abundance of an element per H nucleus, all of its gas-phase species together"""
    return GAS_PHASE_FRACTION[element] * TOTAL_ABUNDANCE[element]
