"""The 23 reactions of Table 2 of Nelson & Langer (1999, NL99; ApJ 524, 923), with their rate coefficients (cm^3 s^-1,
or s^-1 for the cosmic-ray and photoreactions), T in K.

Species names: jaco reads the element content of a species from its name, which the composites CHx, OHx and M do not
have. They are written as CH, OH and Mg (Mg+): CH and OH stand for the radical composites (CH + CH2; OH + H2O + O2),
Mg for the low-ionization metals (Mg, Fe, Ca, Na). With these names the conservation conditions of the paper come out
of jaco's own reductions: charge neutrality gives e-, and C, O, He and M (here Mg) follow from their totals.

Hydrogen: as in the paper the cloud is molecular (H_2 fixed in the model, per H nucleus 1/2), so H and H+ products of
reactions, which the paper does not follow, are left out. The photoreactions use G_0 (Habing units) e^(-gamma A_V); the
CO photodissociation also carries S(N(CO), N(H_2)) from van Dishoeck & Black (1988) Table 5 (shielding.py).
The host supplies both column densities. With both columns at 1 cm^-2, the shielding factor is unity.
"""

import sympy as sp
from jaco.processes import Reaction
from jaco.symbols import T
from .shielding import co_shielding

NL99  = "1999ApJ...524..923N"
# The rates are transcribed from NL99 Table 2, which credits the UMIST catalog and cites Bergin et al. (1995).
UMIST = "1995ApJ...441..222B"
VDB88 = "1988ApJ...334..771V"

G_0 = sp.Symbol("G_0")
A_V = sp.Symbol("A_V")
S_CO = co_shielding(sp.Symbol("N_CO"), sp.Symbol("N_H2"))


def photoreaction(equation, alpha, gamma, shielding=1, extra_bibliography=()):
    return Reaction(
        equation,
        alpha * G_0 * sp.exp(-gamma * A_V) * shielding,
        name=f"Photo: {equation}",
        bibliography=[NL99, UMIST, *extra_bibliography]
    )


def reactions():
    """The reactions of NL99 Table 2, in the table's order"""
    R = []

    def add(equation, k, name=None):
        R.append(Reaction(equation, k, name=name or equation, bibliography=[NL99, UMIST]))

    # cosmic-ray ionization (per H_2 and per He; the H product is not followed)
    add("H_2 -> H3+ + e-", 1.2e-17, "cr + H_2 -> H3+ + e- + H")
    add("He -> He+ + e-",  6.8e-18, "cr + He -> He+ + e-")

    # ion-molecule
    add("H3+ + C -> CH + H_2",       2e-9)
    add("H3+ + O -> OH + H_2",       8e-10)
    add("H3+ + CO -> HCO+ + H_2",  1.7e-9)
    # NL99: He+ + H_2 -> He + H + H+; the H+ is not followed. Without it the reaction does not conserve charge, which
    # only matters for e-, and that follows from charge neutrality
    add("He+ + H_2 -> He",           7e-15, "He+ + H_2 -> He + H + H+")
    add("He+ + CO -> C+ + O + He", 1.6e-9)
    R.append(Reaction("C+ + H_2 -> CH", 4e-16, name="C+ + H_2 -> CH + H", bibliography=[NL99, UMIST]))
    add("C+ + OH -> HCO+",           1e-9)

    # neutral-neutral
    add("O + CH -> CO",   2e-10,          "O + CH -> CO + H")
    add("C + OH -> CO", 5.8e-12 * T**0.5, "C + OH -> CO + H")

    # electron recombination
    add("He+ + e- -> He",    9e-11 * T**-0.64)
    add("H3+ + e- -> H_2", 1.9e-6  * T**-0.54, "H3+ + e- -> H + H_2")
    add("C+ + e- -> C",    1.4e-10 * T**-0.61)
    add("HCO+ + e- -> CO", 3.3e-5  * T**-1.0, "HCO+ + e- -> CO + H")
    add("Mg+ + e- -> Mg",  3.8e-10 * T**-0.65)

    # charge transfer
    add("H3+ + Mg -> Mg+ + e- + H_2", 2e-9)

    # photoreactions
    R += [
        photoreaction("C -> C+ + e-",     3e-10, 3.0),
        photoreaction("CH -> C",          1e-9,  1.5),
        photoreaction("CO -> C + O",      1e-10, 3.0, S_CO, extra_bibliography=[VDB88]),
        photoreaction("OH -> O",          5e-10, 1.7),
        photoreaction("Mg -> Mg+ + e-",   2e-10, 1.9),
        photoreaction("HCO+ -> CO",     1.5e-10, 2.5),
    ]
    return R
