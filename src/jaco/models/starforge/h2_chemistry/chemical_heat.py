"""Heat of H2 formation and of collisional dissociation.

Formation releases the 4.48 eV binding energy partly as kinetic energy and partly as rovibrational excitation, which
heats the gas only where collisions de-excite it faster than it radiates, i.e. above the critical density of H2:
Hollenbach & McKee 1979 for formation on grains (0.2 eV kinetic + 4.2 eV internal) and via H- (3.53 eV), Omukai 2000
for three-body formation (4.48 eV), as compiled by Nickerson, Teyssier & Rosdahl 2018 (Eqs. 46-47). Collisional
dissociation takes the full binding energy from the gas.
"""

import sympy as sp
from ..symbols import T, n_Htot, x_, EV_CGS, H2_BINDING_ENERGY

BIBLIOGRAPHY = ["1979ApJS...41..555H", "2000ApJ...534..809O", "2018MNRAS.479.3206N"]


def deexcitation_fraction():
    """1 / (1 + n_cr / n_H) with HM79's n_cr = 1e6 T^-1/2 / (1.6 x_H exp(-(400/T)^2) + 1.4 x_H2 exp(-12000/(T + 1200)))
    cm^-3, written without the division by the collider abundances"""
    w = n_Htot * (1.6 * x_("H") * sp.exp(-((400 / T) ** 2)) + 1.4 * x_("H_2") * sp.exp(-12000 / (T + 1200)))
    return w / (w + 1e6 / sp.sqrt(T))


def formation_heat(channel, on=True):
    """Heat given to the gas per H2 formed through channel ("grain", "H-" or "3-body"), in erg"""
    if not on:
        return 0
    f = deexcitation_fraction()
    return {"grain": 0.2 + 4.2 * f, "H-": 3.53 * f, "3-body": 4.48 * f}[channel] * EV_CGS


def dissociation_heat(on=True):
    """Heat given to the gas per collisional dissociation of H2, in erg"""
    return -H2_BINDING_ENERGY if on else 0
