"""Heat of H2 formation and of collisional dissociation.

Formation releases the 4.48 eV binding energy partly as kinetic energy and partly as rovibrational excitation, which
heats the gas only where collisions de-excite it faster than it radiates, i.e. above the critical density of H2:
Hollenbach & McKee 1979 for formation on grains (0.2 eV kinetic + 4.2 eV internal) and via H- (3.53 eV), Omukai 2000
for three-body formation (4.48 eV), as compiled by Nickerson, Teyssier & Rosdahl 2018 (Eqs. 46-47). Collisional
dissociation takes the full binding energy from the gas.

The binding energy formation does not give the gas (formation_radiated: the H2 lines, on the H- route also the 0.754 eV
attachment photon, on grains also what the grain keeps) is radiation. A formation reaction gives the gas the whole 4.48
eV, and formation_radiation() takes the remainder back as a separate cooling term, which a model's optically-thick limit
caps like the rest of its cooling: optically thin gas gets formation_heat per molecule, and gas thick to its own cooling
keeps the remainder, so that a dissociate-and-reform cycle conserves its energy there.
"""

import sympy as sp
from jaco.processes import ThermalTerm
from ..symbols import T, n_Htot, x_, EV_CGS, H2_BINDING_ENERGY

BIBLIOGRAPHY = ["1979ApJS...41..555H", "2000ApJ...534..809O", "2018MNRAS.479.3206N"]
RADIATION_OF = "Radiation of "  # formation_radiation's name: this and the reaction's


def deexcitation_fraction():
    """1 / (1 + n_cr / n_H) with HM79's n_cr = 1e6 T^-1/2 / (1.6 x_H exp(-(400/T)^2) + 1.4 x_H2 exp(-12000/(T + 1200)))
    cm^-3, written without the division by the collider abundances"""
    w = n_Htot * (1.6 * x_("H") * sp.exp(-((400 / T) ** 2)) + 1.4 * x_("H_2") * sp.exp(-12000 / (T + 1200)))
    return w / (w + 1e6 / sp.sqrt(T))


def formation_heat(channel, on=True):
    """Heat given to the gas per H2 formed through channel ("grain", "H-" or "3-body") where its radiation escapes, in
    erg: 0.2 + 4.2 f eV on grains and 3.53 f eV by associative detachment (HM79), 4.48 f eV for three-body formation, f
    the de-excitation fraction"""
    if not on:
        return 0
    f = deexcitation_fraction()
    return {"grain": 0.2 + 4.2 * f, "H-": 3.53 * f, "3-body": 4.48 * f}[channel] * EV_CGS


def formation_energy(on=True):
    """Energy a formation reaction gives the gas per H2 formed, in erg: the binding energy, of which
    formation_radiation() takes back what is radiated"""
    return H2_BINDING_ENERGY if on else 0


def formation_radiated(channel):
    """Binding energy radiated per H2 formed through channel, in erg"""
    return H2_BINDING_ENERGY - formation_heat(channel)


def formation_radiation(reaction, channel):
    """Cooling term of the energy radiated by the H2 that reaction forms through channel"""
    return ThermalTerm(-formation_radiated(channel) * reaction.unclumped().rate, name=RADIATION_OF + reaction.name,
                       bibliography=BIBLIOGRAPHY, clumping=reaction.clumping)


def dissociation_heat(on=True):
    """Heat given to the gas per collisional dissociation of H2, in erg"""
    return -H2_BINDING_ENERGY if on else 0
