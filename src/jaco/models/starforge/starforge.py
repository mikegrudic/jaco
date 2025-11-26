"""Specification of STARFORGE network including radiation, thermo, cosmic rays, and dust"""

from jaco.processes import CollisionalIonization, GasPhaseRecombination, FreeFreeEmission, LineCoolingSimple
from ..model import Model
import sympy as sp
# import h2_chemistry


def make_model():
    atoms = "H", "He", "C"
    ions = "H+", "He+", "He++", "C+"
    molecules = ("H2",)

    processes = (
        [CollisionalIonization(s) for s in ("H", "He", "He+")]
        + [GasPhaseRecombination(i) for i in ("H+", "He+", "He++")]
        + [FreeFreeEmission(i) for i in ("H+", "He+", "He++")]
        + [LineCoolingSimple(i) for i in ("H", "He+", "C+")]
    )


#     processes += sum(h2_chemistry.reactions)
# #    processes += CO_cooling(prescription="Whitworth+2018")

#     processes += [photon_absorption(band) for band in "EUV", "FUV", "NUV", "OPT", "FIR"]
#     processes += [dust_emission(band) for band in "EUV", "FUV", "NUV", "OPT", "FIR"]

# assumptions: x_D = 2.527e-5 x_H, x_D+ = 2.527e-5 x_H+, zero out associated collisional dissociation rates
