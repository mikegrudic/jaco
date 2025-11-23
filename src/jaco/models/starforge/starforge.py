"""Specification of STARFORGE network including radiation, thermo, cosmic rays, and dust"""

from jaco.processes import CollisionalIonization, GasPhaseRecombination, FreeFreeEmission, LineCoolingSimple
from .CO_cooling import CO_cooling
from ..model import Model


def make_model():
    atoms = "H", "He", "C"
    ions = "H+", "He+", "He++", "C+"
    molecules = "H2",

    processes = (
        [CollisionalIonization(s) for s in ("H", "He", "He+")]
        + [GasPhaseRecombination(i) for i in ("H+", "He+", "He++")]
        + [FreeFreeEmission(i) for i in ("H+", "He+", "He++")]
        + [LineCoolingSimple(i) for i in ("H", "He+", "C+")]
    )

    processes += gizmo_H2_chemistry
    processes += 
    processes += CO_cooling(prescription="Whitworth+2018")

    processes += [photon_absorption(band) for band in "EUV", "FUV", "NUV", "OPT", "FIR"]
    processes += [dust_emission(band) for band in "EUV", "FUV", "NUV", "OPT", "FIR"]
    
