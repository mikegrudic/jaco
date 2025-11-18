from jaco.processes import CollisionalIonization, GasPhaseRecombination
import numpy as np
from matplotlib import pyplot as plt


def test_CIE():
    """Solve for CIE of a H-He mixture and check that the He+ abundance looks reasonable"""
    chianti_data = np.load("tests/chianti_He_abundances.npy")
    N = len(chianti_data)
    Tmin, Tmax = chianti_data[:, 0].min(), chianti_data[:, 0].max()

    processes = [CollisionalIonization(s) for s in ("H", "He", "He+")] + [
        GasPhaseRecombination(i) for i in ("H+", "He+", "He++")
    ]
    system = sum(processes)

    Tgrid = np.logspace(np.log10(Tmin), np.log10(Tmax), N)
    ngrid = np.ones_like(Tgrid)

    knowns = {"T": Tgrid, "n_Htot": ngrid}

    guesses = {"H": 0.5 * np.ones_like(Tgrid), "He": 1e-5 * np.ones_like(Tgrid), "He+": 1e-5 * np.ones_like(Tgrid)}
    sol = system.solve(knowns, guesses, tol=1e-3)
    Hep_frac = sol["He+"] / (sol["He"] + sol["He+"] + sol["He++"])
    assert (
        np.abs(Hep_frac - chianti_data[:, 2]).max() < 0.2
    )  # just a loose test because the codes don't use the same rates


def generate_CIE_testdata():
    import ChiantiPy.core as ch

    Tgrid = np.logspace(3, 6, 10**4)
    he = ch.ioneq(2)
    he.load()
    he.calculate(Tgrid)
    np.save("chianti_He_abundances.npy", np.c_[Tgrid, he.Ioneq])
