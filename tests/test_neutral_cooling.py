from pism.processes import (
    FreeFreeEmission,
    LineCoolingSimple,
    CollisionalIonization,
    GasPhaseRecombination,
    Ionization,
)
import numpy as np
import matplotlib

matplotlib.use("Agg")
from matplotlib import pyplot as plt
import sympy as sp
from pism.data import SolarAbundances


def test_neutral_cooling():
    """Solve for thermochemical structure of neutral ISM at a range of densities. This can potentially change
    as we update the data or process implementations...
    """

    processes = (
        [CollisionalIonization(s) for s in ("H", "He", "He+")]
        + [GasPhaseRecombination(i) for i in ("H+", "He+", "He++")]
        + [FreeFreeEmission(i) for i in ("H+", "He+", "He++")]
        + [LineCoolingSimple(i) for i in ("H", "He+")]
    )
    system = sum(processes)

    # throw in some order-of-magntiude numbers to model photoelectric heating and cosmic-ray ionization
    heat_per_H = 1e-27
    zeta_CR = 2e-16
    system.heat += heat_per_H * sp.Symbol("n_Htot")
    system += Ionization(species="H", rate_per_volume=zeta_CR * sp.Symbol("n_Htot"))
    T = sp.Symbol("T")

    # C+ cooling - assumes x_C+ = x_e- = 3e-4
    x_C = 3e-4
    system.heat -= (
        1e-27
        * sp.Symbol("n_Htot")
        * sp.exp(-91.211 / T)
        * (4890 / sp.sqrt(T) * (x_C * sp.Symbol("n_H")) + 0.47 * T**0.15 * sp.Symbol("n_Htot"))
    )

    ngrid = np.logspace(-3, 3, 10**4)

    knowns = {"n_Htot": ngrid}
    y = SolarAbundances.x("He")
    guesses = {
        "T": 100.0 * np.ones_like(ngrid),
        "H": ngrid * 0.99,  # density
        "He": y * ngrid * 0.99,  # density
        "He+": y * ngrid * 0.01,  # density
    }

    sol = system.steadystate(
        knowns, guesses, tol=1e-3, input_abundances=False, output_abundances=True, careful_steps=30
    )

    fig, ax = plt.subplots(figsize=(3, 3))
    ax.loglog(ngrid, sol["T"], label=r"T (K)", color="black")
    ax.loglog(ngrid, sol["H"], label=r"$x_{\rm H0}$", ls="dotted", color="black")
    ax.loglog(ngrid, sol["e-"], label=r"$x_{e-}$", ls="dashed", color="black")
    ax.set_xlabel(r"$n_{\rm H}\,\left(\rm cm^{-3}\right)$")
    ax.legend(labelspacing=0)
    ax.set_xlim(ngrid[0], ngrid[-1])
    ax.set_ylim(1e-4, 3e4)
    ax.set_yticks(10.0 ** np.arange(-4, 5))
    plt.savefig("tests/neutral_cooling.png", bbox_inches="tight")
    #    np.save("neutral_cooling_testdata.npy", np.c_[ngrid, sol["T"]])
    T_test = np.load("tests/neutral_cooling_testdata.npy")[:, 1]
    assert np.all(np.abs((T_test - sol["T"]) / sol["T"]) < 0.1)
