from .eos import species_heat_capacity, species_energy, k_B_cgs
import numpy as np
import os
import sympy as sp


def test_H2_energy_fit():
    """Checks that the fit for H_2 energy is with 1% of the exactly-computed value"""
    Tgrid, ebeta_exact, cv_exact = np.load(os.path.dirname(os.path.abspath(__file__)) + "/T_vs_ebeta_cv.npy").T
    T = sp.Symbol("T")
    numfunc = sp.lambdify(T, species_energy("H_2") / k_B_cgs / T)
    ebeta = numfunc(Tgrid)
    assert np.all(np.isclose(ebeta, ebeta_exact, atol=1e-2))

    numfunc = sp.lambdify(T, species_heat_capacity("H_2") / k_B_cgs)
    cv = numfunc(Tgrid)
    assert np.all(np.isclose(cv, cv_exact, atol=0.03))
