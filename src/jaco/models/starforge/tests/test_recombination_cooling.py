"""H/He recombination rates and cooling against GIZMO's MakeCoolingTable() + find_abundances_and_rates()."""

import numpy as np
import pytest
import sympy as sp
from jaco.processes import GasPhaseRecombination
from jaco.symbols import T


def gizmo_rates(Tv):
    aHp = 7.982e-11 / (np.sqrt(Tv / 3.148) * (1 + np.sqrt(Tv / 3.148)) ** 0.252 * (1 + np.sqrt(Tv / 7.036e5)) ** 1.748)
    aHep = 9.356e-10 / (np.sqrt(Tv / 4.266e-2) * (1 + np.sqrt(Tv / 4.266e-2)) ** 0.2108 * (1 + np.sqrt(Tv / 3.676e7)) ** 1.7892)
    aHepp = 2 * 7.982e-11 / (np.sqrt(Tv / (4 * 3.148)) * (1 + np.sqrt(Tv / (4 * 3.148))) ** 0.252
                            * (1 + np.sqrt(Tv / (4 * 7.036e5))) ** 1.748)
    ad = 1.9e-3 * Tv**-1.5 * np.exp(-470000 / Tv) * (1 + 0.3 * np.exp(-94000 / Tv)) if 470000.0 / Tv < 70 else 0.0
    rates = {"H+": aHp, "He+": aHep + ad, "He++": aHepp}
    # LambdaRecHp, LambdaRecHep + LambdaRecHepd, LambdaRecHepp per n_e n_ion
    cooling = {"H+": 1.036e-16 * Tv * aHp, "He+": 1.036e-16 * Tv * aHep + 6.526e-11 * ad, "He++": 1.036e-16 * Tv * aHepp}
    return rates, cooling


@pytest.mark.parametrize("ion", ["H+", "He+", "He++"])
@pytest.mark.parametrize("Tv", [1e3, 1e4, 3e4, 1e5, 3e5, 1e6, 1e7])
def test_matches_gizmo(ion, Tv):
    process = GasPhaseRecombination(ion)
    rates, cooling = gizmo_rates(Tv)
    assert float(process.rate_coefficient.subs(T, Tv)) == pytest.approx(rates[ion], rel=1e-10, abs=0)
    assert -float(process.heat_rate_coefficient.subs(T, Tv)) == pytest.approx(cooling[ion], rel=1e-10, abs=0)
