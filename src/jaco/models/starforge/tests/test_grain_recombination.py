"""Grain-assisted recombination coefficient against GIZMO's alpha_recomb_grain (simple_chemistry.cc, Weingartner &
Draine 2001): per ion and H nucleus, with GIZMO's grain charging parameter G0 sqrt(T)/n_e + 50."""

import pytest
import sympy as sp
from ..grain_assisted_recombination import alpha_grain
from ..symbols import T, G_0, Z_dust, f_dust, n_Htot
from . import gizmo_electrons as G

x_e = sp.Symbol("x_e")


@pytest.mark.parametrize("ion", ["H+", "C+"])
@pytest.mark.parametrize("Tv", [10.0, 100.0, 1e3, 8e3])
@pytest.mark.parametrize("G0,nH,xe", [(1.7, 1.0, 1e-3), (0.3, 30.0, 2e-4), (1e-4, 1e4, 1e-7), (0.0, 1e2, 1e-6)])
@pytest.mark.parametrize("Zsol", [1.0, 0.1])
def test_matches_gizmo(ion, Tv, G0, nH, xe, Zsol):
    alpha = sp.lambdify((T, G_0, Z_dust, f_dust, n_Htot, x_e), alpha_grain(ion, x_e), modules="numpy")
    expected = G.alpha_recomb_grain(ion, Tv, xe, nH, G0, [Zsol * z for z in G.SOLAR])
    assert alpha(Tv, G0, Zsol, 1.0, nH, xe) == pytest.approx(expected, rel=1e-12)
