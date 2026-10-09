"""GIZMO's optically-thick cooling limit in the starforge model (jaco.models.starforge.optically_thick)."""

import numpy as np
import pytest
import sympy as sp
from jaco.models.starforge import optically_thick as ot
from jaco.models.starforge.symbols import T, n_Htot, X_H
from jaco.models.starforge.dust_balance import N_H_cool

PROTONMASS = 1.67262178e-24


def gizmo_limit(n_Htot, X, T, N_H, kappa):
    """cooling.cc 1460-1483: Lambda_Thick_BlackBody nHcgs^2, nHcgs = 0.76 rho / m_p, Sigma = N_H m_p"""
    nHcgs = 0.76 * n_Htot / X  # 0.76 rho / m_p with rho = n_Htot m_p / X
    Sigma = N_H * PROTONMASS
    effective_area = 2.3 * PROTONMASS / Sigma
    return 5.67e-5 * T**4 * effective_area / ((1 + kappa * Sigma) * nHcgs) * nHcgs**2


def test_limit_is_gizmos():
    L = sp.lambdify((n_Htot, X_H, T, N_H_cool, ot.kappa_thick), ot.thick_limit(), "numpy")
    grid = np.meshgrid(np.logspace(0, 18, 7), [0.7155], np.logspace(0.5, 4, 6), np.logspace(18, 27, 5),
                       np.logspace(-1, 2, 4), indexing="ij")
    assert np.allclose(L(*grid), gizmo_limit(*grid), rtol=1e-12, atol=0)


def test_smooth_cap_tracks_the_hard_cap():
    Q, L = sp.symbols("Q L", positive=True)
    capped = sp.lambdify((Q, L), Q * ot.smooth_cap_factor(Q, L), "numpy")
    r = np.logspace(-6, 100, 20001)
    ratio = np.array([capped(x, 1.0) for x in r]) / np.minimum(r, 1)
    assert np.all(ratio <= 1) and np.all(np.isfinite(ratio))
    assert ratio.min() >= 2 ** (-1 / 16) - 1e-12  # the worst case is |Q| = L
    assert capped(1.0, 1.0) == pytest.approx(2 ** (-1 / 16), rel=1e-12)
    far = (r < 1 / 1.5) | (r > 1.5)
    assert ratio[far].min() > 0.999
    assert capped(-3.0, 1.0) == pytest.approx(-1.0, rel=1e-3)  # heating is capped too, as in GIZMO


def test_limit_spares_chemistry_pdv_and_the_gas_dust_exchange():
    """The H2/H- reactions, PdV work and the gas-dust collisions (heat moved between gas and dust) are not capped; the
    gas's radiation and the dust's (absorption and emission, through its dust heat row) are"""
    from jaco.models.starforge.starforge import make_model
    from jaco.models.starforge.h2_chemistry import h2_chemistry_processes
    from jaco.models.starforge.dust_balance import DUST_RADIATION
    m = make_model()
    chem = {p.name for p in h2_chemistry_processes()}
    capped = set(ot.capped_ids(m, chem | {"Gas-dust collisions"}))
    assert ot.LIMIT_ID in {r.name for r in m.rules} and "kappa_thick" in {p.name for p in m.parameters}
    assert not capped & (chem | {"PdV work", "Gas-dust collisions"})
    assert {"H2 + HD Line Cooling", "Collisional Ionization of H", *DUST_RADIATION} <= capped
    m_thin = make_model(optically_thick=False)
    assert ot.LIMIT_ID not in {r.name for r in m_thin.rules} and "kappa_thick" not in {p.name for p in m_thin.parameters}
    # every capped energy row carries the factor, nothing else does
    for i, p in {p.name: p for p in m.subprocesses}.items():
        for row in ot.ENERGY_ROWS:
            e = sp.sympify(p.network[row].rhs) if row in p.network else sp.S.Zero
            assert e == 0 or (ot.F_THICK in e.free_symbols) == (i in capped), (i, row)
    Q = dict(m.intermediates)[ot.Q_THICK]
    assert sp.Symbol("A_dust") in Q.free_symbols and sp.Symbol("Td") in Q.free_symbols
