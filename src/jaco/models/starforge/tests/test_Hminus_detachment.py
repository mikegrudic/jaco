"""H- collisional detachment (GA08 / Abel+97 fits in T_eV) against GIZMO's k15/k16 in update_explicit_molecular_fraction."""

import numpy as np
import pytest
import sympy as sp
from ..h2_chemistry.collisional_detachment import Hminus_collisional_detachment
from ..symbols import T

EXPMAX = 90.0


def gizmo_k(collider, Tval):
    lnTeV = np.log(Tval) - 9.35915
    if collider == "e-":  # k15
        c = [-1.801849334e1, 2.36085220e0, -2.827443e-1, 1.62331664e-2, -3.36501203e-2, 1.17832978e-2,
             -1.65619470e-3, 1.06827520e-4, -2.63128581e-6]
        return np.exp(max(-EXPMAX, sum(ci * lnTeV**i for i, ci in enumerate(c))))
    if Tval < 1160.45:  # k16
        return 1.46629e-16 * Tval**1.78186
    c = [-2.0372609e1, 1.13944933e0, -1.4210135e-1, 8.4644554e-3, -1.4327641e-3, 2.0122503e-4, 8.6639632e-5,
         -2.5850097e-5, 2.4555012e-6, -8.0683825e-8]
    return np.exp(max(-EXPMAX, sum(ci * lnTeV**i for i, ci in enumerate(c))))


@pytest.mark.parametrize("collider", ["H", "e-"])
@pytest.mark.parametrize("Tval", [30.0, 100.0, 300.0, 1000.0, 1160.0, 1161.0, 3000.0, 1e4, 3e4, 1e5])
def test_matches_gizmo(collider, Tval):
    expected = gizmo_k(collider, Tval)
    if expected <= np.exp(-EXPMAX):
        pytest.skip("below GIZMO's exp(-90) floor")
    k = sp.lambdify(T, Hminus_collisional_detachment(collider).rate_coefficient, modules="numpy")
    assert k(Tval) == pytest.approx(expected, rel=1e-10, abs=0)
