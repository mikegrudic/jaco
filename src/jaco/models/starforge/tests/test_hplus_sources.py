"""H+ is made only as in GIZMO's ionization balance (find_abundances_and_rates: nH0 = aHp / (aHp + geH0 + gJH0ne)):
by collisions, and not by cosmic rays, whose electrons GIZMO counts as metal and molecular ions (metal_electrons.py).
GIZMO's photoionization term (UV background, RT) has no counterpart in this model."""

import sympy as sp
from jaco.symbols import n_
from ..starforge import make_model


def test_only_collisions_make_Hplus():
    sources = [p.name for p in make_model().subprocesses
               if "H+" in p.network and sp.simplify(p.network["H+"].rhs.subs(n_("H+"), 0)) != 0]
    assert sources == ["Collisional Ionization of H"]
