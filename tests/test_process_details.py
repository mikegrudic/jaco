"""Small invariants of individual processes and helpers."""

import sympy as sp

from jaco.processes.recombination import Recombination
from jaco.symbols import n_


def test_recombination_colliders_are_ordered():
    process = Recombination("H+")
    assert process.colliding_species == ("H+", "e-")
    assert process.nprod == n_("H+") * n_("e-")
    assert process.clumping_factor == sp.Symbol("C_2")
