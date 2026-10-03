"""Heating and cooling terms with no chemistry"""

import sympy as sp
from ..process import Process
from ..symbols import c_s, G, ρ, T, n_e, z, n_


class ThermalTerm(Process):
    """A heating (positive) or cooling (negative) term of the gas energy, immutable after construction.

    Parameters
    ----------
    heat: expression
        Energy given to the gas per unit volume and time (erg cm^-3 s^-1), before clumping.
    name: str, optional
    bibliography: sequence of str, optional
    clumping: expression, optional
        Factor <n^k>/<n>^k multiplying heat (e.g. C_2 for a collisional rate); 1 by default.
    reservoir: str, optional
        Energy reservoir (e.g. "dust heat") that receives what the gas loses, so the two rows sum to zero.
    """

    def __init__(self, heat, name="", bibliography=(), *, clumping=1, reservoir=None):
        self.clumping = clumping
        self.reservoir = reservoir
        self._heat_unclumped = heat
        total = heat if clumping == 1 else heat * clumping
        rows = {"heat": total}
        if reservoir is not None:
            rows[reservoir] = -total
        super().__init__(name, bibliography, rows)

    def unclumped(self):
        if self.clumping == 1:
            return self
        return ThermalTerm(self._heat_unclumped, self.name, self.bibliography, reservoir=self.reservoir)


def collisional_thermal_term(colliders, heat_rate_coefficient, *, clumping=None, name="", bibliography=()):
    """Heat ``heat_rate_coefficient * prod(n_collider) * C_k`` of k colliders (a sequence, repeated for like
    colliders); ``heat_rate_coefficient`` is negative for cooling. clumping defaults to C_k, or 1 for one collider."""
    colliders = tuple(colliders)
    if clumping is None:
        clumping = sp.Symbol(f"C_{len(colliders)}") if len(colliders) > 1 else 1
    return ThermalTerm(heat_rate_coefficient * sp.prod([n_(c) for c in colliders]), name, bibliography,
                       clumping=clumping)


ThermalProcess = ThermalTerm

PdV_heating = ThermalTerm(
    sp.Symbol("C_1") * c_s**2 * sp.sqrt(4 * sp.pi * G * ρ),
    name="Grav. Compression",
    bibliography=["1998ApJ...495..346M"],
)

inv_compton_cooling = ThermalTerm(
    -5.41e-36 * n_e * T * (1 + z) ** 4, name="Inverse Compton Cooling", bibliography=["1986ApJ...301..522I"]
)
