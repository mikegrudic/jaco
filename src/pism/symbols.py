"""Definition of various sympy symbols used throughout the module"""

import sympy as sp
from typing import Union

T = sp.Symbol("T", nonnegative=True)  # temperature
T5 = T / 10**5
T6 = T / 10**6
T3 = T / 10**3
T4 = T / 10**4
c_s = sp.Symbol("c_s", nonnegative=True)  # sound speed
G = sp.Symbol("G", nonnegative=True)  # gravitational constant
ρ = sp.Symbol("ρ", nonnegative=True)  # total mass density
n_e = sp.Symbol("n_e-", nonnegative=True)  # electron number density
z = sp.Symbol("z", nonnegative=True)  # cosmological redshift
t = sp.Symbol("t")  # time
dt = sp.Symbol("Δt", nonnegative=True)


def d_dt(species: Union[str, sp.core.symbol.Symbol]):
    if isinstance(species, str):
        return sp.diff(sp.Function(n_(species))(t), t)
    else:
        return sp.diff(sp.Function(species)(t), t)


def n_(species: str):
    match species:
        case "heat":
            return sp.Symbol(f"⍴u", nonnegative=True)
        case _:
            return sp.Symbol(f"n_{species}", nonnegative=True)


def BDF(species):
    return (n_(species) - sp.Symbol(str(n_(species)) + ",0")) / dt
