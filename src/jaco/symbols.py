"""Definition of various sympy symbols used throughout the module"""

import sympy as sp
from typing import Union
import sympy as sp
from astropy.constants import k_B, m_p
from astropy import units as u

T = sp.Symbol("T")  # temperature
T5 = T / 10**5
T6 = T / 10**6
T3 = T / 10**3
T4 = T / 10**4
c_s = sp.Symbol("c_s")  # sound speed
G = sp.Symbol("G")  # gravitational constant
ρ = sp.Symbol("ρ")  # total mass density
n_e = sp.Symbol("n_e-")  # electron number density
z = sp.Symbol("z")  # cosmological redshift
t = sp.Symbol("t")  # time
dt = sp.Symbol("Δt")
n_Htot = sp.Symbol("n_Htot")


boltzmann_cgs = k_B.to(u.erg / u.K).value
protonmass_cgs = m_p.to(u.g).value
# write down internal energy density in terms of number densities - this defines the EOS


def d_dt(species: Union[str, sp.core.symbol.Symbol]):
    if isinstance(species, str):
        return sp.diff(sp.Function(n_(species))(t), t)
    else:
        return sp.diff(sp.Function(species)(t), t)


def n_(species: str):
    match species:
        case "heat":
            return sp.Symbol(f"⍴u")
        case _:
            return sp.Symbol(f"n_{species}")


egy_density = boltzmann_cgs * T * (1.5 * (n_("e-") + n_("H") + n_("H+") + n_("He") + n_("He+") + n_("He++")))
rho = protonmass_cgs * (n_("H") + n_("H+") + 4 * (n_("He") + n_("He+") + n_("He++")))
internal_energy = sp.factor(egy_density / rho)


def x_(species: str):
    return sp.Symbol(f"x_{species}")


def BDF(species):
    if species in ("T", "u"):  # this is the heat equation
        return rho * (internal_energy - sp.Symbol("u_0")) / dt
    else:
        return (n_(species) - sp.Symbol(str(n_(species)) + "_0")) / dt
