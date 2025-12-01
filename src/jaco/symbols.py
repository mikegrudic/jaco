"""Definition of various sympy symbols used throughout the module"""

import sympy as sp
from typing import Union
import sympy as sp
from astropy.constants import k_B, m_p
from astropy import units as u
import numpy as np

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
    """Returns number density symbol for species"""
    match species:
        case "heat":
            return sp.Symbol(f"⍴u")
        case "dust heat":
            return sp.Symbol(f"⍴u_dust")
        case _:
            return sp.Symbol(f"n_{species}")


def x_(species: str):
    """Returns abundance symbol for species"""
    return sp.Symbol(f"x_{species}")


def BDF(species):
    # if species in ("T", "u"):  # this is the heat equation
    #     return rho * (internal_energy - sp.Symbol("u_0")) / dt
    # else:
    return (n_(species) - sp.Symbol(str(n_(species)) + "_0")) / dt


def sanitize_symbols(expr):
    """
    Given a symbolic expression, replace all symbols with + or - in them with plus or minus
    to avoid ambiguity or syntax errors in code generation.
    """
    if hasattr(expr, "__iter__"):  # if iterable call recursively until we get down to actual expressions
        expr = [sanitize_symbols(e) for e in expr]
    else:
        for s in expr.free_symbols:
            if "+" in str(s):
                expr = expr.subs(s, sp.Symbol(str(s).replace("+", "plus")))
            if "-" in str(s):
                expr = expr.subs(s, sp.Symbol(str(s).replace("-", "minus")))
    return expr


def piecewise_linear(X, Y, x, extrapolate=False):
    """
    Returns a symbolic expression describing a piecewise-linear interpolant of data X and Y

    Parameters
    ----------
    X: array_like
        array of X values
    Y: array_like
        array of Y values
    x: sympy.Symbol
        Symbol for the independent variable associated with X
    extrapolate: Boolean, optional
        Whether to use a linear extrapolant past the limits, or simply limit to the outermost values
    """
    if not np.all(np.diff(X) > 0):
        raise ValueError("X values passed to piecewise_linear must be monotonic.")
    if len(X) != len(Y):
        raise ValueError("X and Y must have the same length.")

    slopes = np.diff(Y) / np.diff(X)
    cases = []
    if extrapolate:
        cases.append((Y[0] + slopes[0] * (x - X[0]), x < X[0]))
        cases.append((Y[-1] + slopes[-1] * (x - X[-1]), x >= X[-1]))
    else:
        cases.append((Y[0], x < X[0]))
        cases.append((Y[-1], x >= X[-1]))

    for i in range(len(X) - 1):
        cases.append((Y[i] + slopes[i] * (x - X[i]), (x >= X[i]) & (x < X[i + 1])))
    print(cases)
    return sp.Piecewise(*cases)


def piecewise_powerlaw(X, Y, x, extrapolate=False):
    """
    Returns a symbolic expression describing a piecewise-powerlaw interpolant of data X and Y

    Parameters
    ----------
    X: array_like
        array of X values
    Y: array_like
        array of Y values
    x: sympy.Symbol
        Symbol for the independent variable associated with X
    extrapolate: Boolean, optional
        Whether to use a powerlaw extrapolant past the limits, or simply limit to the outermost values
    """
    return sp.exp(piecewise_linear(np.log(X), np.log(Y), sp.log(x), extrapolate))
