"""One-zone closures, abundances and initial state of Menon et al. 2025's collapse test (KROME driver
popsicle_semenov_photo_cr_full_collapseTest/test.f90), for jaco.onezone.OneZone with this model.

    zone = OneZone(make_model(), params=parameters(Z), closures=PopsicleClosures.DEFAULT, host="plain")
    track = zone.free_fall(n0=0.1, T0=300.0, n_stop=1e17, x0=INITIAL_ABUNDANCES)

The driver takes every column as num2col, half a Jeans length sqrt(pi k T / (G rho mu m_p)) of the species, the
shielding column N_H = n_Htot lambda_J / 2, A_V = Z N_H / 1.87e21, and the dust optical depth over a full lambda_J."""

from jaco.onezone import Closures
from ..starforge.symbols import x_solar

# x(H2) = 2e-3 n, x(H+) = x(e-) = 1e-4 n (plus C+ = 1.6e-4 Z n, which starforge's carbon partition sets), He = 0.1 n
INITIAL_ABUNDANCES = {"H+": 1e-4, "H_2": 2e-3}
HE_PER_H = 0.1
C_GAS_SOLAR, O_GAS_SOLAR = 1.6e-4, 3.2e-4


def parameters(Z):
    """Parameter overrides at metallicity Z (solar units): dust-to-gas Z, gas-phase C and O as the driver's (starforge's
    gas-phase C is x_C,tot times 1.6e-4 / x_solar(C)), no other metals, He/H = 0.1"""
    p = {"Z_d": Z, "x_C,tot": x_solar("C") * Z, "x_O,tot": O_GAS_SOLAR * Z, "y": HE_PER_H,
         "X": 1.0 / (1.0 + 4.0 * HE_PER_H)}
    for el in ("N", "Ne", "Mg", "Si", "S", "Ca", "Fe"):
        p[f"x_{el}"] = 0.0
    return p


class PopsicleClosures(Closures):
    """The bare closures (columns of half a Jeans length, PdV P dln rho/dt) plus the Jeans length itself"""

    @staticmethod
    def lambda_J(st):
        return st["lambda_J"]


PopsicleClosures.DEFAULT = {**Closures.DEFAULT, "lambda_J": PopsicleClosures.lambda_J}

