"""Definition of symbols that are used throughout the model"""

import sympy as sp
from ...symbols import x_, n_

T = sp.Symbol("T")  # Gas temperature
sqrt_T = sp.sqrt(T)
log_T = sp.log(T, 10)
n_Htot = n_("Htot")  # Total number density of H nuclei
X_H = sp.Symbol("X")  # mass fraction of hydrogen
T_dust = sp.Symbol("Td")  # Dust temperature
f_dust = sp.Symbol("f_d")  # Factor accounting for sublimation
Z_dust = sp.Symbol("Z_d")  # Solar-normalized dust abundance. Value of 1 corresponds to Solar neighborhood dust.
G_0 = sp.Symbol("G_0")  # UV radiation field normalized to Habing
f_shield = sp.Symbol("f_shield")  # Lyman-Werner self-shielding factor
grad_v = sp.Symbol("∇v")  # velocity gradient Frobenius norm
NH = sp.Symbol("N_H")  # column density of H nuclei
dx = sp.Symbol("Δx")  # effective cell size
ISRF = sp.Symbol("ISRF")  # scaling factor for ISRF and cosmic ray background
H2_formation_heat_cgs = 7.2e-12
