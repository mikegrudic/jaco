"""CO rotational cooling, Whitworth & Jaffa (2018, A&A 611, A20; WJ18), Eqs. 15-18.

WJ18 calibrate on Goldsmith & Langer (1978) for gas with all of its H in H2: their n_H2 = rho / m_H2 stands for the mass
density (Sec. 1) and their X_CO is n_CO / n_H2. Here that n_H2 is n_Htot / 2 in the rates and in X_CO alike, whatever
the molecular fraction of the gas being cooled; its CO abundance is the model's (carbon_abundances). WJ18's |div v| is
the model's velocity-gradient norm.
"""

import sympy as sp
from jaco.processes import ThermalTerm
from jaco.symbols import n_
from .symbols import T, grad_v, x_, n_Htot, cmb_bath_factor, lowtemp_truncation

KMS_PER_PC = 1e5 / 3.085678e18  # 1 km/s/pc in s^-1
# Eq. 14
LAMBDA_LO = 2.16e-27  # erg s^-1
LAMBDA_HI = 2.21e-28  # erg s^-1
BETA_0, BETA_NH2, BETA_T = 1.23, 0.0533, 0.164

n_H2_wj18 = n_Htot / 2  # WJ18's n_H2: all H nuclei in H2
X_CO = x_("CO") * n_Htot / n_H2_wj18

# Eqs. 15-17 divided by n_H2, so per n_CO n_H2 (cm^3 erg s^-1)
lambda_CO_lo = LAMBDA_LO * T**1.5
lambda_CO_hi = LAMBDA_HI * (X_CO * KMS_PER_PC / grad_v) ** -1 * T**4 / n_H2_wj18**2
beta = BETA_0 * n_H2_wj18**BETA_NH2 * T**BETA_T
lambda_CO = (lambda_CO_lo ** (-1 / beta) + lambda_CO_hi ** (-1 / beta)) ** -beta  # Eq. 18

CO_cooling = ThermalTerm(
    -lambda_CO * n_("CO") * n_H2_wj18 * cmb_bath_factor * lowtemp_truncation,
    name="CO Cooling",
    bibliography=["2018A&A...611A..20W"],
    clumping=sp.Symbol("C_2"),
)
