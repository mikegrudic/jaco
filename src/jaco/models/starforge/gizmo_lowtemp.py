"""GIZMO's low-temperature metal and molecular cooling (COOL_LOW_TEMPERATURES block of CoolingRate, detailed branch
for COOL_METAL_LINES_BY_SPECIES and FIRE-3), for STARFORGE_LEGACY.

The rates are written as GIZMO's per-n_H^2 coefficients times n_Htot^2 with n_Htot as the density; the model builder
replaces n_Htot by GIZMO's nHcgs (switch rate_density). Abundances are per H nucleus: x_H0 = 1 - x_H+ counts H2 nuclei
as neutral H, as GIZMO's nH0 does, and GIZMO's f_molec is x_H2.
"""

import sympy as sp
from jaco.processes import ThermalProcess
from .symbols import (T, log_T, sqrt_T, x_, n_Htot, G_0, NH, grad_v, X_H, x_solar, cmb_bath_factor,
                      lowtemp_truncation, PROTONMASS_CGS)
from .H2_cooling import lambda_H2_thin, Lambda_HD_thin, Lambda_H2_LTE_per_molecule

xH0 = 1 - x_("H+")
x_e = x_("e-")
KMS_PER_PC = 1e5 / 3.085678e18  # 1 km/s/pc in s^-1
MIN_REAL_NUMBER = 1e-56  # GIZMO constants.h


def f_Cplus_CCO(n=n_Htot):
    """C+ share of the C+ / CO cooling interpolation: fco/(1 - fco) ~ (n / 340 G0)^2 / sqrt(T) (Tielens)"""
    return 1 / (1 + (n / (340 * sp.Max(sp.Rational(1, 10), G_0))) ** 2 / sqrt_T)


def gizmo_carbon_cooling_coefficient():
    """Lambda_Metals = f_Cplus_CCO Lambda_Cplus + Lambda_CO per n_H0 n_H (erg cm^3 s^-1): C+ fine structure (Barinovs+05,
    Wilson & Bell 02) plus [CI] 609 um (Hocuk+16), both weighted by the C+ share, and HM79 CO with GIZMO's
    recalibration and column correction, capped by the Whitworth & Jaffa LVG rate"""
    Z_C = sp.Max(1e-6, sp.Symbol("x_C,tot") / x_solar("C"))
    column = PROTONMASS_CGS * NH  # g cm^-2
    Lambda_Cplus = Z_C * (4.7e-28 * (T**0.15 + 1.04e4 * x_e / sqrt_T) * sp.exp(-91.211 / T)
                          + 2.08e-29 * sp.exp(-23.6 / T))
    ncrit_CO, Sigma_crit_CO = 1.9e4 * sqrt_T, 3.0e-5 * T / Z_C
    Lambda_CCO = Z_C * T * sqrt_T * 2.73e-31 / (1 + (n_Htot / ncrit_CO) * (1 + sp.Max(column, 0.017) / Sigma_crit_CO))
    Lambda_CO_HI = 4.42e-28 * (grad_v / KMS_PER_PC) * T**4 / n_Htot**2
    f = f_Cplus_CCO()
    return f * Lambda_Cplus + sp.Min(Lambda_CO_HI, (1 - f) * Lambda_CCO)


gizmo_carbon_cooling = ThermalProcess(
    -n_Htot**2 * xH0 * gizmo_carbon_cooling_coefficient() * cmb_bath_factor * lowtemp_truncation,
    name="GIZMO C+, [CI] and CO cooling",
    bibliography=["2005ApJ...620..537B", "2002MNRAS.337.1027W", "2016MNRAS.456.2586H", "1979ApJS...41..555H",
                  "2018A&A...611A..20W"],
)


def gizmo_H2_cooling_coefficient():
    """Lambda_H2 + Lambda_HD per n_H^2 (erg cm^3 s^-1): GA08 thin rates with GIZMO's collider weights (abundance times
    the fixed mass fraction X_H, and Y_He for He), the HM79 LTE limit, and HD/H2 = min(0.00126, 4e-5 x_H0 / x_H2)"""
    x_H2 = x_("H_2")
    Y_He = 4 * sp.Symbol("y") * X_H  # y = n_He / n_H = Y / (4 X)
    thin = (sp.Max(xH0 - 2 * x_H2, 0) * X_H * lambda_H2_thin("H") + Y_He * lambda_H2_thin("He")
            + x_H2 * X_H * lambda_H2_thin("H_2") + x_("H+") * X_H * lambda_H2_thin("H+") + x_e * X_H * lambda_H2_thin("e-"))
    n_over_ncrit = thin * n_Htot / Lambda_H2_LTE_per_molecule()
    f_HD = sp.Min(0.00126 * x_H2, 4.0e-5 * xH0)
    return x_H2 * thin / (1 + n_over_ncrit) + f_HD * Lambda_HD_thin() / (1 + f_HD / (x_H2 + MIN_REAL_NUMBER) * n_over_ncrit)


gizmo_H2_cooling = ThermalProcess(
    -n_Htot**2 * gizmo_H2_cooling_coefficient() * cmb_bath_factor * lowtemp_truncation,
    name="GIZMO H2 + HD cooling",
    bibliography=["2008MNRAS.388.1627G", "1998A&A...335..403G", "1979ApJS...41..555H"],
)
