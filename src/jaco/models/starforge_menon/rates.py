"""Rate coefficients and shielding functions of the POPSICLE KROME network of Menon et al. (2025), as its reaction file
(react_popsicle_semenov_photo_cr_full) and its collapse-test driver write them. Numbers in the comments are that
network's reaction indices. T in K."""

import numpy as np
import sympy as sp
from ..starforge.symbols import T, log_T, n_Htot, x_
from ..starforge.h2_chemistry.associative_detachment import k2_new as k10_kreckel

EV = 1.602176634e-12
LN10 = sp.log(10.0)
ln_T = sp.log(T)
T_eV = T / 11604.5
lnTe = sp.log(T_eV)

# 9: H + e- -> H- + photon (de Jong 1972)
k9 = 1.4e-18 * T**0.928 * sp.exp(-sp.Min(T, 5e5) / 16200.0)  # exp held below 1e-13 above 5e5 K
# 10: H- + H -> H2 + e- (Kreckel et al. 2010)
k10 = k10_kreckel
# 21: H- + H+ -> 2H (Stenrup et al. 2009); clamped at 0 where the fit turns negative
k21 = sp.Max(2.96e-6 / sp.sqrt(T) - 1.73e-9 + 2.50e-10 * sp.sqrt(T) - 7.77e-13 * T, 0)
# 22: H- + H+ -> H2+ + e- (Poulaert 1978)
k22 = 1e-8 * T**-0.4
# 11, 12: H + H+ -> H2+ + photon (Ramaker & Peek 1976, Coppola)
k11 = sp.Piecewise((2.10e-20 * (T / 30.0) ** -0.15, T < 30),
                   (10 ** (-18.20 - 3.194 * log_T + 1.786 * log_T**2 - 0.2072 * log_T**3), True))
# 13: H2+ + H -> H2 + H+
k13 = 6.4e-10
# 14: H2 + H+ -> H2+ + H (Savin et al. 2004), 100 K <= T <= 3e4 K
_lnTs = sp.log(sp.Min(sp.Max(T, 100.0), 3e4))
k14 = sp.Max(-3.3232183e-7 + 3.3735382e-7 * _lnTs - 1.4491368e-7 * _lnTs**2 + 3.4172805e-8 * _lnTs**3
             - 4.7813728e-9 * _lnTs**4 + 3.9731542e-10 * _lnTs**5 - 1.8171411e-11 * _lnTs**6
             + 3.5311932e-13 * _lnTs**7, 0) * sp.exp(-21237.15 / T)
# 15: H2 + e- -> H + H- (Capitelli et al. 2007)
k15 = 35.5 * T**-2.28 * sp.exp(-46707.0 / T)
# 16: H2 + e- -> 2H + e- (Mitchell & Deveau 1983)
k16 = 4.38e-10 * sp.exp(-102000.0 / T) * T**0.35
# 23: H2+ + e- -> 2H (Coppola et al. 2011), zero above 23642.7 K
_Tc = sp.Min(T, 23642.692081)
k23 = sp.Max(1e6 * (4.2278e-14 - 2.3088e-17 * _Tc + 7.3428e-21 * _Tc**2 - 7.5474e-25 * _Tc**3
                    + 3.3468e-29 * _Tc**4 - 5.528e-34 * _Tc**5), 0)
# 153: H2+ + H2 -> H3+ + H; 162, 163: H3+ + e- -> H2 + H, 3H (Glover et al. 2010)
k153 = 2.24e-9 * (T / 300) ** 0.042 * sp.exp(-T / 46600.0)
k162 = 2.34e-8 * (T / 300) ** -0.52
k163 = 4.36e-8 * (T / 300) ** -0.52
# 25: 3H -> H2 + H (Forrey 2013); 26: 2H + H2 -> 2H2, the same / 8
k25 = 6e-32 * T**-0.25 + 2e-31 * T**-0.5
k26 = k25 / 8
# 2, 3: H+ + e- -> H (Abel et al. 1997, from Ferland et al. 1992)
k2 = sp.Piecewise(
    (3.92e-13 * T_eV**-0.6353, T <= 5500),
    (sp.exp(-28.61303380689232 - 0.7241125657826851 * lnTe - 0.02026044731984691 * lnTe**2
            - 0.002380861877349834 * lnTe**3 - 0.0003212605213188796 * lnTe**4 - 0.00001421502914054107 * lnTe**5
            + 4.989108920299513e-6 * lnTe**6 + 5.755614137575758e-7 * lnTe**7 - 1.856767039775261e-8 * lnTe**8
            - 3.071135243196595e-9 * lnTe**9), True))
# Cen 1992 H+ recombination cooling coefficient [erg cm^3 s^-1]
lambda_rec_Hplus = 8.7e-27 * sp.sqrt(T) * (T / 1e3) ** -0.2 / (1 + (T / 1e6) ** 0.7)

# Deuterium fractionation (Galli & Palla 2002; Glover & Jappsen 2007): 33/34 D+ + H2 <-> HD + H+ and 38/40 the charge
# exchange D + H+ <-> D+ + H
k33 = 1e-9 * sp.Max(0.417 + 0.846 * log_T - 0.137 * log_T**2, 1e-3)
k34 = 1e-9 * sp.exp(-457.0 / T)
k38 = sp.Max(2.00e-10 * T**0.402 * sp.exp(-37.1 / T) - 3.31e-17 * T**1.48, 1e-30)
k40 = 2.06e-10 * T**0.396 * sp.exp(-33.0 / T) + 2.03e-9 * T**-0.332


def _logaddexp(a, b):
    return sp.Max(a, b) + sp.log(1 + sp.exp(sp.Min(a, b) - sp.Max(a, b)))


# Martin et al. 1996 high-density H2 + H collision-induced dissociation and dissociative tunnelling (log10 fits)
_K_CID = [(-178.4239, -68.42243, 43.20243, -4.633167, 69.70086, 40870.38, -23705.70),
          (-142.7664, 42.70741, -2.027365, -0.2582097, 21.36094, 27535.31, -21467.79)]


def _ln_k_martin():
    lg = [c1 + c2 * log_T + c3 * log_T**2 + c4 * log_T**3 + c5 * sp.log(1 + c6 / T, 10) + c7 / T
          for c1, c2, c3, c4, c5, c6, c7 in _K_CID]
    return _logaddexp(LN10 * lg[0], LN10 * lg[1])


def _density_weight():
    """n/n_cr / (1 + n/n_cr) with Glover & Jappsen 2007's n_cr from the H and H2 critical densities (Lepp & Shull 1983,
    Shapiro & Kang 1987), the weight of the high-density limit"""
    logT4 = log_T - 4.0
    inv_ncr = (x_("H") * sp.exp(-LN10 * (3.0 - 0.416 * logT4 - 0.327 * logT4**2))
               + x_("H_2") * sp.exp(-LN10 * (4.845 - 1.3 * logT4 + 1.62 * logT4**2)))
    r = n_Htot * inv_ncr
    return r / (1 + r)


def k17_hunter():
    """17: H2 + H -> 3H (Hunter et al. 2023): Mac Low & Shull 1986 at low density, Martin et al. 1996 at high"""
    w = _density_weight()
    ln_low = sp.log(6.67e-12) + 0.5 * ln_T - (1.0 + 63590.0 / T)
    return sp.exp(w * _ln_k_martin() + (1 - w) * ln_low)


def k27_hunter():
    """27: H2 + H2 -> 2H + H2 (Hunter et al. 2023): Gong et al. 2017's low-density form, Martin et al. 1996 / 8 at
    high density"""
    w = _density_weight()
    ln_low = sp.log(5.996e-30) + 4.1881 * ln_T - 5.6881 * sp.log(1 + 6.761e-6 * T) - 54657.4 / T
    return sp.exp(w * (_ln_k_martin() - sp.log(8.0)) + (1 - w) * ln_low)


def h2heatfac():
    """1 / (1 + n_cr / n) with Hollenbach & McKee 1979's n_cr = 1e6 T^-1/2 / (1.6 x_H exp(-(400/T)^2) + 1.4 x_H2
    exp(-12000/(T + 1200))), the fraction of the rovibrational energy of new H2 that collisions turn into heat"""
    w = n_Htot * (1.6 * x_("H") * sp.exp(-((400 / T) ** 2)) + 1.4 * x_("H_2") * sp.exp(-12000 / (T + 1200)))
    return w / (w + 1e6 / sp.sqrt(T))


def f_shield_H2_WG11(N_H2, b5=1.0):
    """Wolcott-Green, Haiman & Bryan 2011 Eq. 12 at a fixed Doppler parameter b5 (km/s), as the driver calls it"""
    x = N_H2 / 5e14
    return 0.965 / (1 + x / b5) ** 1.1 + 0.035 / sp.sqrt(1 + x) * sp.exp(-8.5e-4 * sp.sqrt(1 + x))


def f_shield_H_WG11(N_H, alpha=0.149, delta=1.62, xi=2.85e23):
    """Wolcott-Green & Haiman 2011: Lyman-series shielding of the LW band by atomic H, floored at 1e-15"""
    return sp.Max((1 + N_H / xi) ** -delta * sp.exp(-alpha * N_H / xi), 1e-15)


def cosmic_ray_heat_per_H_ionization(x_e_tot):
    """Bialy & Sternberg 2019 Eq. 39, erg per primary ionization; x_e_tot is n_e over all particles"""
    return 6.43 * (1 + 4.06 * sp.sqrt(x_e_tot / (x_e_tot + 0.07))) * EV


def cosmic_ray_heat_per_H2_ionization(n_H2):
    """KROME's fit (Glassgold et al. 2012 data, arXiv:1502.03380) of the heat per H2 ionization, erg, in the H2 number
    density [cm^-3]"""
    lg = sp.log(sp.Max(n_H2, 1e-40), 10)
    q = sp.Piecewise((0, n_H2 < 1.0366017494884676e-05),
                     (1.4510849135088 * lg + 7.23277032061136, n_H2 < 1e5),
                     (9.1462198642562 - 0.443120145068494 * lg + 0.601271617560071 * lg**2
                      - 0.0710284101055109 * lg**3 + 0.00242079494592405 * lg**4, n_H2 <= 1e10),
                     (18.021719523307954, True))
    return q * EV


def table_data():
    """the H2 escape probability (Omukai) and Semenov et al. 2003 Planck-mean tables of the KROME build"""
    from importlib.resources import files
    with files("jaco.models.starforge_menon").joinpath("data/krome_tables.npz").open("rb") as f:
        d = np.load(f)
        return {k: np.array(d[k]) for k in d.files}
