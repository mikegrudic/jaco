"""GIZMO's IR-band opacities and dust survival fraction under RT_INFRARED, as expressions with derivatives in the dust,
radiation and gas temperatures.

References are to radiation/rt_dust_opacity.cc (dust_planck_mean_opacity), radiation/rt_utilities.cc
(rt_kappa_adaptive_IR_band, rt_kappa) and eos/eos.cc (return_dust_to_metals_ratio_vs_solar) at gizmo_jaco_dev 2df0c6dd.

The Semenov et al. (2003) table switches dust composition discontinuously at T_dust = 160, 275, 425 and 680 K (ice and
refractory sublimation lines). The Planck mean drops 2.4x at 160 K and 3.1x at 425 K, so the dust emission
4 sigma kappa(T_d, T_d) T_d^4 falls across those switches and, for heating rates between its values on either side, the
dust energy balance has two roots, one per composition; GIZMO's bracketed search (rt_eqm_dust_temp, rt_ir_lambdadust)
takes the one its walk from the previous dust temperature reaches first, a hysteresis. Each switch here is a C1
smoothstep in T_dust over [boundary - 5 K, boundary + 5 K], which keeps both roots (and a third, unstable one inside
the window) and gives Newton finite derivatives; the solver starts from the previous dust temperature as GIZMO's walk
does. Smoothing wide enough to make the emission monotonic (23.5 and 97.5 K half-widths at 160 and 425 K) would remove
the hysteresis and move the dust temperature by up to ~100 K near 425 K. Outside the windows the opacity is the table's
exactly.
"""

import numpy as np
import sympy as sp
from jaco.symbols import piecewise_linear, x_
from .symbols import T

T_dust = sp.Symbol("Td")
T_rad = sp.Symbol("T_rad")
rho = sp.Symbol("rho")  # gas mass density [g cm^-3]
Z_metals = sp.Symbol("Z_metals")  # metal mass fraction (GIZMO's Metallicity[0]; 0.014 without METALS)
gamma_eos = sp.Symbol("gamma_eos")  # adiabatic index of the cell at the start of the step
Z_dust = sp.Symbol("Z_d")  # metallicity in solar units (GIZMO's Zfac)
HYDROGEN_MASSFRAC = 0.76  # GIZMO's
PROTONMASS, BOLTZMANN = 1.6726e-24, 1.38066e-16  # GIZMO's PROTONMASS_CGS, BOLTZMANN_CGS

LOG_T_RAD = np.linspace(0.0, 4.0, 15)
LOG_KAPPA = np.array([
    [-1.909515215033498, -1.5017616295543856, -1.2610211141587906, -1.0643751130254193, -0.7661794028048912,
     -0.2485164276172981, 0.3936319485109052, 0.7185651396015793, 1.003536500936941, 1.0703750048744685,
     1.185414318657744, 1.4334392971521788, 1.6154100043688104, 1.833331149478953, 2.2402919406422592],
    [-2.0584711507666156, -1.639523808954471, -1.3937562378746238, -1.217712452271845, -1.0035307263948003,
     -0.6537979484512214, -0.14444336531252724, 0.33781494089162734, 0.6204502553977504, 0.7073792380378322,
     0.817126602196553, 1.1620592165098858, 1.4775329387335934, 1.7314724407802422, 2.122441393841703],
    [-2.092756656625387, -1.6688435106176933, -1.4197982288972344, -1.2491484838263678, -1.0526101054461765,
     -0.7239962293949143, -0.21897626943630635, 0.26791745490238816, 0.5498747399044497, 0.6380011931553505,
     0.7545581854919681, 1.1204828950445178, 1.447717175478172, 1.6975481925473666, 2.066531977172171],
    [-2.466846959771108, -1.9983094494283768, -1.7094791411944892, -1.5454814504892613, -1.4403397374848839,
     -1.2856048236415571, -0.8297121104893529, -0.28151883033918534, 0.004258175582704628, 0.12669094376131687,
     0.31518820522515606, 0.8267768233252955, 1.2296194436753016, 1.4691152621947785, 1.6952767195239868],
    [-3.7448342907548136, -3.298656729562377, -2.8766563714426696, -2.4999445785716605, -2.1658033063307927,
     -1.8743193706624732, -1.5588128589881036, -1.0606192296411485, -0.33965876540303136, 0.5368140402181846,
     1.1915541410857275, 1.4384491935887818, 1.4792273276721044, 1.5197919949174217, 1.6603749111299055],
])
ZONE_BOUNDARIES = (160.0, 275.0, 425.0, 680.0)  # GIZMO's Tdust_zones; above 1500 K RT_INFRARED keeps the last zone
ZONE_HALF_WIDTHS = (5.0, 5.0, 5.0, 5.0)  # K


def _zone_opacity(zone, T_r):
    """Planck mean of one composition zone at radiation temperature T_r [cm^2/g, solar]: linear in log T_rad and
    log kappa on the table, constant beyond it"""
    return 10 ** piecewise_linear(LOG_T_RAD, LOG_KAPPA[zone], sp.log(T_r, 10), name=f"kappa_dust_zone{zone}")


def _smoothstep(x):
    """0 below 0, 1 above 1, 3x^2 - 2x^3 between"""
    s = sp.Min(1, sp.Max(0, x))
    return s * s * (3 - 2 * s)


def semenov_planck_mean(T_r, T_d):
    """dust_planck_mean_opacity(T_rad, T_dust) [cm^2/g at solar metallicity], the zone switches smoothed"""
    kappa = _zone_opacity(0, T_r)
    for i, (b, h) in enumerate(zip(ZONE_BOUNDARIES, ZONE_HALF_WIDTHS)):
        kappa += _smoothstep((T_d - b + h) / (2 * h)) * (_zone_opacity(i + 1, T_r) - _zone_opacity(i, T_r))
    return kappa


def dust_survival(T_d):
    """return_dust_to_metals_ratio_vs_solar under RT_INFRARED: the dust surviving sublimation at 1500 K"""
    s = 9 * (1 - T_d / 1500.0)
    return sp.Max(0.5 * (1 + s / sp.sqrt(1 + s * s)) * sp.exp(-sp.Min(40, (T_d / 1500.0) ** 2 / 9)), 1e-25)


def ir_dust_opacity(T_d, T_r, absorption=True):
    """The dust part of rt_kappa_adaptive_IR_band [cm^2/g]: absorption (flag -1; emission, flag 1, is this at
    T_r = T_d) or extinction (flag 0), times the metallicity and the surviving dust"""
    kappa = semenov_planck_mean(T_r, T_d)
    if absorption:
        kappa *= 1 - 0.5 / (1 + 725.0**2 / (1 + T_r**2))
    return kappa * Z_dust * dust_survival(T_d)


def opacity_gas_temperature():
    """rt_kappa_adaptive_IR_band's estimate of the gas temperature from the cell's specific energy at mu = 0.59:
    1 + 0.59 (gamma - 1) (m_p / k_B) u, here at the start of the step"""
    return 1 + 0.59 * (gamma_eos - 1) * (PROTONMASS / BOLTZMANN) * sp.Symbol("u_initial")


def ir_gas_opacity(T_r, T_d, absorption=True):
    """The non-dust part of rt_kappa_adaptive_IR_band [cm^2/g] (flag -1: absorption, electron scattering excluded;
    flag 0 with it): free-free/bound-free (Kramers), iron line blanketing, molecular lines, H- and Rayleigh scattering,
    with GIZMO's composition inputs: x_e (Ne), x_H+ (HII), x_H0 = 1 - x_H+ (HI), the molecular mass fraction 2 x_H2,
    and the metals not in dust at the dust temperature T_d"""
    X = HYDROGEN_MASSFRAC
    x_e, x_Hp = x_("e-"), x_("H+")
    x_H0 = 1 - x_Hp
    T_g = opacity_gas_temperature()
    f_neutral = sp.Max(0, 1 - x_e)
    f_free_metals = Z_metals * sp.Max(0, 1 - 0.5 * dust_survival(T_d))
    k_electron = 0.4 * X * x_e / ((1 + 2.7e11 * rho / T_g**2) * (1 + (T_r / 4.5e8) ** 0.86))
    k_molecular = 0.1 * (f_free_metals + 3e-9) * f_neutral * 2 * x_("H_2")
    k_kramers = (4.0e25 * (1 + X) * (f_free_metals * sp.exp(-sp.Min(1.5e5 / T_r, 40)) + 0.001 * x_e) * rho
                 / (T_r**3 * sp.sqrt(T_g)))
    k_kramers += (1.5e20 * f_free_metals * rho / T_r**2 * sp.exp(-sp.Min((0.8e4 / T_r) ** 4, 40))
                  * sp.exp(-sp.Min((T_r / 0.7e6) ** 2, 40)))
    k_rayleigh = f_neutral * sp.Min(5e-19 * T_r**4, 0.2 * (1 + X))
    tg = (T_g / 1.3e4) ** 2
    x_Hminus = 4e-10 * T_g * x_e * x_H0 / ((1 + x_Hp * 300 + x_e * 1000 * tg / (1 + tg) + 4e-17) * (1 + T_g / 3e4))
    k_bf = 4.2e7 * (8760 / T_r) ** 1.5 * sp.exp(-sp.Min(8760 / T_r, 40))
    phi = sp.Min(T_g / 5040, 2)
    k_ff = 1.9e6 * (8760 / T_r) ** 2 * sp.exp(-sp.Min(8760 / T_r, 40)) * (0.6 - 2.5 * sp.sqrt(phi) + 2.5 * phi
                                                                         + 2.7 * phi * sp.sqrt(phi))
    kappa = k_molecular + k_kramers + x_Hminus * (k_bf + k_ff) + k_rayleigh
    return kappa if absorption else kappa + k_electron


def band_dust_opacity(kappa_dust, floor_metals=0.0):
    """rt_kappa for a dust-absorbed band [cm^2/g]: max(kappa_HHe, kappa_dust Z f_dust), kappa_HHe = 0.02 + 0.35 x_e
    (neutral and free-electron floor); the photoelectric band floors Z f_dust at floor_metals"""
    Zf = Z_dust * dust_survival(T_dust)
    if floor_metals:
        Zf = sp.Max(floor_metals, Zf)
    return sp.Max(0.02 + 0.35 * x_("e-"), kappa_dust * Zf)
