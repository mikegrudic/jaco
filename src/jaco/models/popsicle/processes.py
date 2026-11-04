"""The processes of the popsicle model that differ from starforge's: Menon et al. 2025's KROME network and
thermal terms as their collapse test runs them (see __init__ for the inventory)."""

import numpy as np
import sympy as sp
from jaco.processes import Reaction, ThermalTerm, Transfer
from jaco.symbols import n_, x_, table_interp_2d
from ..starforge.symbols import T, log_T, n_Htot, Z_dust, f_dust, NH, z
from ..starforge.h2_chemistry.grain_formation import H2_dust_formation_rate
from . import rates as r
from .rates import EV, h2heatfac

KROME = "2026arXiv261017577M"  # Menon et al. 2025 (KROME POPSICLE network)
SIGMA_SB = 5.670374e-5

chi_0 = sp.Symbol("chi_0")  # field strength in Draine units
Gamma_ion0 = sp.Symbol("Gamma_ion0")  # the driver's unattenuated H photoionization rate [s^-1]
zeta = sp.Symbol("zeta_CR")  # KROME's user_crate, the total CR ionization rate of H2 [s^-1]
lambda_J = sp.Symbol("lambda_J")  # Jeans length sqrt(pi k T / (G rho mu m_p)) [cm]
rho = sp.Symbol("rho")
x_D_tot = sp.Symbol("x_D,tot")
Td = sp.Symbol("Td")
T_cmb = 2.73 * (1 + z)

# A_V of the shielding column N_H = n_Htot lambda_J / 2 (KROME's num2col), with dust-to-gas Z_d
A_V = Z_dust * NH / 1.87e21
chi_LW = chi_0 * sp.exp(-2.805 * A_V)
chi_PE = chi_0 * sp.exp(-1.608 * A_V)
chi_FUV = (1.8e-4 * chi_PE + 3e-5 * chi_LW) / 2.1e-4
G_0_habing = 1.69 * chi_FUV

# all particles (KROME's ntot, which the CR heating's electron fraction and the H2 escape column use)
n_particles = n_Htot * sum(x_(s) for s in ("H", "H_2", "H+", "He", "He+", "He++", "e-"))


def photochemistry():
    return [
        # The driver's H-ionizing field: 2.19e-12 chi s^-1 attenuated by exp(-1.14e4 A_V) only, A_V proportional to
        # Z, so unshielded at Z = 0; 3.4 eV per photoionization (heat_photoAv)
        Reaction("H -> H+ + e-", Gamma_ion0 * chi_0 * sp.exp(-1.14e4 * A_V), heat_per_reaction=3.4 * EV,
                 name="Photoionization of H by the driver's ionizing field", bibliography=[KROME]),
        Reaction("H_2 -> H + H",
                 5.6e-11 * chi_LW * r.f_shield_H2_WG11(x_("H_2") * NH) * r.f_shield_H_WG11(x_("H") * NH),
                 heat_per_reaction=(0.4 + 16 * h2heatfac()) * EV,
                 name="Photodissociation of H_2 (WG11 shielding)",
                 bibliography=[KROME, "2011MNRAS.418..838W", "2011MNRAS.412.2603W", "1977ApJS...34..405B"]),
        Reaction("H- -> H + e-", 7.1e-7 * chi_0 * sp.exp(-0.5 * A_V), name="Photodetachment of H-",
                 bibliography=["2007ApJ...666....1G"]),
    ]


def cosmic_rays():
    x_e_tot = n_("e-") / n_particles
    Q_H = r.cosmic_ray_heat_per_H_ionization(x_e_tot)
    return [
        Reaction("H -> H+ + e-", (0.5 + 3.731 * x_("H_2")) * zeta, heat_per_reaction=Q_H,
                 name="Cosmic-ray ionization of H", bibliography=[KROME, "2019ApJ...878..120B"]),
        Reaction("He -> He+ + e-", 0.55 * zeta, heat_per_reaction=Q_H, name="Cosmic-ray ionization of He",
                 bibliography=[KROME]),
        Reaction("H_2 -> H + H", 0.7 * zeta, heat_per_reaction=20 * EV, name="Cosmic-ray dissociation of H_2",
                 bibliography=[KROME]),
        Reaction("H_2 -> H + H+ + e-", 0.023 * zeta, heat_per_reaction=20 * EV,
                 name="Cosmic-ray dissociative ionization of H_2", bibliography=[KROME]),
        Reaction("H- -> H + e-", 924 * x_("H_2") * zeta, heat_per_reaction=20 * EV,
                 name="Cosmic-ray detachment of H-", bibliography=[KROME]),
        # H2 -> H2+ + e-: the electrons are the molecular ions of starforge's electron budget; here only its heat
        ThermalTerm(0.976 * zeta * n_("H_2") * r.cosmic_ray_heat_per_H2_ionization(n_("H_2")),
                    name="Cosmic-ray heating by H_2 ionization", bibliography=[KROME, "2012ApJ...756..157G"]),
    ]


def photoelectric_heating():
    """Weingartner & Draine 2001 Eq. 44 (R_V = 3.1, b_C = 4e-5, ISRF), psi + 50 as Kim et al. 2023"""
    psi = G_0_habing * sp.sqrt(T) / sp.Max(n_("e-"), 1e-40) + 50.0
    eps = (5.22 + 2.25 * T**0.147) / (1 + 0.04996 * psi**0.431 * (1 + 0.00430 * psi**0.692))
    return ThermalTerm(1e-26 * G_0_habing * n_Htot * eps * Z_dust, name="Photoelectric heating (WD01)",
                       bibliography=["2001ApJS..134..263W", "2023ApJS..264...10K"])


def h2_network():
    """H-, H2+, H3+ and H2 reactions of the network with their chemical heat (KROME's heatingChem). H2+ is not a
    species: its four fates (H2 + H+ with H, H3+ + H with H2, 2H with e-, H + H+ by photodissociation) split each
    H2+-forming reaction into effective reactions in proportion to their rates. H3+ is a species; it recombines."""
    n_e_ions = n_Htot * (x_("H+") + x_("He+") + 2 * x_("He++") + x_("H_3+"))
    gamma_H2p = 1.1e-9 * chi_0 * sp.exp(-1.9 * A_V)
    d = r.k13 * n_("H") + r.k23 * n_e_ions + gamma_H2p + r.k153 * n_("H_2")
    f_H, f_e, f_ph = r.k13 * n_("H") / d, r.k23 * n_e_ions / d, gamma_H2p / d
    f_3 = r.k153 * n_("H_2") / d
    heat13 = 1.83 * EV * h2heatfac()
    form = lambda q: q * h2heatfac() * EV
    return [
        Reaction("H + e- -> H- + photon_assoc,H", r.k9, name="Radiative association of H with e-",
                 bibliography=["1972ApL....11..165D"]),
        Reaction("H + H- -> H_2 + e-", r.k10, heat_per_reaction=form(3.53),
                 name="Associative detachment of H with H-", bibliography=["2010Sci...329...69K"]),
        Reaction("H- + H+ -> H + H", r.k21, name="Mutual neutralization of H+ and H-",
                 bibliography=["2009PhRvA..79a2713S"]),
        Reaction("H- + H+ + H -> H_2 + H+ + e-", rate=r.k22 * n_("H-") * n_("H+") * f_H, heat_per_reaction=heat13,
                 name="H- + H+ -> H2+ + e-, then H2+ + H -> H2 + H+", bibliography=["1978JPhB...11L.671P"]),
        Reaction("H- + H+ -> H + H", rate=r.k22 * n_("H-") * n_("H+") * f_e,
                 name="H- + H+ -> H2+ + e-, then H2+ + e- -> 2H", bibliography=["1978JPhB...11L.671P"]),
        Reaction("H- -> H + e-", rate=r.k22 * n_("H-") * n_("H+") * f_ph,
                 name="H- + H+ -> H2+ + e-, then H2+ photodissociation", bibliography=["1978JPhB...11L.671P"]),
        Reaction("H + H+ + H -> H_2 + H+", rate=r.k11 * n_("H") * n_("H+") * f_H, heat_per_reaction=heat13,
                 name="H + H+ -> H2+, then H2+ + H -> H2 + H+", bibliography=["1976JPhB....9L.239R"]),
        Reaction("H_2 + H+ + e- -> H + H + H", rate=r.k14 * n_("H_2") * n_("H+") * f_e,
                 name="H2 + H+ -> H2+ + H, then H2+ + e- -> 2H", bibliography=["2004ApJ...606L.167S"]),
        Reaction("H_2 + H+ -> H + H + H+", rate=r.k14 * n_("H_2") * n_("H+") * f_ph,
                 name="H2 + H+ -> H2+ + H, then H2+ photodissociation", bibliography=["2004ApJ...606L.167S"]),
        ThermalTerm(r.k14 * n_("H_2") * n_("H+") * f_H * heat13, name="Heat of H2+ + H -> H2 + H+ after H2 + H+",
                    bibliography=[KROME]),
        Reaction("2H_2 + H+ -> H_3+ + 2H", rate=r.k14 * n_("H_2") * n_("H+") * f_3,
                 name="H2 + H+ -> H2+ + H, then H2+ + H2 -> H3+ + H", bibliography=["2004ApJ...606L.167S"]),
        Reaction("H+ + H_2 -> H_3+", rate=r.k11 * n_("H") * n_("H+") * f_3,
                 name="H + H+ -> H2+, then H2+ + H2 -> H3+ + H", bibliography=["1976JPhB....9L.239R"]),
        Reaction("H- + H+ + H_2 -> H_3+ + H + e-", rate=r.k22 * n_("H-") * n_("H+") * f_3,
                 name="H- + H+ -> H2+ + e-, then H2+ + H2 -> H3+ + H", bibliography=["1978JPhB...11L.671P"]),
        Reaction("H_3+ + e- -> H_2 + H", r.k162, name="Dissociative recombination of H3+ to H2 + H",
                 bibliography=["2010MNRAS.404....2G"]),
        Reaction("H_3+ + e- -> 3H", r.k163, name="Dissociative recombination of H3+ to 3H",
                 bibliography=["2010MNRAS.404....2G"]),
        Reaction("H_2 + e- -> H + H-", r.k15, name="Dissociative attachment of H_2",
                 bibliography=["2007JPhB...40..727C"]),
        Reaction("H_2 + e- -> H + H + e-", r.k16, heat_per_reaction=-4.48 * EV,
                 name="Collisional dissociation of H_2 by e-", bibliography=["1983ppam.book..237M"]),
        Reaction("H_2 + H -> H + H + H", r.k17_hunter(), heat_per_reaction=-4.48 * EV,
                 name="Collisional dissociation of H_2 by H",
                 bibliography=["2023MNRAS.519.4152H", "1996ApJ...461..265M"]),
        Reaction("H_2 + H_2 -> H + H + H_2", r.k27_hunter(), heat_per_reaction=-4.48 * EV,
                 name="Collisional dissociation of H_2 by H_2",
                 bibliography=["2023MNRAS.519.4152H", "1996ApJ...461..265M"]),
        Reaction("H + H + H -> H_2 + H", r.k25, heat_per_reaction=form(4.48), name="3-body formation of H_2 (H)",
                 bibliography=["2013ApJ...773L..25F"]),
        Reaction("H + H + H_2 -> H_2 + H_2", r.k26, heat_per_reaction=form(4.48), name="3-body formation of H_2 (H_2)",
                 bibliography=["2008MNRAS.388.1627G"]),
        Reaction("H + H -> H_2", rate=H2_dust_formation_rate * n_Htot * n_("H"), clumping=1,
                 heat_per_reaction=(0.2 + 4.2 * h2heatfac()) * EV, name="Formation of H_2 on dust grains",
                 bibliography=["1979ApJS...41..555H"]),
        Reaction("H+ + e- -> H", r.k2, heat_rate_coefficient=-r.lambda_rec_Hplus, name="Gas-phase recombination of H+",
                 bibliography=["1997NewA....2..181A", "1992ApJS...78..341C"]),
    ]


def hd_abundance():
    """HD per H nucleus with the deuterium in fractionation equilibrium between D and HD through D+ (D + H+ <-> D+ + H,
    D+ + H2 <-> HD + H+): HD/D = (x_H2/x_H) (k33 k38)/(k34 k40), total x_D,tot. The KROME network follows D, D+, HD
    in time; this is its equilibrium."""
    R = x_("H_2") / (x_("H") + 1e-30) * (r.k33 * r.k38) / (r.k34 * r.k40)
    return x_D_tot * R / (1 + R)


def _sigmoid(x, x0, s):
    return 10 / (10 + sp.exp(sp.Min(-s * (x - x0), 600)))


def _wcool(lt, lo, hi):
    """KROME's smooth window for the H2 cooling fits"""
    x = (lt - lo) / (hi - lo)
    return 10 ** (200 * (_sigmoid(x, -0.2, 50) * _sigmoid(-x, -1.2, 50) - 1))


def h2_cooling():
    """KROME's cooling_H2: Glover & Abel 2008 low-density rates with Glover 2015's updates (o:p = 3), the HM79 / Glover
    2015 LTE limit, combined harmonically and times Omukai's escape probability on the column of all particles"""
    lt = log_T
    t3 = T / 1e3
    l3 = sp.log(t3, 10)
    w14, w24 = _wcool(lt, 1, 4), _wcool(lt, 2, 4)
    pw = lambda c: sum(ci * l3**i for i, ci in enumerate(c))
    f_H = sp.Piecewise(
        (10 ** pw([-16.818342, 37.383713, 58.145166, 48.656103, 20.159831, 3.8479610]), T <= 100),
        (10 ** pw([-24.311209, 3.5692468, -11.332860, -27.850082, -21.328264, -4.2519023]), T <= 1e3),
        (10 ** pw([-24.311209, 4.6450521, -3.7209846, 5.9369081, -5.5108049, 1.5538288]), T <= 6e3),
        (1.862314467912518e-22 * _wcool(lt, 1, np.log10(6e3)), True)) * n_("H")
    f_Hp = sp.Piecewise(
        (10 ** pw([-22.089523, 1.5714711, 0.015391166, -0.23619985, -0.51002221, 0.32168730]), (T > 10) & (T <= 1e4)),
        (1.182509139382060e-21 * w14, True)) * n_("H+")
    f_H2 = w24 * 10 ** pw([-23.962112, 2.09433740, -0.77151436, 0.43693353, -0.14913216, -0.033638326]) * n_("H_2")
    f_e = w24 * sp.Piecewise(
        (10 ** sp.Min(pw([-21.928796, 16.815730, 96.743155, 343.19180, 734.71651, 983.67576, 801.81247, 364.14446,
                          70.609154]), 30), T <= 500),
        (10 ** pw([-22.921189, 1.6802758, 0.93310622, 4.0406627, -4.7274036, -8.8077017, 8.9167183, 6.4380698,
                   -6.3701156]), True)) * n_("e-")
    f_He = sp.Piecewise(
        (10 ** pw([-23.689237, 2.1892372, -0.81520438, 0.29036281, -0.16596184, 0.19191375]), (T > 10) & (T <= 1e4)),
        (1.002560385050777e-22 * w14, True)) * n_("He")
    LDL = f_H + f_Hp + f_H2 + f_e + f_He
    HDL = sp.Piecewise(
        ((9.5e-22 * t3**3.76) / (1 + 0.12 * t3**2.1) * sp.exp(-((0.13 / t3) ** 3)) + 3e-24 * sp.exp(-0.51 / t3)
         + 6.7e-19 * sp.exp(-5.86 / t3) + 1.6e-18 * sp.exp(-11.7 / t3), T < 2e3),
        (10 ** pw([-20.584225, 5.0194035, -1.5738805, -4.7155769, 2.4714161, 5.4710750, -3.9467356, -2.2148338,
                   1.8161874]), T <= 1e4),
        (5.531333679406485e-19 / (1 + sp.exp(sp.Min((T - 3e4) * 2e-4, 300))), True))
    tab = r.table_data()
    axes = [np.logspace(20, 29, 10), np.logspace(np.log10(200), np.log10(5000), 20)]  # log N, log T: uniform
    log_esc = table_interp_2d("h2_escape_omukai", tab["h2esc_log10"], axes,
                              sp.Max(n_particles * lambda_J / 2, 1e10), T, log_axes=[True, True])
    return ThermalTerm(-n_("H_2") * LDL * HDL / (LDL + HDL) * 10**log_esc, name="H2 line cooling (KROME)",
                       bibliography=["2008MNRAS.388.1627G", "2015MNRAS.451.2082G", "1979ApJS...41..555H",
                                     "2010ApJ...722.1793O"])


_LIPOVKA = [[-42.56788, 0.92433, 0.54962, -0.07676, 0.00275],
            [21.93385, 0.77952, -1.06447, 0.11864, -0.00366],
            [-10.19097, -0.54263, 0.62343, -0.07366, 0.002514],
            [2.19906, 0.11711, -0.13768, 0.01759, -0.00066631],
            [-0.17334, -0.00835, 0.0106, -0.001482, 0.00006192]]


def hd_cooling():
    """Lipovka, Nunez-Lopez & Avila-Reese 2005, with KROME's density argument (atomic H only, clamped to 1e-2..1e10)
    and T clamped at 1e4 K"""
    lT = sp.log(sp.Min(T, 1e4), 10)
    ln = sp.log(sp.Min(sp.Max(n_("H"), 1e-2), 1e10), 10)
    logW = sum(_LIPOVKA[i][j] * lT**i * ln**j for i in range(5) for j in range(5))
    return ThermalTerm(-10**logW * n_("HD") * sp.Piecewise((1, T > T_cmb), (0, True)), name="HD line cooling (Lipovka)",
                       bibliography=["2005MNRAS.361..850L"])


def continuum_cooling():
    """Omukai 2000's gas continuum (Lenzuni et al. 1996 Planck mean, zero below 1e-12 g cm^-3) with escape
    probability min(1, tau^-2) over the Jeans length"""
    kappa = sp.Piecewise((10 ** (1.000042 * sp.log(sp.Min(rho, 0.5), 10) + 2.14989), rho > 1e-12), (0, True))
    tau = lambda_J * kappa * rho
    return ThermalTerm(-4 * SIGMA_SB * T**4 * kappa * rho / sp.Max(1, tau**2), name="Gas continuum cooling",
                       bibliography=["2000ApJ...534..809O", "1996ApJS..106..273L"])


def dust_processes():
    """KROME's compute_Semenov_Tdust: dust emission 4 sigma rho kappa_P(rho, Td) Z_d beta Td^4 against the CMB, the
    attenuated FUV (chi_FUV 2.1e-4 erg cm^-2 s^-1 sr^-1 on 1e-21 cm^2 per H times Z_d) and the gas (alpha_gd,
    3.2e-34 erg cm^3 s^-1 K^-3/2 in KROME, Krumholz 2014); beta = 1 below tau = 1 and tau^-2 above (here
    (1 + tau^4)^-1/2), tau over the Jeans length"""
    tab = r.table_data()
    log_rho = np.linspace(tab["semenov_log_rho"][0], tab["semenov_log_rho"][-1], len(tab["semenov_log_rho"]))
    Td_axis = np.geomspace(tab["semenov_Td"][0], tab["semenov_Td"][-1], len(tab["semenov_Td"]))
    # Where an evaporating component makes kappa_P Td^4 fall with Td (ice at ~140 K, organics at ~420 K) the dust
    # balance has several roots, which KROME's Newton from the previous Td picks between; here the emission is held at
    # its running maximum in Td so the root is unique
    kappa_T4 = np.maximum.accumulate(tab["semenov_kappaP"] * Td_axis**4, axis=1)
    log_kappa = table_interp_2d("semenov_planck_krome", np.log10(kappa_T4 / Td_axis**4), [10**log_rho, Td_axis],
                                sp.Min(sp.Max(rho, 10 ** log_rho[0]), 10 ** log_rho[-1]),
                                sp.Min(sp.Max(Td, Td_axis[0]), Td_axis[-1]), log_axes=[True, True])
    kappa = 10**log_kappa * Z_dust * f_dust
    tau = rho * kappa * lambda_J
    emissivity = 4 * SIGMA_SB * rho * kappa / sp.sqrt(1 + tau**4)
    uv = chi_FUV * 2.1e-4 * 4 * np.pi * 1e-21 * Z_dust * f_dust * n_Htot
    return [
        Transfer(emissivity * Td**4, {"dust heat": -1}, name="Dust thermal emission",
                 bibliography=["2003A&A...410..611S", KROME]),
        Transfer(emissivity * T_cmb**4 + uv, {"dust heat": 1}, name="Dust absorption of the CMB and the FUV field",
                 bibliography=[KROME]),
        ThermalTerm(sp.Symbol("alpha_gd") * sp.sqrt(T) * (Td - T) * Z_dust * f_dust * n_Htot**2,
                    name="Gas-dust collisions",
                    bibliography=["2014MNRAS.437.1662K"], clumping=1, reservoir="dust heat"),
    ]


def fine_structure_with_escape():
    """starforge's [CII], [CI] and [OI] lines times KROME's optical-depth factors (DESPOTIC, LTE, tabulated in T and
    the coolant's column x n_Htot lambda_J / 2; linear in the factor, clamped at the table edges)"""
    from ..starforge.fine_structure import FineStructureCooling, coolant_abundances
    tab = r.table_data()
    x = coolant_abundances()
    out = []
    for coolant, key in (("C+", "Cp"), ("C", "C"), ("O", "O")):
        lT, lN = tab[f"fs_{key}_logT"], tab[f"fs_{key}_logN"]
        lT_axis, lN_axis = np.linspace(lT[0], lT[-1], len(lT)), np.linspace(lN[0], lN[-1], len(lN))
        column = sp.Min(sp.Max(x[coolant] * n_Htot * lambda_J / 2, 10 ** lN_axis[0]), 10 ** lN_axis[-1])
        factor = table_interp_2d(f"fs_escape_{key}", tab[f"fs_{key}_factor"], [10**lT_axis, 10**lN_axis],
                                 sp.Min(sp.Max(T, 10 ** lT_axis[0]), 10 ** lT_axis[-1]), column, log_axes=[True, True])
        p = FineStructureCooling(coolant, x[coolant])
        out.append(ThermalTerm(p.network["heat"].rhs * factor, name=p.name + " (optical depth)",
                               bibliography=list(p.bibliography) + ["2014MNRAS.437.1662K", KROME]))
    return out
