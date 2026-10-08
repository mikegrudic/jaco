"""GIZMO's matter-radiation coupling as processes on solved radiation bands (M1 RADTRANSFER with RT_CHEM_PHOTOION,
RT_PHOTOELECTRIC, RT_NUV, RT_OPTICAL_NIR and RT_INFRARED).

Every term GIZMO's RT kick (rt_update_driftkick: absorption, donation, dust re-emission), its dust temperature
(rt_eqm_dust_temp, rt_ir_lambdadust) and its cooling module (photoionization and photoheating, the cooling-radiation
return to the bands) apply is a process here; the host transports the bands only. Line numbers refer to
radiation/rt_utilities.cc, cooling/cooling.cc, radiation/rt_dust_opacity.cc and eos/eos.cc at gizmo_jaco_dev 2df0c6dd.

Bands and their normalization (abundances per H nucleus, n_photon_X = n_Htot x_photon_X):

- photon_EUV: ionizing photons per H nucleus (13.6-500 eV; GIZMO counts them at rt_nu_eff_eV, the Param hnu_EUV);
- photon_FUV (GIZMO's PHOTOELECTRIC band, 8-13.6 eV), photon_NUV (3.444-8 eV), photon_ONIR (0.4133-3.444 eV), photon_IR
  (0.001-0.4133 eV, with the radiation temperature T_rad GIZMO tracks per cell): energy per H nucleus in eV.

Speeds of light: matter rates run at the true c, as GIZMO's cooling module and dust balance do; the bands are
transported at c_tilde = rsol c, so every band row is rsol times the matter-side rate (Reaction row_factors,
Transfer factors). Both are physical: GIZMO renormalizes the sources so that a band's energy density is the physical
one.
"""

import numpy as np
import sympy as sp
from jaco.processes import Reaction, ThermalTerm, Transfer
from jaco.symbols import n_, dt
from .symbols import T, n_Htot, nH_gizmo_cooling, sqrt_T
from .dust_opacity import (T_dust, T_rad, rho, Z_dust, ir_dust_opacity, ir_gas_opacity, band_dust_opacity,
                           dust_survival)

C_LIGHT = 2.9979e10  # GIZMO's C_LIGHT_CGS
EV = 1.60217733e-12  # GIZMO's ELECTRONVOLT_IN_ERGS: the energy unit of the bands counted in eV
SIGMA_SB = 5.67e-5  # as GIZMO writes the dust emission
U_HABING = 1.6e-3 / C_LIGHT  # Habing energy density [erg cm^-3]: HABING_FLUX_CGS / C_LIGHT_CGS

rsol = sp.Symbol("rsol")  # c_tilde / c
sigma_HI = sp.Symbol("sigma_HI")
eps_HI = sp.Symbol("eps_HI")
hnu_EUV = sp.Symbol("hnu_EUV")
T_CMB = sp.Symbol("T_CMB")
G_LW_bg = sp.Symbol("G_LW_bg")
gamma_12_UVB = sp.Symbol("gamma_12_UVB")
eps_H0_UVB = sp.Symbol("eps_H0_UVB")

EUV, FUV, NUV, ONIR, IR = "photon_EUV", "photon_FUV", "photon_NUV", "photon_ONIR", "photon_IR"
BANDS = (EUV, FUV, NUV, ONIR, IR)
PHOTON_DENSITIES = frozenset(n_(b) for b in BANDS)
EUV_SINKS = ("Photoionization of H by the ionizing band",)


def band_energy_eV(band):
    """Energy density of a band [eV cm^-3]"""
    return n_(band) * hnu_EUV if band == EUV else n_(band)


def blackbody_fraction(E_lower, E_upper, T_eff):
    """blackbody_lum_frac (rt_utilities.cc 1593): the fraction of a blackbody's energy between E_lower and E_upper [eV],
    with GIZMO's fits and its fall-backs where their difference rounds to <= 0 (both ends on the large-x fit, whose
    exponential GIZMO caps at e^-40)"""
    k_B = 8.617e-5
    x1, x2 = E_lower / (k_B * T_eff), E_upper / (k_B * T_eff)

    def large_x(x):  # 1 - cumulative fraction on the large-x fit
        return 0.15398973382026504 * (6 + x * (6 + x * (3 + x))) * sp.exp(-sp.Min(x, 40))

    def cumulative(x):
        return sp.Piecewise((131.4045728599595 * x**3 / (2560 + x * (960 + x * (232 + 39 * x))), x < 3.40309),
                            (1 - large_x(x), True))
    df_large = large_x(x1) - large_x(x2)  # both ends at x >= 3.40309
    tail = 0.15398973382026504 * (6 + x1 * (6 + x1 * (3 + x1))) * sp.exp(-sp.Min(x1, 120))
    return sp.Piecewise((sp.Max(cumulative(x2) - cumulative(x1), 0), x1 < 3.40309),
                        (tail, (df_large <= 0) & (x1 > 4) & (x1 < 120)),
                        (2e-47, (df_large <= 0) & (x1 >= 120)),
                        (sp.Max(df_large, 0), True))


# --- photoionization ---------------------------------------------------------------------------------------------

def photoionization(donation=None):
    """H + photon_EUV -> H+ + e-: photoionizations c sigma_HI n_gamma n_H0 per unit volume at the true c
    (find_abundances_and_rates, cooling.cc 846-904), heat eps_HI = rt_ion_G_HI each (Heat_Ion_from_RHD, 1225-1266),
    one photon of the band taken per event at c_tilde (the kick's absorption, kappa rho = sigma_HI n_H0 at
    rt_utilities.cc 141-151, with its albedo 1e-6 neglected). With donation (a band counted in eV), the photon's energy
    hnu_EUV goes to it, as the kick donates the absorbed ionizing energy to the optical band (rt_get_donation_target_bin,
    1208-1212). He is not photoionized without RT_CHEM_PHOTOION_HE.

    The rate acts on atomic H, as the model's ionization balance does. GIZMO's neutral H, HI = 1 - HII, includes the nuclei
    in H2: its kick absorbs the band on them, its H balance photoionizes them, and its H2 then stays a fixed fraction of
    the neutral H, so H2 survives in gas in photoionization equilibrium. A channel H_2 + 2 photon_EUV -> 2 H+ + 2 e- at the
    same rate per nucleus absorbs as GIZMO does but destroys H2 at the gross rate (recombinations make atomic H): on
    HII_region_simple it removed the H2 coolant from the partially ionized front and moved r_IF 3-8% beyond GIZMO's,
    against 0.99-1.02 of it without the channel. Molecular gas therefore absorbs the band less than in GIZMO, by
    x_H / (x_H + 2 x_H2).

    With the band's photons solved in the same backward-Euler step, a cell takes at most the photons it has, so GIZMO's
    rate law needs no photon-limited time average (slab_averaging_function, commented out at cooling.cc 875 and 1244):
    on the Iliev test 1 R-type front at c_tilde = 0.1 c, with tau ~ 40 per cell and step, r_I/analytic is 0.976-1.020
    over 0.5-4 t_rec either way."""
    rows = {EUV: rsol}
    equation = "H + photon_EUV -> H+ + e-"
    if donation:
        equation += f" + {donation}"
        rows[donation] = rsol * hnu_EUV
    bib = ["GIZMO cooling.cc find_abundances_and_rates and Heat_Ion_from_RHD", "GIZMO rt_utilities.cc rt_kappa",
           "GIZMO rt_chem.cc rt_get_sigma"]
    name = "Photoionization of H by the ionizing band"
    return Reaction(equation, C_LIGHT * sigma_HI, eps_HI, clumping=1, row_factors=rows, name=name, bibliography=bib)


def ir_tail_photoionization():
    """Photoionization of H by the IR band's blackbody tail above 13.6 eV at T_rad (cooling.cc 1235-1238 and 859-878,
    rt_irband_egydensity_in_band), photons counted at max(hnu_EUV, T_rad/2959.81) eV; GIZMO's kick does not take them
    from the IR band, so neither does this. Atomic H only: the tail matters only where T_rad >~ 2e4 K"""
    n_tail = n_(IR) * blackbody_fraction(13.6, 500.0, T_rad) / sp.Max(hnu_EUV, T_rad / 2959.81)
    return Reaction("H -> H+ + e-", rate=C_LIGHT * sigma_HI * n_("H") * n_tail, heat_per_reaction=eps_HI, clumping=1,
                    name="Photoionization of H by the IR band's tail",
                    bibliography=["GIZMO cooling.cc CoolingRate and find_abundances_and_rates"])


# --- dust ----------------------------------------------------------------------------------------------------------

DUST_BAND_OPACITY = {FUV: (720.0, 1e-4), NUV: (480.0, 0.0), ONIR: (180.0, 0.0)}  # rt_kappa 209-221: kappa, Z floor


KICK_EXPONENT_CAP = 10.0


def kick_absorption_factor(a_dt):
    """expm1(a dt) / (a dt): the factor on an absorption rate a for which the backward-Euler step gives the kick's
    exponential, 1 / (1 + factor a dt) = exp(-a dt), for a dt up to KICK_EXPONENT_CAP: beyond it the band keeps e^-10
    of its energy over the step rather than e^-(a dt). GIZMO's kick caps the exponent at 50; a factor of e^50 on the
    rows of a band absorbed in a dense cell made a few percent of Newton's tier-1 attempts fail their line search
    (shu_M120: 2.7% of the solves to tier 2 or 3, a mean 23 evaluations per solve, against 0.08% and 3 at 10)"""
    x = sp.Min(a_dt, KICK_EXPONENT_CAP)
    return sp.Piecewise((1 + x * (sp.Rational(1, 2) + x / 6), x < 1e-3), ((sp.exp(x) - 1) / a_dt, True))


def dust_band_rate(band):
    """c f_abs kappa rho [s^-1]: the dust's absorption rate of a non-ionizing band at the true c, with rt_kappa's
    opacity (max of the neutral/electron floor and the dust's) and albedo f_abs = 1/2 (rt_absorb_frac_albedo), the
    surviving dust at the start-of-step dust temperature, as GIZMO's kick takes rt_kappa from the cell before the dust
    temperature is updated. The band's absorption is then linear in its energy at fixed rate over the solve: through
    the exponential kick factor, a rate depending on the solved Td makes the rows of a band absorbed in a dense cell
    exponentially sensitive to Td, and Newton's line search fails there"""
    kappa, floor = DUST_BAND_OPACITY[band]
    return C_LIGHT * 0.5 * band_dust_opacity(kappa, floor, T_d=T_dust_initial) * rho


def kick_factor(band):
    """The intermediate holding a dust-absorbed band's kick_absorption_factor (kick_absorption_intermediates)"""
    return sp.Symbol(f"kick_{band}")


T_dust_initial = sp.Symbol("Td_initial")  # dust temperature at the start of the step [K]


def kick_absorption_intermediates(bands=(FUV, NUV, ONIR)):
    """(Symbol, expression) per dust-absorbed band: kick_absorption_factor at the band's absorption over the step,
    a dt = c_tilde f_abs kappa rho dt, with the dust's opacity only (an intermediate cannot read the electron
    abundance, which rt_kappa's neutral/electron floor does; the floor exceeds the dust's opacity only where the dust
    is gone or the metallicity is below ~1e-3, and there the absorption falls back toward backward Euler), at the
    start-of-step dust temperature (dust_band_rate)"""
    out = []
    for band in bands:
        kappa, floor = DUST_BAND_OPACITY[band]
        Zf = Z_dust * dust_survival(T_dust_initial)
        a = C_LIGHT * 0.5 * kappa * (sp.Max(floor, Zf) if floor else Zf) * rho
        out.append((kick_factor(band), kick_absorption_factor(rsol * a * dt)))
    return out


def band_energy_over_step(band):
    """A dust-absorbed band's energy density averaged over the step [eV cm^-3], under the kick's exponential absorption
    (kick_absorption_factor times the energy at the step's end): what the dust absorbs and the gas sees. The RT kicks'
    transport brackets the absorption here, so the band's energy at the step's end has been absorbed over the whole
    step; GIZMO's cooling sees the band mid-step (after the first half-kick), as the average does"""
    return kick_factor(band) * band_energy_eV(band)


def dust_band_power(band):
    """What the dust absorbs from a non-ionizing band, per unit volume per unit time [erg cm^-3 s^-1]: the rate on
    the band's energy over the step (band_energy_over_step), so that the band ends the step at the kick's
    exp(-a dt) of its initial energy, a = c_tilde f_abs kappa rho (the backward-Euler steady state at the step's end
    would be (1 + a dt / 2) times the exponential's, a dt reaching ~2 in a molecular cloud)"""
    return dust_band_rate(band) * band_energy_over_step(band) * EV


def dust_band_absorption(band):
    """Dust absorption of a non-ionizing band (dust_band_power), the band losing it at c_tilde (the kick,
    rt_update_driftkick); the kick passes it to the IR band once, as a source of the IR band's update that the dust
    re-emits (E_abs_tot_toIR), as the dust balance here re-emits it"""
    return Transfer(dust_band_power(band), {band: -rsol / EV, "dust heat": 1}, name=f"Dust absorption of {band}",
                    bibliography=["GIZMO rt_utilities.cc rt_kappa, rt_absorb_frac_albedo, rt_update_driftkick and "
                                  "dust_dE_cooling"])


def dust_ir_absorption():
    """Dust absorption of the IR band at the dust's absorption opacity kappa(T_dust, T_rad) (dust_dE_cooling, 1335;
    rt_kappa_adaptive_IR_band flags -1, 1)"""
    power = C_LIGHT * ir_dust_opacity(T_dust, T_rad) * rho * n_(IR) * EV
    return Transfer(power, {IR: -rsol / EV, "dust heat": 1}, name="Dust absorption of photon_IR",
                    bibliography=["GIZMO rt_utilities.cc dust_dE_cooling", "2003A&A...410..611S"])


def kick_gas_share():
    """The share 2 (1 - exp(-x/2)) / x of the gas IR absorption rate GIZMO's two half-step kicks give the gas: each
    absorbs at most the band's energy, e0 (1 - exp(-x/2)), and gives the gas its opacity share of that, straight into
    its internal energy (rt_update_driftkick). x = c_tilde
    kappa rho dt at the dust's absorption opacity (which dominates the band's) and the start-of-step dust temperature,
    so a function of parameters only"""
    a = rsol * C_LIGHT * ir_dust_opacity(T_dust_initial, T_rad) * rho * dt
    x = sp.Min(a, 100)
    return sp.Piecewise((1 - x / 4 + x * x / 24, x < 1e-3), (2 * (1 - sp.exp(-x / 2)) / a, True))


def gas_ir_absorption():
    """Gas absorption of the IR band at the non-dust absorption opacity (rt_kappa_adaptive_IR_band flags -1, -1), as
    GIZMO's kick and dust balance take it: each half-step kick heats the gas with its share of the energy the band
    absorbs, at c_tilde/c of the physical rate (1:1 with the band's loss) and at most the band's energy per kick
    (kick_gas_share); rt_eqm_dust_temp counts the whole absorption at the true c as the dust's heating (its absorbed
    power is the band's at the gas and dust absorption opacity, flags -1, 0), while the kick re-emits only the dust's
    share into the band: the dust balance's re-emission of the gas share is taken back from the band. Reproduced,
    though the dust heating creates energy. In an optically thick cell the band is absorbed and re-emitted many times
    over a step; uncapped, the gas share would drain it at every pass and the kick's per-kick cap matters"""
    power = C_LIGHT * ir_gas_opacity(T_rad, T_dust) * rho * n_(IR) * EV
    share = kick_gas_share()
    return Transfer(power, {IR: -rsol * (share + 1) / EV, "heat": rsol * share, "dust heat": 1},
                    name="Gas absorption of photon_IR",
                    bibliography=["GIZMO rt_utilities.cc rt_update_driftkick, rt_eqm_dust_temp and "
                                  "rt_kappa_adaptive_IR_band"])


def dust_ir_emission():
    """Thermal emission of the dust into the IR band: 4 sigma kappa_P(T_dust) rho T_dust^4 (dust_dEdt and dust_dE_cooling,
    the emission opacity rt_kappa_adaptive_IR_band(Td, Td, 1, 1))"""
    power = 4 * SIGMA_SB * ir_dust_opacity(T_dust, T_dust) * rho * T_dust**4
    return Transfer(power, {"dust heat": -1, IR: rsol / EV}, name="Dust emission into photon_IR",
                    bibliography=["GIZMO rt_utilities.cc dust_dEdt and dust_dE_cooling"])


def gas_dust_collisions():
    """Gas-dust heat exchange as GIZMO's dust balance writes it (rt_ir_lambdadust, dust_dE_cooling): the coefficient
    gas_dust_heating_coeff (eos.cc 255-263) with the surviving dust at T_dust, without CoolingRate's high-temperature
    truncations, which its RT_INFRARED branch overwrites"""
    coeff = 1.116e-32 * sqrt_T * (1.0 - 0.8 * sp.exp(-75.0 / T)) * Z_dust * dust_survival(T_dust)
    return ThermalTerm(n_Htot * n_Htot * coeff * (T_dust - T), name="Gas-dust collisions",
                       bibliography=["1979ApJS...41..555H", "GIZMO eos.cc gas_dust_heating_coeff"],
                       clumping=sp.Symbol("C_2"), reservoir="dust heat")


# --- Compton -------------------------------------------------------------------------------------------------------

def compton_teff(band):
    """evaluate_Compton_heating_cooling_rate's effective photon temperature of a band (cooling.cc 2387-2440)"""
    return {EUV: 2340 * hnu_EUV, FUV: 24400.0, NUV: 12000.0, ONIR: 2800.0, IR: T_rad}[band]


def compton_off_bands(bands):
    """Compton heating and cooling of the gas off the RT bands: 2.16e-35 n u_k[eV cm^-3] (T - T_eff,k) per unit
    volume, with n the free electrons where T_eff < 3e4 K and the H nuclei otherwise; GIZMO does not take it from the
    bands"""
    heat = sp.S.Zero
    for b in bands:
        teff = compton_teff(b)
        colliders = sp.Piecewise((n_("e-"), teff < 3e4), (n_Htot, True)) if isinstance(teff, sp.Expr) and \
            teff.free_symbols else (n_("e-") if float(teff) < 3e4 else n_Htot)
        heat -= 2.16e-35 * colliders * band_energy_eV(b) * (T - teff)
    return ThermalTerm(heat, name="Inverse Compton cooling (RT bands)",
                       bibliography=["GIZMO cooling.cc evaluate_Compton_heating_cooling_rate"])


# --- inputs of the legacy processes, now expressions of the bands ---------------------------------------------------

def G0_of_band(over_step=False):
    """get_FUV_G0 under RT_PHOTOELECTRIC and M1 RADTRANSFER: the photoelectric band in Habing units, at most 1e8, plus
    1e-56 for GIZMO's floor of MIN_REAL_NUMBER (added rather than a maximum, so that G_0 is smooth at an empty band and
    the photoelectric efficiency's G_0^0.73 differentiable); at the step's end, or over the step
    (band_energy_over_step) for a rate integrated over it"""
    u = band_energy_over_step(FUV) if over_step else band_energy_eV(FUV)
    return 1e-56 + sp.Min(u * EV / U_HABING, 1e8)


def G_LW_of_bands():
    """update_explicit_molecular_fraction's G_LW (cooling.cc 1997-2008): the photoelectric band (over the step,
    band_energy_over_step) as the Lyman-Werner proxy plus the UV background's G_LW_bg, within [1e-10, 1e10], plus the
    IR band's tail above 11.2 eV"""
    G = sp.Min(sp.Max(band_energy_over_step(FUV) * EV / U_HABING + G_LW_bg, 1e-10), 1e10)
    return G + band_energy_eV(IR) * EV * blackbody_fraction(11.2, 500.0, T_rad) / U_HABING


def background_temperature():
    """get_background_radiation_temperature_for_emission_corrections under RT_INFRARED: the IR band and the CMB,
    energy-weighted"""
    e_cmb = 0.262 * (T_CMB / 2.73) ** 4
    e_ir = band_energy_eV(IR)
    return (e_ir * T_rad + e_cmb * T_CMB) / (e_ir + e_cmb)


def ir_self_absorption():
    """CoolingRate's IR self-absorption factor 1/(1 + tau^2), tau = kappa_gas(T, T) rho dx/2 (cooling.cc 1308-1313)"""
    tau = ir_gas_opacity(T, T) * rho * sp.Symbol("Δx") / 2
    return 1 / (1 + tau**2)


def uvb_shielding(nH, gamma_12):
    """return_uvb_shieldfac (Rahmati et al. 2012 form, cooling.cc 2571-2586)"""
    log_T = sp.log(T, 10)
    nss = 0.0123 * sp.Piecewise((gamma_12**0.66, gamma_12 > 0), (1, True)) * 10 ** (0.173 * (log_T - 4))
    q = nH / nss
    return 0.98 / (1 + q**1.64) ** 2.28 + 0.02 / (1 + q * (1 + 1e-4 * nH**4)) ** 0.84


def recombination_return_fraction():
    """The share of the recombination cooling CoolingRate returns to the NUV band (1322): (1 - shieldfac) times the
    band's share of the photoheating, Heat_Ion_from_RHD / (Heat_Ion_from_UVB + Heat_Ion_from_RHD), here per neutral H
    with the UV background's eps_H0_UVB. GIZMO's MIN_REAL_NUMBER guard against no heating at all makes the share jump
    from 0 to 1 at an empty band without a background; here the guard (1e-30 erg/s) is added to both heatings, so the
    share is 1 rather than 0 where neither heats, and smooth"""
    S = uvb_shielding(nH_gizmo_cooling, gamma_12_UVB)
    heat_rhd = eps_HI * C_LIGHT * sigma_HI * n_(EUV)
    return (1 - S) * (heat_rhd + 1e-30) / (eps_H0_UVB * S + heat_rhd + 1e-30)


def ir_radiation_temperature(absorbers, dust_emitters, gas_emitters):
    """(Output expression, sums) of the IR band's radiation temperature after the step, by GIZMO's photon-number
    weighting of what the band keeps and gains (the kick's absorption/re-emission update in rt_update_driftkick and
    the cooling return's in rt_cooling_radiation_to_bands, in one): the band's initial photons that survive absorption
    at T_rad, the rest of the final band at the dust and gas temperatures in proportion to their (positive) emission.
    The surviving share is the kick's exp(-a dt), a the absorption rate at the end of the step: the gross absorption
    and emission of the implicit step can exceed the band's energy many times over in an optically thick cell, so they
    cannot weight it directly. Not reproduced: the kick also counts the gas's opacity share of the absorbed energy,
    which goes to the gas, as photons at max(T_rad, T) (a share below 1e-2 outside dense ionized gas)"""
    A, D, G = sp.symbols("A_IR D_IR G_IR")
    n0 = n_Htot * sp.Symbol("x_photon_IR_initial")
    n1 = n_(IR)
    a_dt = sp.Min(sp.Max(-dt * A, 0) / (n1 + 1e-300), 700)
    unabsorbed = sp.Min(n0 * sp.exp(-a_dt), n1)
    d, g = sp.Max(D, 0), sp.Max(G, 0)
    count = unabsorbed / T_rad + (n1 - unabsorbed) * (d / T_dust + g / T) / (d + g + 1e-300)
    T_new = n1 / (count + 1e-300)
    expr = sp.Max(sp.Min(T, T_dust, T_rad), sp.Min(sp.Max(T, T_dust, T_rad), T_new))
    sums = {"A_IR": {(p, IR): 1 for p in absorbers}, "D_IR": {(p, IR): 1 for p in dust_emitters},
            "G_IR": {(p, IR): 1 for p in gas_emitters}}
    return expr, sums
