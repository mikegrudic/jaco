"""The dust temperature without transported radiation, as a steady-state unknown of the solve.

GIZMO's non-RT cooling (SINGLE_STAR_SINK_DYNAMICS, as every STARFORGE configuration) takes the dust temperature from
get_equilibrium_dust_temperature_estimate (cooling.cc 2335-2380) at each temperature it evaluates, in CoolingRate and in
update_explicit_molecular_fraction: rt_eqm_dust_temp (rt_utilities.cc 1605) finds the root in Td of dust_dEdt
(rt_utilities.cc 1580-1597),

    Lambda_gd nHcgs^2 (T - Td) + A - 4 sigma rho kappa_P(Td, Td) Td^4 / (1 + tau^2),   tau = Sigma kappa_P(Td, Td),

with A the dust's absorption of the CMB and of the attenuated interstellar radiation field (a host input here: it
depends on neither temperature), kappa_P the dust's Planck mean (rt_kappa_adaptive_IR_band(Td, Td, 1, 1): the Semenov
table at T_rad = Td, times 1 - albedo under RADTRANSFER only, the dust-to-metals ratio 1 under SINK_DYNAMICS), and Sigma
the column CoolingRate takes (the cell's limited density gradient; N_H_cool). Without RT_INFRARED the table has no dust
above 1500 K (dust_planck_mean_opacity returns ~0 there), a switch smoothed here like the composition switches
(dust_opacity.py). Td is the solve variable whose equation is the dust heat row, solved with the gas rather than taken
from the start of the step.
"""

import sympy as sp
from jaco.declarations import Parameter, Variable
from jaco.processes.thermal_process import Transfer
from .dust_opacity import semenov_planck_mean, _smoothstep, T_dust, rho
from .symbols import Z_dust, f_dust, n_Htot, sqrt_T, T, lowtemp_truncation, dust_sputtering_truncation

SIGMA_SB = 5.67e-5  # as GIZMO writes the dust emission
PROTONMASS = 1.67262178e-24  # GIZMO's PROTONMASS_CGS, with which the host converts Sigma to N_H_cool
SUBLIMATION_T, SUBLIMATION_HALF_WIDTH = 1500.0, 5.0  # dust_planck_mean_opacity's last zone boundary [K]
DUST_TEMPERATURE_FLOOR = 2.73  # rt_eqm_dust_temp's: the CMB at z = 0
MAX_DUST_TEMP = 1.0e4  # GIZMO's

A_dust = sp.Symbol("A_dust")
N_H_cool = sp.Symbol("N_H_cool")
DUST_PARAMETERS = [
    Parameter("A_dust", "erg cm^-3 s^-1", 0.0,
              "the dust's absorption of the CMB and the interstellar radiation field (get_equilibrium_dust_temperature_"
              "estimate's absorption_rate)"),
    Parameter("rho", "g cm^-3", 2.3e-22, "gas density"),
]
N_H_COOL = Parameter("N_H_cool", "cm^-2", 1e20,
                     "CoolingRate's column (the cell's limited density gradient) in nucleons, Sigma/m_p: the "
                     "optically-thick limit's and the dust emission's optical depth")
TD_VARIABLE = Variable("Td", "dust heat", floor=DUST_TEMPERATURE_FLOOR, ceiling=MAX_DUST_TEMP, units="K",
                       doc="dust temperature: the steady state of the dust's energy balance")


def emission_opacity(T_d=T_dust, albedo=False):
    """kappa_P(Td, Td) per gram of gas [cm^2/g]: the Planck mean at the dust temperature, none above sublimation;
    albedo: times rt_kappa_adaptive_IR_band's 1 - albedo, which it applies under RADTRANSFER"""
    sublimated = _smoothstep((T_d - SUBLIMATION_T + SUBLIMATION_HALF_WIDTH) / (2 * SUBLIMATION_HALF_WIDTH))
    kappa = semenov_planck_mean(T_d, T_d) * (1 - sublimated) * Z_dust * f_dust
    return kappa * (1 - 0.5 / (1 + 725.0**2 / (1 + T_d**2))) if albedo else kappa


def dust_absorption():
    return Transfer(A_dust, {"dust heat": 1}, name="Dust absorption of the background radiation",
                    bibliography=["GIZMO get_equilibrium_dust_temperature_estimate (SINGLE_STAR_SINK_DYNAMICS)"])


def dust_emission(albedo=False):
    """4 sigma rho kappa_P Td^4, reduced by 1 + tau^2 (dust_dEdt's estimate of the outward transport where thick)"""
    kappa = emission_opacity(albedo=albedo)
    tau = kappa * N_H_cool * PROTONMASS
    return Transfer(4 * SIGMA_SB * rho * kappa * T_dust**4 / (1 + tau**2), {"dust heat": -1},
                    name="Dust thermal emission", bibliography=["2003A&A...410..611S", "GIZMO dust_dEdt"])


DUST_RADIATION = ("Dust absorption of the background radiation", "Dust thermal emission")


def dust_side_beyond_truncation():
    """What dust_dEdt's gas-dust term takes beyond CoolingRate's: the dust balance has gas_dust_heating_coeff without
    the high-temperature truncations the gas side carries (gas_dust_collisions.py), so above ~3e4 K the dust gets more
    than the gas loses"""
    coeff = 1.116e-32 * sqrt_T * (1.0 - 0.8 * sp.exp(-75.0 / T)) * Z_dust * f_dust
    return Transfer(n_Htot**2 * coeff * (1 - lowtemp_truncation * dust_sputtering_truncation) * (T - T_dust),
                    {"dust heat": 1}, name="Gas-dust collisions, dust side beyond the truncation",
                    bibliography=["GIZMO dust_dEdt"])


def with_dust_balance(model, extra=(), albedo=False):
    """model with Td a steady-state solve variable on the dust heat row, the dust's absorption and emission added (with
    extra processes; albedo: GIZMO under RADTRANSFER); its Td parameter becomes the variable"""
    params = [p for p in model.parameters if p.name != "Td"]
    names = {p.name for p in params}
    params += [p for p in (*DUST_PARAMETERS, N_H_COOL) if p.name not in names]
    return model.evolve(list(model.processes.values()) + [dust_absorption(), dust_emission(albedo), *extra],
                        solve_vars=(*model.solve_vars, "Td"), parameters=params,
                        variables=(*model.variables, TD_VARIABLE))
