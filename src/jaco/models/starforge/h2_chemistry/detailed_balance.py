"""Three-body H2 formation by detailed balance with the collisional dissociation rates, for the colliders H, H2 and He.

In local thermodynamic equilibrium formation and dissociation balance at the Saha ratio
n_H^2 / n_H2 = K(T) = (pi m_H k T / h^2)^(3/2) g_H^2 / Z_int(T) exp(-D / kT), with g_H = 4 (electron and nuclear spin),
D the binding energy the model's heat takes (4.48 eV), and Z_int the internal partition function of H2 consistent with
the model's EOS: d ln Z_int / dT = e_int / (k T^2), e_int the rotational-vibrational energy of jaco's H2 fit
(jaco.eos.H2_partition_function), normalized to the exact equilibrium (3:1 nuclear-spin weighted) sum at 3000 K. The
formation coefficient for collider c is then k_3b,c = k_diss,c / K above ~1000 K, so the equilibrium the network relaxes
to there is Saha's. The default network (Forrey 2013 with H as the only third body, against GA08's dissociation) puts its
equilibrium 7.6x above Saha at 2000 K."""

import numpy as np
import sympy as sp
from scipy.integrate import quad
from jaco.processes import Reaction
from jaco.symbols import piecewise_linear
from jaco.eos.H2_partition_function import H2_energy_over_kbT
from ..symbols import T, log_T, H2_BINDING_ENERGY
from .chemical_heat import formation_energy
from . import collisional_dissociation, three_body

K_B, H_PLANCK, M_H = 1.380649e-16, 6.62607015e-27, 1.6735575e-24
THETA_ROT, THETA_VIB, T_NORM = 85.4, 6332.0, 3000.0
LOG_T_GRID = np.linspace(1.0, 5.0, 401)


def _ln_Z_exact(T_):
    J = np.arange(0, 80)
    g = np.where(J % 2 == 0, 1, 3)
    return np.log(np.sum(g * (2 * J + 1) * np.exp(-J * (J + 1) * THETA_ROT / T_)) / (1 - np.exp(-THETA_VIB / T_)))


def ln_Z_int_table():
    """ln Z_int on LOG_T_GRID from the model's H2 energy fit, normalized at T_NORM"""
    e_over_kT = sp.lambdify(T, H2_energy_over_kbT(T), "numpy")
    ref = _ln_Z_exact(T_NORM)
    return np.array([ref + quad(lambda t: (e_over_kT(t) - 1.5) / t, T_NORM, 10**lt, limit=400)[0] for lt in LOG_T_GRID])


def ln_saha_K():
    """ln K(T), n_H^2 / n_H2 in LTE [cm^-3]"""
    ln_Z = piecewise_linear(LOG_T_GRID, ln_Z_int_table(), log_T, name="ln_Z_int_H2")
    return (1.5 * sp.log(sp.pi * M_H * K_B / H_PLANCK**2) + 1.5 * sp.log(T) + sp.log(16) - ln_Z
            - H2_BINDING_ENERGY / (K_B * T))


def ln_dissociation_coefficient(collider):
    """ln k of collisional_dissociation.py's coefficient for collider (GA08's interpolation between the v = 0 and LTE
    limits, at the network's density and composition)"""
    k = collisional_dissociation.H2_collisional_dissociation(collider).rate_coefficient
    return sp.expand_log(sp.log(k), force=True)


# Detailed balance with GA08's rates holds only where their LTE limits do; extrapolated to cold gas, the H-collider LTE
# limit (exp(-43900/T)) against K (exp(-52000/T)) makes the coefficient grow as exp(8100/T). Below SWITCH_T the
# default (Forrey) coefficient is used, blended in log k over a logistic of width SWITCH_DEX in log T.
SWITCH_T, SWITCH_DEX = 1000.0, 0.05


def H2_3body_formation_detailed_balance(chemical_heat=True):
    """2H + c -> H_2 + c for c = H, H_2, He with k = k_diss,c / K above ~1000 K; below, H as the third body at
    Forrey's coefficient and the others off"""
    ln_K = ln_saha_K()
    s = 1 / (1 + sp.exp(-(log_T - sp.log(SWITCH_T, 10)) / SWITCH_DEX))
    ln_cold = {"H": sp.log(three_body.threebody_rate(5)[0]), "H_2": sp.log(1e-60), "He": sp.log(1e-60)}
    return [Reaction(f"H + H + {c} -> H_2 + {c}", sp.exp((1 - s) * ln_cold[c] + s * (ln_dissociation_coefficient(c) - ln_K)),
                     heat_per_reaction=formation_energy(chemical_heat),
                     name=f"3-body formation of H_2 with {c} (detailed balance)",
                     bibliography=["2008MNRAS.388.1627G", "2013ApJ...773L..25F"])
            for c in ("H", "H_2", "He")]
