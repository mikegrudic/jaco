"""Three-body H2 formation by detailed balance with the collisional dissociation rates, for the colliders H, H2 and He.

In local thermodynamic equilibrium formation and dissociation balance at the Saha ratio
n_H^2 / n_H2 = K(T) = (pi m_H k T / h^2)^(3/2) g_H^2 / Z_int(T) exp(-D / kT), with g_H = 4 (electron and nuclear spin),
D the binding energy the model's heat takes (4.48 eV), and Z_int the internal partition function of H2 consistent with
the model's EOS: d ln Z_int / dT = e_int / (k T^2), e_int the rotational-vibrational energy of jaco's H2 fit
(jaco.eos.H2_partition_function), normalized to the exact equilibrium (3:1 nuclear-spin weighted) sum at 3000 K. The
formation coefficient for collider c is then k_3b,c = k_diss,c^LTE / K, so the equilibrium the network relaxes to is
Saha's. The default network (Forrey 2013 with H as the only third body, against GA08's dissociation) puts its
equilibrium 7.6x above Saha at 2000 K."""

import numpy as np
import sympy as sp
from scipy.integrate import quad
from jaco.processes import Reaction
from jaco.symbols import piecewise_linear
from jaco.eos.H2_partition_function import H2_energy_over_kbT
from ..symbols import T, log_T, H2_BINDING_ENERGY
from .chemical_heat import formation_heat

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


LN10 = sp.log(10.0)
LN_K_DISS_LTE = {  # GA08's LTE limits, as collisional_dissociation.py writes them
    "H": sp.log(3.52e-9) - 43900.0 / T,
    "H_2": sp.log(1.3e-9) - 53300.0 / T,
    "He": LN10 * (-2.729 - 1.75 * log_T - 23474.0 / T),
}


def H2_3body_formation_detailed_balance(chemical_heat=True):
    """2H + c -> H_2 + c for c = H, H_2, He with k = k_diss,c^LTE / K"""
    ln_K = ln_saha_K()
    return [Reaction(f"H + H + {c} -> H_2 + {c}", sp.exp(LN_K_DISS_LTE[c] - ln_K),
                     heat_per_reaction=formation_heat("3-body", chemical_heat),
                     name=f"3-body formation of H_2 with {c} (detailed balance)",
                     bibliography=["2008MNRAS.388.1627G"])
            for c in ("H", "H_2", "He")]
