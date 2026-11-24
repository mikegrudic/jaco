"""Grain-assisted recombination (Weingartner & Draine 2001), as GIZMO's alpha_recomb_grain (simple_chemistry.cc).

WD01 Eq. 5 defines the coefficient per ion and H nucleus (dn_i/dt = -alpha_g n_i n_H) and Eq. 8 fits it with ln T, in
terms of psi = G sqrt(T) / n_e with G the FUV field in Habing units; the fit holds to 20% for 10 <= T <= 1e3 K and
1e2 <= psi <= 1e6. GIZMO (after Kim+23) adds 50 to psi to keep it inside that range in shielded gas.
"""

import sympy as sp
from jaco.processes import Reaction
from jaco.symbols import n_
from .symbols import T, G_0, Z_dust, f_dust, n_Htot, x_, dust_sputtering_truncation

coeffs = {  # fit parameters from 2001ApJ...563..842W
    "H+": [12.25, 8.074e-6, 1.378, 5.087e2, 1.586e-2, 0.4723, 1.102e-5],  #  H+
    "He+": [5.572, 3.185e-7, 1.512, 5.115e3, 3.903e-7, 0.4956, 5.494e-7],  #  He+
    "C+": [45.58, 6.089e-3, 1.128, 4.331e2, 4.845e-2, 0.8120, 1.333e-4],  #  C+
    "Na+": [2.178, 1.732e-7, 2.133, 1.029e4, 1.859e-6, 1.0341, 3.223e-5],  #  Na+
    "Mg+": [2.510, 8.116e-8, 1.864, 6.170e4, 2.169e-6, 0.9605, 7.232e-5],  #  Mg+
    "Si+": [2.166, 5.678e-8, 1.874, 4.375e4, 1.635e-6, 0.8964, 7.538e-5],  #  Si+
    "S+": [3.064, 7.769e-5, 1.319, 1.087e2, 3.475e-1, 0.4790, 4.689e-2],  #  S+
    "K+": [1.596, 1.907e-7, 2.123, 8.138e3, 1.530e-5, 1.0380, 4.550e-5],  #  K+
    "Ca+": [1.636, 8.208e-9, 2.289, 1.254e5, 1.349e-9, 1.1506, 7.204e-4],  #  Ca+
    "Mn+": [2.029, 1.433e-6, 1.673, 1.403e4, 1.865e-6, 0.9358, 4.339e-9],  #  Mn+
    "Fe+": [1.701, 9.554e-8, 1.851, 5.763e4, 4.116e-8, 0.9456, 2.198e-5],  #  Fe+
    "Ca++": [8.270, 2.051e-4, 1.252, 1.590e2, 6.072e-2, 0.5980, 4.497e-7],  #  Ca++
}


def grain_charge_psi(x_e, n_H=n_Htot):
    """GIZMO's grain_charge_psi: G0 sqrt(T) / n_e, plus 50 to keep it finite in shielded gas (Kim+2023)"""
    return G_0 * sp.sqrt(T) / (n_H * x_e) + 50


def alpha_grain(ion, x_e, n_H=n_Htot):
    """Rate coefficient (cm^3 s^-1) per ion and H nucleus: the recombination rate per volume is alpha n_ion n_H"""
    if ion not in coeffs:
        raise NotImplementedError(f"idk the grain-assisted recombination coefficient for {ion}.")
    C = coeffs[ion]
    psi = grain_charge_psi(x_e, n_H)
    return (Z_dust * f_dust * 1e-14 * C[0]
            / (1 + C[1] * psi ** C[2] * (1 + C[3] * T ** C[4] * psi ** (-C[5] - C[6] * sp.log(T)))))


def GrainAssistedRecombination(ion):
    """ion + grain -> neutral at alpha n_ion n_H C_2 per volume, the grain supplying the electron. The released energy
    goes to the grain, so the gas loses no heat. Grains are taken to be sputtered away above 3e5 K, as in the gas-dust
    coupling."""
    return Reaction(f"{ion} -> {ion.rstrip('+')}",
                    rate=alpha_grain(ion, x_("e-")) * n_(ion) * n_Htot * dust_sputtering_truncation,
                    clumping=sp.Symbol("C_2"), name=f"Grain-assisted recombination of {ion}",
                    bibliography=["2001ApJ...563..842W"])
