"""Implementation of molecular hydrogen chemistry following Hopkins+22"""

import sympy as sp
from . import grain_formation
from ....processes import ChemicalReaction

# reactions = []

#    "H + H -> H_2": (H2_dust_formation_rate,,
#    "H_2 + H -> 3H": 6.67e-12 * sqrt_T * sp.exp(1 + 63590 / T),

reactions = [grain_formation.process]


# # Collisional dissociation
# collisional_dissociation = ChemicalReaction(
#     "H_2 + H -> 3H",
#     6.67e-12 * sqrt_T * sp.exp(1 + 63590 / T),
#     name="Collisional dissociation of H_2 by H",
#     bibliography=["1986ApJ...302..585M", "1983ApJ...270..578L"],
# )

# # k1 H + e- -> H- + 0.754195eV 1979MNRAS.187P..59W


# if(T<=6000.) {k1=-17.845 + 0.762*log_T + 0.1523*log_T*log_T - 0.03274*log_T*log_T*log_T;} else {k1=-16.420 + 0.1998*log_T*log_T - 5.447e-3*log_T*log_T*log_T*log_T + 4.0415e-5*log_T*log_T*log_T*log_T*log_T*log_T;}

# k5 = 5.7e-6 / sqrt_T + 6.3e-8 - 9.2e-11 * sqrt_T + 4.4e-13 * T

# lnTeV = sp.ln(T / 8.617e-5)
# k15 = sp.exp(
#     -1.801849334e1
#     + 2.36085220e0 * lnTeV
#     - 2.827443e-1 * lnTeV * lnTeV
#     + 1.62331664e-2 * lnTeV * lnTeV * lnTeV
#     - 3.36501203e-2 * lnTeV * lnTeV * lnTeV * lnTeV
#     + 1.17832978e-2 * lnTeV * lnTeV * lnTeV * lnTeV * lnTeV
#     - 1.65619470e-3 * lnTeV * lnTeV * lnTeV * lnTeV * lnTeV * lnTeV
#     + 1.06827520e-4 * lnTeV * lnTeV * lnTeV * lnTeV * lnTeV * lnTeV * lnTeV
#     - 2.63128581e-6 * lnTeV * lnTeV * lnTeV * lnTeV * lnTeV * lnTeV * lnTeV * lnTeV
# )
