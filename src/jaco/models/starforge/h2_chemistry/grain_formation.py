"""Implementation of H_2 formation on dust grain surfaces"""

from ....processes import Reaction
from ....symbols import n_
from ..symbols import sp, T, T_dust, f_dust, Z_dust, n_Htot
from .chemical_heat import formation_heat, BIBLIOGRAPHY

# Formation on dust grains from Hollenbach & McKee 1979
H2_dust_formation_rate = (
    3.0e-18
    * sp.sqrt(T)
    / ((1.0 + 4.0e-2 * sp.sqrt(T + T_dust) + 2.0e-3 * T + 8.0e-6 * T * T) * (1.0 + 1.0e4 / sp.exp(600.0 / T_dust)))
    * f_dust
    * Z_dust
)


def grain_formation(chemical_heat=True):
    """H + H -> H_2 on grains: the rate per volume is R n_H,tot n_HI C_2, since the dust abundance scales with all of
    the gas (Hollenbach & McKee 1979; Glover & Jappsen 2007), not R n_HI^2"""
    return Reaction(
        "H + H -> H_2",
        rate=H2_dust_formation_rate * n_Htot * n_("H"),
        clumping=sp.Symbol("C_2"),
        heat_per_reaction=formation_heat("grain", chemical_heat),
        name="Formation of H_2 on dust grains",
        bibliography=["1979ApJS...41..555H", "2007ApJ...666....1G"] + (BIBLIOGRAPHY if chemical_heat else []),
    )
