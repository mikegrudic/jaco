"""Implementation of H_2 formation on dust grain surfaces"""

from ....processes import ChemicalReaction
from ..symbols import sp, T, T_dust, f_dust, Z_dust, H2_formation_heat_cgs

# Formation on dust grains from Hollenbach & McKee 1979
H2_dust_formation_rate = (
    3.0e-18
    * sp.sqrt(T)
    / ((1.0 + 4.0e-2 * sp.sqrt(T + T_dust) + 2.0e-3 * T + 8.0e-6 * T * T) * (1.0 + 1.0e4 / sp.exp(600.0 / T_dust)))
    * f_dust
    * Z_dust
)

grain_formation = ChemicalReaction(
    "H + H -> H_2",
    H2_dust_formation_rate,
    heat_per_reaction=H2_formation_heat_cgs,
    name="Formation of H_2 on dust grains",
    bibliography=["1979ApJS...41..555H"],
)

model_process = grain_formation
