import sympy as sp
from jaco.model import Model
from jaco.processes import Reaction, ThermalTerm
from jaco.declarations import Species
from jaco.symbols import T, n_, n_Htot, piecewise_powerlaw
from jaco.math import logistic

EV = 1.602e-12
cooling = ThermalTerm(-piecewise_powerlaw([10, 1e4, 1e5, 1e8, 1e10], [1e-30, 1e-23, 1e-21, 1e-22, 1e-21], T,
                                          extrapolate=True, name="lambda_toy") * n_Htot ** 2,
                      name="toy cooling", bibliography=["made up"])
heating = ThermalTerm(2e-26 * n_("H") * logistic(-(T - 1.5e4) / 1e3), name="toy heating",
                      bibliography=["made up"])  # a heating floor that switches off above ~1.5e4 K
pdv = ThermalTerm(sp.Symbol("pdv_work"), name="PdV work")  # required: the host's compressional heating
ionization = Reaction("H + e- -> H+ + e- + e-", 5.85e-11 * sp.sqrt(T) * sp.exp(-157809.1 / T),
                      heat_per_reaction=-13.6 * EV, clumping=1,
                      name="H collisional ionization", bibliography=["1992ApJS...78..341C"])
recombination = Reaction("H+ + e- -> H", 2.59e-13 * (T / 1e4) ** -0.7, clumping=1,
                         name="H+ recombination", bibliography=["1992ApJS...78..341C"])
cosmic_rays = Reaction("H -> H+ + e-", 1e-16, heat_per_reaction=10 * EV,
                      name="cosmic-ray ionization of H", bibliography=["made up"])  # an electron source besides H+ itself


def make_model():
    return Model([cooling, heating, pdv, ionization, recombination, cosmic_rays],
                 solve_vars=["u", "T", "H+"], time_dependent=["T"],
                 species=[Species("H"), Species("H+"), Species("e-")])
