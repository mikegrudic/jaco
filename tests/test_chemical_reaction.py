from jaco.processes.chemical_reaction import ChemicalReaction
from jaco.processes.recombination import GasPhaseRecombination, gasphase_recombination_rates, mean_kinetic_energy
import sympy as sp
from jaco.symbols import T, n_


def test_chemical_reaction():
    equation = "2Tom + 3Dick++ + Harry -> 3Gerald"
    rate = 12345 * sp.exp(-69.0 / T)
    reaction = ChemicalReaction(equation, rate, bibliography=["The Great Escape"])
    print(reaction.network["Gerald"].rhs)
    assert reaction.rate == rate * n_("Tom") ** 2 * n_("Dick++") ** 3 * n_("Harry") * sp.Symbol("C_6")
    assert reaction.network["Gerald"].rhs == 3 * reaction.rate
    assert reaction.network["Tom"].rhs == -2 * reaction.rate
    assert reaction.network["Dick++"].rhs == -3 * reaction.rate
    assert reaction.network["Harry"].rhs == -1 * reaction.rate

    reaction2 = ChemicalReaction(
        "H+ + e- -> H",
        rate_coefficient=gasphase_recombination_rates["H+"],
        heat_per_reaction=-mean_kinetic_energy,
        bibliography=["Osterbrock 1987"],
    )
    assert reaction2.network == GasPhaseRecombination("H+").network


if __name__ == "__main__":
    test_chemical_reaction()
