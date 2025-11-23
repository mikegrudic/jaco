from .nbody_process import NBodyProcess


class ChemicalReaction(NBodyProcess):
    """Process that implements a chemical reaction described by an equation and a rate coefficient"""

    def __init__(self, equation, rate_coefficient, heat_per_reaction=0, bibliography=[]):
        """
        Chemical reaction process

        Parameters
        ----------
        equation: string
            Chemical equation for the reaction. Syntax: <c1>reactant1 + <c2>reactant2... -> <c3>product1 + <c4>production2...
            examples: '3H -> H_2 + H', 'H+ + e- -> H'
        """
        self.name = equation
        self.bibliography = bibliography
        self.lhs_coeffs, self.rhs_coeffs = self.species_and_coeffs(equation)
        self.colliding_species = sum(
            [c * [s] for s, c in self.lhs_coeffs.items()], start=[]
        )  # repeating in the list based on multiplcity so we get the correct powers in the rate expression
        self.rate_coefficient = rate_coefficient
        self.initialize_network()
        self.update_network()
        self.heat_rate_coefficient = rate_coefficient * heat_per_reaction

    def update_network(self):
        """Sets up rate terms in the associated chemistry network for each species involved"""
        rate = self.rate
        for s, coeff in self.lhs_coeffs.items():
            self.network[s] -= rate * coeff
        for s, coeff in self.rhs_coeffs.items():
            self.network[s] += rate * coeff

    @staticmethod
    def species_and_coeffs(equation) -> list[dict]:
        """Returns the lists of species and corresponding stoichiometric coefficients from the equation"""
        if "->" not in equation:
            raise ValueError(f"Chemical equation {equation} has no ->")

        coeffs_dicts = []
        for idx, side in enumerate(("lhs", "rhs")):
            terms = equation.split("->")[idx].split(" + ")
            terms = [t.strip() for t in terms]
            coefficients = len(terms) * [1]
            species = terms.copy()

            # strip off leading coefficients on the species
            for i, species_string in enumerate(terms):
                for j in range(1, len(species_string)):
                    substr = species_string[:j]
                    if substr.isnumeric():
                        coefficients[i] = int(substr)
                        species[i] = species_string[j:]
                    else:
                        break

            coeffs_dicts.append(dict(zip(species, coefficients)))

        return coeffs_dicts
