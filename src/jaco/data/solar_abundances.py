from .atomic_weights import atomic_weights


class SolarAbundancesClass:
    """Container for solar abundances with methods to convert between mass fraction and abundance per H"""

    bibliography = ["2009ARA&A..47..481A"]

    @property
    def mass_fraction(self):
        """Returns a hard-coded dict of Solar abundance mass fractions"""
        return {
            "Z": 0.0142,
            "He": 0.27030,
            "C": 2.53e-3,
            "N": 7.41e-4,
            "O": 6.13e-3,
            "Ne": 1.34e-3,
            "Mg": 7.57e-4,
            "Si": 7.12e-4,
            "S": 3.31e-4,
            "Ca": 6.87e-5,
            "Fe": 1.38e-3,
        }

    @property
    def abundance_per_H(self):
        """Returns dictionary of abundances per H nucleon"""
        return {species: f / (1 - f) / atomic_weights[species] for species, f in self.mass_fraction.items()}

    def x(self, species):
        return self.get_abundance(species)

    def get_mass_fraction(self, species: str) -> float:
        """Returns the mass fraction of a given species"""
        if species in self.mass_fraction:
            return self.mass_fraction[species]
        else:
            raise NotImplementedError(f"Solar abundance of {species} not available.")

    def get_abundance(self, species: str) -> float:
        """Returns the abundance per H nuclear of an input species"""
        if species in self.mass_fraction:
            return self.abundance_per_H[species]
        else:
            raise NotImplementedError(f"Solar abundance of {species} not available.")


SolarAbundances = SolarAbundancesClass()
