"""Implementation of ionization process"""

from .chemical_reaction import Reaction
from ..species_strings import remove_electron
from ..symbols import T, T5
import sympy as sp
from astropy import units as u


def Ionization(species: str, rate=None, *, rate_coefficient=None, collider=None, heat_per_reaction=0, clumping=None,
               name="", bibliography=()) -> Reaction:
    """An ionization: species -> species+ + e-, or species + collider -> species+ + e- + collider.

    Collisional, photo- or cosmic-ray ionization alike; rate is the events per unit volume and time (before clumping),
    or rate_coefficient the mass-action k. See :class:`Reaction` for the other arguments.
    """
    ion = remove_electron(species)
    if collider is None:
        equation = f"{species} -> {ion} + e-"
    elif collider == "e-":
        equation = f"{species} + e- -> {ion} + 2e-"
    else:
        equation = f"{species} + {collider} -> {ion} + e- + {collider}"
    return Reaction(equation, rate_coefficient, heat_per_reaction, rate=rate, clumping=clumping,
                    name=name or f"Ionization of {species}", bibliography=bibliography)


def ionization_energy(species, unit=u.erg):
    """Return the energy in erg required to ionize a species"""
    # NOTE: come back and get this from a proper datafile
    energies_eV = {"H": 13.6, "He": 24.59, "He+": 54.42}
    return energies_eV[species] * u.eV.to(unit)


collisional_ionization_cooling_rates = {
    "H": 1.27e-21 * sp.sqrt(T) * sp.exp(-157809.1 / T) / (1 + sp.sqrt(T5)),  # 1996ApJS..105...19K
    "He": 9.38e-22 * sp.sqrt(T) * sp.exp(-285335.4 / T) / (1 + sp.sqrt(T5)),  # 1996ApJS..105...19K
    "He+": 4.95e-22 * sp.sqrt(T) * sp.exp(-631515 / T) / (1 + sp.sqrt(T5)),  # 1996ApJS..105...19K
}

collisional_ionization_rates = {
    "H": 5.85e-11 * sp.sqrt(T) * sp.exp(-157809.1 / T) / (1 + sp.sqrt(T5)),  # 1996ApJS..105...19K
    "He": 2.38e-11 * sp.sqrt(T) * sp.exp(-285335.4 / T) / (1 + sp.sqrt(T5)),  # 1996ApJS..105...19K
    "He+": 5.68e-12 * sp.sqrt(T) * sp.exp(-631515 / T) / (1 + sp.sqrt(T5)),  # 1996ApJS..105...19K
}


def CollisionalIonization(species=None, clumping=None):
    """Collisional ionization of species by electrons, species + e- -> species+ + 2e-.

    Parameters
    ----------
    species: str, optional
        Species being collisionally ionized. If None, a list of the processes for every species with a known rate.
    clumping: optional
        Factor <n^2>/<n>^2 multiplying the two-body rate; C_2 by default.

    Returns
    -------
    process: Reaction
    """
    if species is None:
        return [CollisionalIonization(s, clumping) for s in collisional_ionization_rates]
    if species not in collisional_ionization_rates:
        raise NotImplementedError(f"{species} does not have an available collisional ionization coefficient.")
    return Ionization(species, rate_coefficient=collisional_ionization_rates[species], collider="e-",
                      heat_per_reaction=-ionization_energy(species), clumping=clumping,
                      name=f"Collisional Ionization of {species}", bibliography=["1996ApJS..105...19K"])
