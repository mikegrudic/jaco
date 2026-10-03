"""Implementation of recombination process"""

from .chemical_reaction import Reaction
from ..species_strings import add_electron
from ..symbols import T
import sympy as sp


def hydrogenic_recombination_rate(Z):
    """Verner & Ferland 1996"""
    return (
        Z
        * 7.982e-11
        / (
            sp.sqrt(T / 3.148 / Z**2)
            * sp.Pow((1.0 + sp.sqrt(T / 3.148 / Z**2)), 0.252)
            * sp.Pow((1.0 + sp.sqrt(T / 7.036e5 / Z**2)), 1.748)
        )
    )


# Radiative rates: Verner & Ferland 1996. He+ dielectronic: Aldrovandi & Pequignot 1973 (as in KWH96).
He_plus_radiative_recombination_rate = 9.356e-10 / (
    sp.sqrt(T / 4.266e-2) * sp.Pow((1.0 + sp.sqrt(T / 4.266e-2)), 0.2108) * sp.Pow((1.0 + sp.sqrt(T / 3.676e7)), 1.7892)
)
He_plus_dielectronic_recombination_rate = 1.9e-3 * T**-1.5 * sp.exp(-4.7e5 / T) * (1 + 0.3 * sp.exp(-9.4e4 / T))
gasphase_recombination_rates = {
    "H+": hydrogenic_recombination_rate(1),
    "He+": He_plus_radiative_recombination_rate + He_plus_dielectronic_recombination_rate,
    "He++": hydrogenic_recombination_rate(2),
}
# radiative recombination removes the mean kinetic energy of the captured electron
mean_kinetic_energy = 1.036e-16 * T
# dielectronic recombination of He+ removes ~ the n=2 excitation energy of He+ (40.7 eV)
gasphase_recombination_cooling = {
    "H+": mean_kinetic_energy * gasphase_recombination_rates["H+"],
    "He+": mean_kinetic_energy * He_plus_radiative_recombination_rate + 6.526e-11 * He_plus_dielectronic_recombination_rate,
    "He++": mean_kinetic_energy * gasphase_recombination_rates["He++"],
}


def Recombination(ion: str, rate_coefficient=0.0, heat_rate_coefficient=None, *, clumping=None, name="",
                  bibliography=()) -> Reaction:
    """A recombination ion + e- -> neutral + photon at mass-action rate coefficient rate_coefficient, with heat
    heat_rate_coefficient * n_ion * n_e (negative for cooling); see :class:`Reaction`."""
    return Reaction(f"{ion} + e- -> {add_electron(ion)}", rate_coefficient, heat_rate_coefficient=heat_rate_coefficient,
                    clumping=clumping, name=name or f"Recombination of {ion}", bibliography=bibliography)


def GasPhaseRecombination(ion=None) -> Reaction:
    """Gas-phase (radiative and dielectronic) recombination of ion

    Parameters
    ----------
    ion: str, optional
        Ionic species getting recombined. If None, a list of the processes for every ion with a known rate.

    Returns
    -------
    process: Reaction
    """
    if ion is None:
        return [GasPhaseRecombination(s) for s in gasphase_recombination_rates]
    if ion not in gasphase_recombination_rates:
        raise NotImplementedError(f"{ion} does not have an available gas-phase recombination coefficient.")
    return Recombination(ion, gasphase_recombination_rates[ion], heat_rate_coefficient=-gasphase_recombination_cooling[ion],
                         name=f"Gas-phase recombination of {ion}",
                         bibliography=["1996ApJS..103..467V", "1973A&A....25..137A", "1996ApJS..105...19K"])
