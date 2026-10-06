"""What a process declares about the photons it absorbs or emits, independently of any band set."""

from dataclasses import dataclass

from .spectra import as_spectral, support

HEAT_ROWS = ("heat", "dust heat")


@dataclass(frozen=True)
class Absorber:
    """A process absorbing photons with cross-section sigma(E) [cm^2] per absorber, or opacity kappa(E) [cm^2 g^-1].

    Each absorbed photon of energy E gives E_th to chemistry (e.g. the ionization energy) and E - E_th to heat_to;
    with a yield Y(E), Y (E - E_th) goes to heat_to and (1 - Y) (E - E_th) to remainder_to (photoelectric heating:
    the gas takes what the ejected electrons carry, the grain the rest).

    Parameters
    ----------
    name: str
        Id of the process.
    cross_section: number, callable of E, PowerLaw or Spectral
        sigma(E) or kappa(E); zero below its E_min.
    absorber: str, optional
        The absorbing species (sigma per particle), or None for an opacity per unit mass.
    E_th: float
        Energy per absorbed photon to chemistry [eV]; also a lower bound of the absorption.
    heat_to: str
        Row taking (the yield's share of) E - E_th: "heat" (the gas) or "dust heat".
    heat_yield: number or callable of E, optional
        Y(E) in [0, 1].
    remainder_to: str, optional
        Row taking (1 - Y)(E - E_th); required with heat_yield.
    scattering: number or callable of E, optional
        Scattering cross-section or opacity, in the flux mean (momentum) only.
    """

    name: str
    cross_section: object
    absorber: str = None
    E_th: float = 0.0
    heat_to: str = "heat"
    heat_yield: object = None
    remainder_to: str = None
    scattering: object = None

    def __post_init__(self):
        object.__setattr__(self, "cross_section", as_spectral(self.cross_section))
        if self.scattering is not None:
            object.__setattr__(self, "scattering", as_spectral(self.scattering))
        if self.heat_yield is not None:
            object.__setattr__(self, "heat_yield", as_spectral(self.heat_yield))
            if self.remainder_to is None:
                raise ValueError(f"absorber {self.name}: a heat_yield needs a remainder_to row")
        for row in (self.heat_to, self.remainder_to):
            if row is not None and row not in HEAT_ROWS:
                raise ValueError(f"absorber {self.name}: heat rows are {HEAT_ROWS}, not {row!r}")

    @property
    def E_min(self):
        """Lowest absorbed photon energy"""
        return max(self.E_th, support(self.cross_section))

    def sigma(self, E):
        return self.cross_section(E) * (E >= self.E_min)


@dataclass(frozen=True)
class Line:
    """Emission at one photon energy E0 [eV], taken from the row source"""

    name: str
    E0: float
    source: str = "heat"


@dataclass(frozen=True)
class Continuum:
    """Emission with spectrum j(E) or j(E, T) (energy per unit photon energy, any normalization), taken from the row
    source; E_range bounds the integral that normalizes it.

    j: number, callable or Spectral (temperature_dependent for j(E, T))
    """

    name: str
    j: object
    source: str = "heat"
    E_range: tuple = (1e-6, 1e6)

    def __post_init__(self):
        object.__setattr__(self, "j", as_spectral(self.j))

    @property
    def temperature_dependent(self):
        return getattr(self.j, "temperature_dependent", False)


@dataclass(frozen=True)
class ThermalEmission:
    """Kirchhoff emission 4 pi kappa(E) B_E(T) per unit mass (or per absorber with a cross-section) at the emitter's
    temperature T, taken from the row source"""

    name: str
    kappa: object
    source: str = "dust heat"

    def __post_init__(self):
        object.__setattr__(self, "kappa", as_spectral(self.kappa))
