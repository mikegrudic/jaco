"""What a process declares about the photons it absorbs or emits, independently of any band set."""

from dataclasses import dataclass, replace

import numpy as np
import sympy as sp

from .spectra import PowerLaw, Spectral, as_spectral, support

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
    cross_section: number, callable of E, PowerLaw or Spectral, or None
        sigma(E) or kappa(E); zero below its E_min. None: the absorber's band coefficients are all given as overrides
        (e.g. a host code's band-averaged opacities), and those not given are None.
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
    scale: expression, optional
        A state-dependent factor A of the cross-section, in declared symbols (e.g. the dust-to-gas ratio times a
        survival fraction of T_d): the band means sigma_N, chi_E and chi_F are the projection of cross_section times A.
    scattering_scale: expression, optional
        The same for scattering (e.g. x_e for electron scattering); scale by default.
    state: str, optional
        Name of one state variable s (e.g. "Td") on which the cross-section depends otherwise than by a factor:
        cross_section (and scattering) are then functions of (E, s), and the band means are tables in s over
        state_grid (with T_rad, 2-D tables for a tracked band).
    state_grid: sequence of float, optional
        Values of s, log-uniformly spaced, positive; required with state.
    """

    name: str
    cross_section: object
    absorber: str = None
    E_th: float = 0.0
    heat_to: str = "heat"
    heat_yield: object = None
    remainder_to: str = None
    scattering: object = None
    scale: object = 1
    scattering_scale: object = None
    state: str = None
    state_grid: tuple = None

    def __post_init__(self):
        for attr in ("cross_section", "scattering"):
            f = getattr(self, attr)
            if f is None:
                continue
            if self.state is not None and not isinstance(f, (PowerLaw, Spectral)) and callable(f):
                f = Spectral(f, temperature_dependent=True)  # f(E, s)
            object.__setattr__(self, attr, as_spectral(f))
        object.__setattr__(self, "scale", sp.sympify(self.scale))
        if self.scattering_scale is not None:
            object.__setattr__(self, "scattering_scale", sp.sympify(self.scattering_scale))
        if self.state is not None:
            grid = np.asarray(self.state_grid if self.state_grid is not None else (), dtype=float)
            if grid.ndim != 1 or len(grid) < 2 or np.any(grid <= 0):
                raise ValueError(f"absorber {self.name}: state {self.state} needs a positive state_grid")
            steps = np.diff(np.log(grid))
            if not np.allclose(steps, steps[0], rtol=1e-8) or steps[0] <= 0:
                raise ValueError(f"absorber {self.name}: state_grid must be increasing and log-uniform")
            object.__setattr__(self, "state_grid", tuple(float(x) for x in grid))
        elif self.state_grid is not None:
            raise ValueError(f"absorber {self.name}: state_grid without a state")
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

    @property
    def state_symbol(self):
        return sp.Symbol(self.state) if self.state is not None else None

    @property
    def scattering_factor(self):
        return self.scale if self.scattering_scale is None else self.scattering_scale

    def at_state(self, s):
        """The absorber with its state fixed at s: cross-section and scattering functions of E alone, no scale"""
        def fix(f):
            if not isinstance(f, Spectral) or not f.temperature_dependent:
                return f
            return Spectral(lambda E, f=f, s=s: f(E, s), E_min=f.E_min)
        return replace(self, cross_section=fix(self.cross_section), scattering=fix(self.scattering), scale=1,
                       scattering_scale=None, state=None, state_grid=None)


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
class RoutedEmission:
    """Emission whose split among the bands is given rather than derived from a spectrum, taken from the row source:
    overrides (name, band) -> {"fraction": f_b} (any expression) for the bands it reaches; the rest escapes. E.g. a
    host code's fixed routing of a cooling term into bands."""

    name: str
    source: str = "heat"


@dataclass(frozen=True)
class ThermalEmission:
    """Kirchhoff emission 4 pi kappa(E) B_E(T) per unit mass (or per absorber with a cross-section) at the emitter's
    temperature T, taken from the row source"""

    name: str
    kappa: object
    source: str = "dust heat"

    def __post_init__(self):
        object.__setattr__(self, "kappa", as_spectral(self.kappa))
