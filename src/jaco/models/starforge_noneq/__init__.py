"""STARFORGE with time-dependent H and He ionization: jaco.models.starforge with ionization="time-dependent", the
ions integrated with backward-Euler terms like H2 and T instead of solved in equilibrium at each step. Same physics
otherwise; same host-code family."""

from ..starforge import GIZMO_FAMILY  # noqa: F401
from ..starforge.starforge import make_model as _starforge

__all__ = ["make_model", "GIZMO_FAMILY"]


def make_model(**kwargs):
    return _starforge(ionization="time-dependent", **kwargs)
