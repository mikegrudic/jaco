"""STARFORGE (jaco.models.starforge) with the H2 energy ledger's outputs (jaco.models.starforge.ledger); the solve is
starforge's."""

from ..starforge.starforge import GIZMO_FAMILY, make_model as make_starforge  # noqa: F401  (the host-code family)
from ..starforge.ledger import with_ledger


def make_model():
    """The STARFORGE_LEDGER model"""
    return with_ledger(make_starforge())
