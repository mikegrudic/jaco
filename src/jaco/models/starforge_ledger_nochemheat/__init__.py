"""STARFORGE without the heat of H2 formation and dissociation (chemical_heat=False), with the H2 energy ledger's
outputs: the control that isolates the effect of the dissociation heat."""

from ..starforge.starforge import GIZMO_FAMILY, make_model as make_starforge  # noqa: F401  (the host-code family)
from ..starforge.ledger import with_ledger


def make_model():
    """The STARFORGE_LEDGER_NOCHEMHEAT model"""
    return with_ledger(make_starforge(chemical_heat=False))
