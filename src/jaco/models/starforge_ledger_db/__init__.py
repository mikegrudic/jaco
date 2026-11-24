"""starforge_ledger with three-body H2 formation by detailed balance with the collisional dissociation
(h2_chemistry/detailed_balance.py): the arm that isolates the effect of the network's non-Saha equilibrium."""

from ..starforge.starforge import GIZMO_FAMILY, make_model as make_starforge  # noqa: F401  (the host-code family)
from ..starforge.ledger import with_ledger


def make_model():
    """The STARFORGE_LEDGER_DB model"""
    return with_ledger(make_starforge(three_body_detailed_balance=True))
