"""starforge_ledger_closed with three-body H2 formation by detailed balance with the collisional dissociation
(h2_chemistry/detailed_balance.py): where the default network's non-Saha equilibrium puts the dissociation plateau."""

from ..starforge.starforge import GIZMO_FAMILY, make_model as make_starforge  # noqa: F401  (the host-code family)
from ..starforge.ledger import with_ledger


def make_model():
    """The STARFORGE_LEDGER_DB_CLOSED model"""
    return with_ledger(make_starforge(formation_heat_closed=True, three_body_detailed_balance=True))
