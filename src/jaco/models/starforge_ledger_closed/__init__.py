"""STARFORGE whose H2 formed on grains or through H- gives the gas the whole binding energy where collisions
de-excite it (formation_heat_closed), with the H2 energy ledger's outputs: the dissociate-and-reform cycle conserves
the gas's energy, so a second collapse can only come from net dissociation."""

from ..starforge.starforge import GIZMO_FAMILY, make_model as make_starforge  # noqa: F401  (the host-code family)
from ..starforge.ledger import with_ledger


def make_model():
    """The STARFORGE_LEDGER_CLOSED model"""
    return with_ledger(make_starforge(formation_heat_closed=True))
