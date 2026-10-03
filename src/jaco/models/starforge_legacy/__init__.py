"""GIZMO's legacy cooling module as a jaco model: the starforge process library with the STARFORGE_LEGACY switches
(jaco/models/starforge/switches.py)."""

from ..starforge.starforge import SOLVE_VARS, TIME_DEPENDENT, GIZMO_FAMILY
from ..starforge.starforge import make_model as _make_model
from ..starforge.switches import STARFORGE_LEGACY


def make_model():
    return _make_model(STARFORGE_LEGACY)
