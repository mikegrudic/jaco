"""Code generation for GIZMO's JACO microphysics solver. CLI: ``python -m jaco.codegen.gizmo.gizmo``."""

from .generate import generate_funcjac_code, TABLE_FILE

__all__ = ["generate_funcjac_code", "TABLE_FILE"]
