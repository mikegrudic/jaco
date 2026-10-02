"""GIZMO's cooling module carries no heat for H2 formation or collisional dissociation."""

import pytest
from ..h2_chemistry.H2_chemistry import h2_chemistry_processes


@pytest.mark.parametrize("process", h2_chemistry_processes, ids=lambda p: p.name)
def test_no_H2_chemical_heat(process):
    assert process.heat == 0
