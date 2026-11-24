"""Tables are sympy atoms carrying their data: code generation collects a model's tables from its expressions, with no
process-wide registry."""

import numpy as np
import pytest
import sympy as sp

import jaco.interpolation as interpolation
from jaco.interpolation import Table, TableInterp2D, tables_in
from jaco.symbols import table_interp_2d

x, y = sp.symbols("x y")
AXES = [np.logspace(0, 2, 3), np.logspace(0, 3, 4)]


def lookup(name, scale=1.0):
    return table_interp_2d(name, scale * np.arange(12.0).reshape(3, 4), AXES, x, y, log_axes=[True, True])


def test_no_registry():
    assert not hasattr(interpolation, "_TABLE_REGISTRY") and not hasattr(interpolation, "register_table")


def test_tables_are_atoms_with_their_data():
    f = lookup("tab")
    (table,) = f.atoms(Table)
    assert table.name == "tab" and table.info["shape"] == (3, 4)
    assert float(f.subs({x: 10.0, y: 10.0})) == pytest.approx(5.0)  # evaluates from the atom's own data
    assert lookup("tab") == f and lookup("tab", 2.0) != f  # equal only if the contents are
    assert sp.diff(f, x).atoms(Table) == {table}


def test_collected_from_expressions_in_creation_order():
    a, b = lookup("first"), lookup("second")
    assert list(tables_in([b + a])) == ["first", "second"]
    assert tables_in([x * y]) == {}


def test_same_name_different_data_raises():
    with pytest.raises(ValueError, match="two different tables are named clash"):
        tables_in([lookup("clash") + lookup("clash", 2.0)])


def test_metal_line_tables_are_named_by_slice():
    """z = 0 and z = 0.3 are different slices with different names; z = 0.01 is the z = 0 slice (D7)"""
    from jaco.models.starforge import metal_line_cooling as mlc

    if not mlc._hdf5_path().is_file():
        pytest.skip("spcool_tables.hdf5 not available")
    names = {z: {t.name for t in mlc.metal_line_cooling_rate("C", "Carbon_cooling", z).atoms(Table)} for z in (0, 0.01, 0.3)}
    assert names[0] == names[0.01] == {"Carbon_cooling_z0"} and names[0.3] != names[0]
    assert len(tables_in([mlc.metal_line_cooling_rate("C", "Carbon_cooling", z) for z in (0, 0.01, 0.3)])) == 2


def test_tables_of_another_model_do_not_leak(tmp_path):
    """Building a model with tables first must not add them to a table-free model's generated code"""
    from jaco.models.starforge import metal_line_cooling as mlc
    from jaco.models.wind_comparison import make_model
    from jaco.codegen.gizmo import generate_funcjac_code

    if mlc._hdf5_path().is_file():
        mlc.metal_line_cooling()
    lookup("stray")
    generate_funcjac_code(make_model(), output_dir=str(tmp_path), language="c", source_ext=".cc")
    assert "jaco_tables.h" not in (tmp_path / "microphysics_func_jac.cc").read_text()
    assert "extern JacoTable" not in (tmp_path / "jaco_tables.h").read_text()
    assert not (tmp_path / "jaco_tables.hdf5").exists()
