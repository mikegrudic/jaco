"""Declared parameters and species: kinds drive the reduction, and code generation refuses undeclared symbols."""

import pytest
import sympy as sp

from jaco.codegen.gizmo import generate_funcjac_code
from jaco.declarations import Parameter, Species
from jaco.model import Model
from jaco.processes import Reaction, ThermalTerm, CollisionalIonization, GasPhaseRecombination
from jaco.symbols import n_, x_

n_Htot, T = sp.Symbol("n_Htot"), sp.Symbol("T")
HYDROGEN = [Species(s) for s in ("H", "H+", "e-")]
pdv = ThermalTerm(sp.Symbol("pdv_work"), name="PdV work")


def hydrogen_model(*extra, species=(), parameters=()):
    return Model([CollisionalIonization("H"), GasPhaseRecombination("H+"), pdv, *extra],
                 solve_vars=["u", "T", "H+"], time_dependent=["T"], species=HYDROGEN + list(species),
                 parameters=parameters)


def test_species_kinds_are_checked():
    with pytest.raises(ValueError, match="kind"):
        Species("H", "plasma")
    with pytest.raises(ValueError, match="no element composition"):
        Species("photon_EUV")  # material by default; lowercase p is no element
    assert Species("photon_EUV", "radiation").kind == "radiation"
    assert Species("HD", "trace").kind == "trace"


def test_conflicting_declarations_raise():
    with pytest.raises(ValueError, match="conflicting"):
        Model(parameters=[Parameter("G_0", "Habing", 1.0), Parameter("G_0", "Habing", 2.0)])
    with pytest.raises(ValueError, match="conflicting"):
        hydrogen_model() + Model(species=[Species("e-", "trace")])


def test_kinds_drive_eos_conservation_clumping_and_n_to_x():
    photo = Reaction("H + photon_EUV -> H+ + e-", sp.Symbol("k"), name="photoionization", bibliography=["t"])
    model = hydrogen_model(photo, ThermalTerm(-sp.Symbol("L") * n_("HD") * n_Htot, name="HD cooling"),
                           species=[Species("photon_EUV", "radiation"), Species("HD", "trace")],
                           parameters=[Parameter("k"), Parameter("L")])
    net = model.network
    assert set(net.chemical_species) == {"H", "H+", "e-"}  # EOS and conservation: material only
    effective = {p.name: p for p in model.subprocesses}
    assert effective["photoionization"].clumping == 1  # radiation is not a collider
    assert model.processes["photoionization"].clumping == sp.Symbol("C_2")
    red = net.reduced({"n_Htot"}, [])
    assert x_("photon_EUV") in red["H+"].rhs.free_symbols and n_("photon_EUV") not in red["H+"].rhs.free_symbols
    assert x_("HD") in red["heat"].rhs.free_symbols and n_("HD") not in red["heat"].rhs.free_symbols
    assert dict(red.substitutions)[x_("H")] == 1 - x_("H+")  # HD's H is not in the budget


def test_undeclared_rows_raise():
    model = hydrogen_model(Reaction("H + e- -> H-", sp.Symbol("k"), name="attachment", bibliography=["t"]),
                           parameters=[Parameter("k")])
    with pytest.raises(ValueError, match="undeclared species"):
        model.network


def test_removing_metal_cooling_leaves_the_eos_alone():
    """Metals are declared species, not in the EOS by accident of a rate mentioning them (D5)"""
    from jaco.models.starforge import make_model

    model = make_model()
    assert model.without("Metal line cooling").network.chemical_species == model.network.chemical_species


@pytest.mark.parametrize("extra,undeclared", [
    (ThermalTerm(1e-25 * sp.Symbol("G0") * n_Htot, name="typo"), "G0"),
    (Reaction("H+ + e- + 2H -> 3H", sp.Symbol("k"), name="four-body", bibliography=["t"]), "C_4"),
])
def test_codegen_refuses_undeclared_symbols(tmp_path, extra, undeclared):
    model = hydrogen_model(extra, parameters=[Parameter("k")])
    with pytest.raises(ValueError, match=f"undeclared symbols.*{undeclared}"):
        generate_funcjac_code(model, output_dir=str(tmp_path))


def test_header_declares_fields(tmp_path):
    model = hydrogen_model(ThermalTerm(1e-25 * sp.Symbol("G_0") * n_Htot, name="heating"),
                           parameters=[Parameter("G_0", "Habing", 1.0, "FUV field")])
    result = generate_funcjac_code(model, output_dir=str(tmp_path))
    header = (tmp_path / "microphysics_func_jac.h").read_text()
    for v in result["var_names"]:
        assert f"#define JACO_HAS_VAR_{v}\n" in header
    for p in result["param_names"]:
        assert f"#define JACO_HAS_PARAM_{p}\n" in header
    assert set(result["param_names"]) == {"Delta_t", "G_0", "n_Htot", "pdv_work", "u_initial", "C_2"}
    assert "double G_0;  /* FUV field [Habing] */" in header
    assert '#define JACO_PARAM_NAMES {"' in header


def test_codegen_needs_a_model_with_species(tmp_path):
    with pytest.raises(TypeError):
        generate_funcjac_code(CollisionalIonization("H") + pdv, output_dir=str(tmp_path))
    with pytest.raises(ValueError, match="no species"):
        generate_funcjac_code(Model([pdv], solve_vars=["u", "T"], time_dependent=["T"]), output_dir=str(tmp_path))


def test_solve_assumes_declared_defaults():
    model = Model([CollisionalIonization("H"), GasPhaseRecombination("H+"),
                   Reaction("H -> H+ + e-", rate=sp.Symbol("zeta") * n_("H"), name="CR ionization", bibliography=["t"])],
                  species=HYDROGEN, parameters=[Parameter("zeta", "s^-1", 1e-16, "ionization rate")])
    sol = model.solve({"T": [100.0], "n_Htot": [1.0]}, {"H+": [1e-3]}, time_dependent=[])
    ref = model.solve({"T": [100.0], "n_Htot": [1.0], "zeta": [1e-16]}, {"H+": [1e-3]}, time_dependent=[])
    assert float(sol["H+"][0]) == pytest.approx(float(ref["H+"][0]), rel=1e-6)
