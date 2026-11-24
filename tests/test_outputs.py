"""Outputs: the rows of the radiation and energy species that are not solved for, and declared Outputs, evaluated by a
generated Jacobian-free function at a state, without touching the solved system."""

import ctypes
import subprocess

import numpy as np
import pytest
import sympy as sp

from jaco.codegen.gizmo import generate_funcjac_code
from jaco.declarations import Output, Parameter, Species
from jaco.model import Model, Rule
from jaco.processes import CollisionalIonization, GasPhaseRecombination, Reaction, ThermalTerm
from jaco.symbols import n_

n_Htot, T = sp.Symbol("n_Htot"), sp.Symbol("T")
pdv = ThermalTerm(sp.Symbol("pdv_work"), name="PdV work")
SPECIES = [Species("H"), Species("H+"), Species("e-"), Species("dust heat", "energy", "energy the gas gives the dust"),
           Species("photon_EUV", "radiation", "ionizing photons")]
PARAMS = [Parameter("Td", "K", 15.0, "dust temperature"), Parameter("sigma_c", "cm^3 s^-1", 1e-8, "c sigma")]
dust = ThermalTerm(-1e-33 * n_Htot**2 * sp.sqrt(T) * (T - sp.Symbol("Td")), name="gas-dust", reservoir="dust heat")
photo = Reaction("H + photon_EUV -> H+ + e-", sp.Symbol("sigma_c"), 2e-12, name="photoionization", bibliography=["t"])


def model(*extra, outputs=(), rules=()):
    return Model([CollisionalIonization("H"), GasPhaseRecombination("H+"), dust, photo, pdv, *extra],
                 solve_vars=["u", "T", "H+"], time_dependent=["T"], species=SPECIES, parameters=PARAMS,
                 outputs=outputs, rules=rules)


L_rec = Output("L_rec", units="erg cm^-3 s^-1", doc="recombination cooling",
               heat_of={"Gas-phase recombination of H+": -1})


def generate(m, path):
    return generate_funcjac_code(m, output_dir=str(path))


def test_species_rows_and_declared_outputs(tmp_path):
    result = generate(model(outputs=[L_rec]), tmp_path)
    assert result["output_names"] == ["dust_heat", "photon_EUV", "L_rec"]
    assert result["discarded"]["dust heat"] == "output"
    header = (tmp_path / "microphysics_func_jac.h").read_text()
    for name in result["output_names"]:
        assert f"#define JACO_HAS_OUTPUT_{name}\n" in header
    assert ("enum OutputIndex { IDX_OUT_dust_heat = 0, IDX_OUT_photon_EUV = 1, IDX_OUT_L_rec = 2, N_OUTPUTS = 3 };"
            in header)
    assert "double L_rec;  /* recombination cooling [erg cm^-3 s^-1] */" in header
    assert "void microphysics_outputs(const SolveVars *vars, const Params *params, Outputs *out);" in header
    src = (tmp_path / "microphysics_outputs.c").read_text()
    assert "out->L_rec = " in src and "out->dust_heat = " in src


def test_outputs_leave_the_solved_system_alone(tmp_path):
    generate(model(), tmp_path / "plain")
    generate(model(outputs=[L_rec, Output("twice_T", 2 * T)]), tmp_path / "outputs")
    for f in ("microphysics_func_jac.c", "jaco_eos.c"):
        assert (tmp_path / "plain" / f).read_text() == (tmp_path / "outputs" / f).read_text()


def test_model_without_outputs_declares_an_empty_set(tmp_path):
    m = Model([CollisionalIonization("H"), GasPhaseRecombination("H+"), pdv], solve_vars=["u", "T", "H+"],
              time_dependent=["T"], species=SPECIES[:3])
    result = generate(m, tmp_path)
    assert result["output_names"] == []
    header = (tmp_path / "microphysics_func_jac.h").read_text()
    assert "enum OutputIndex { N_OUTPUTS = 0 };" in header and "JACO_HAS_OUTPUTS" not in header


def test_heat_of_reads_the_processes_after_the_rules():
    halve = Rule("halve recombination",
                 lambda p: p.transformed(lambda e: e / 2) if p.name.startswith("Gas-phase") else p)
    m = model(outputs=[L_rec], rules=[halve])
    expr = {o.name: o.expr for o in m.network.outputs}["L_rec"]
    assert sp.simplify(expr - GasPhaseRecombination("H+").heat * sp.Rational(-1, 2)) == 0


def test_output_declarations_are_checked():
    with pytest.raises(ValueError, match="C identifier"):
        Output("L NUV")
    with pytest.raises(ValueError, match="not in the model"):
        model(outputs=[Output("L", heat_of=["no such process"])]).network
    with pytest.raises(ValueError, match="conflicting"):
        model(outputs=[Output("L", T)]) + model(outputs=[Output("L", 2 * T)])
    merged = model(outputs=[L_rec]) + Model(outputs=[Output("twice_T", 2 * T)])
    assert [o.name for o in merged.outputs] == ["L_rec", "twice_T"]
    assert Output("a", heat_of=["x", "y"]).heat_of == Output("a", heat_of={"x": 1, "y": 1}).heat_of


def test_output_name_clash_raises(tmp_path):
    with pytest.raises(ValueError, match="share a name"):
        generate(model(outputs=[Output("dust_heat", T)]), tmp_path)


def test_solved_radiation_species_is_no_output(tmp_path):
    m = model().evolve(solve_vars=["u", "T", "H+", "photon_EUV"], time_dependent=["T", "photon_EUV"])
    assert generate(m, tmp_path)["output_names"] == ["dust_heat"]


def test_generated_outputs_match_the_rows(tmp_path):
    """The C outputs at a state equal the rows (and the declared expression) evaluated there by sympy"""
    m = model(outputs=[L_rec])
    result = generate(m, tmp_path)
    driver = tmp_path / "driver.c"
    driver.write_text('#include "microphysics_func_jac.h"\n'
                      "void eval(const double *v, const double *p, double *o) {\n"
                      "    microphysics_outputs((const SolveVars *)v, (const Params *)p, (Outputs *)o);\n}\n")
    so = tmp_path / "liboutputs.so"
    subprocess.run(["gcc", "-O1", "-shared", "-fPIC", f"-I{tmp_path}", "-o", str(so), str(driver),
                    str(tmp_path / "microphysics_outputs.c"), "-lm"], check=True)
    lib = ctypes.CDLL(str(so))
    state = {"u": 0.0, "T": 3000.0, "x_Hplus": 0.2}
    params = {"n_Htot": 50.0, "Td": 20.0, "sigma_c": 3e-9, "x_photon_EUV": 1e-6, "C_2": 1.0, "Delta_t": 1.0,
              "pdv_work": 0.0, "u_initial": 0.0}
    v = np.array([state[n] for n in result["var_names"]])
    p = np.array([params[n] for n in result["param_names"]])
    o = np.zeros(len(result["output_names"]))
    dp = ctypes.POINTER(ctypes.c_double)
    lib.eval(v.ctypes.data_as(dp), p.ctypes.data_as(dp), o.ctypes.data_as(dp))
    xHp, nH = state["x_Hplus"], params["n_Htot"]
    vals = {T: state["T"], sp.Symbol("Td"): params["Td"], sp.Symbol("sigma_c"): params["sigma_c"], n_Htot: nH,
            n_("H"): nH * (1 - xHp), n_("H+"): nH * xHp, n_("e-"): nH * xHp,
            n_("photon_EUV"): params["x_photon_EUV"] * nH, sp.Symbol("C_2"): 1.0}
    rows = m.network
    expected = [float(dict.__getitem__(rows, "dust heat").rhs.subs(vals)),
                float(dict.__getitem__(rows, "photon_EUV").rhs.subs(vals)),
                float(-GasPhaseRecombination("H+").heat.subs(vals))]
    assert o == pytest.approx(expected, rel=1e-12)
    assert o[0] > 0 and o[1] < 0  # the dust gains what the gas loses; the photons are absorbed
