"""Radiation bands as solved species, per-row rate factors, energy transfers between rows, algebraic variables (a dust
temperature determined by its row's steady state), rules restricted to named processes, and outputs that read a
process's row: the framework pieces a model of the matter-radiation coupling is built from."""

import subprocess

import pytest
import sympy as sp

from jaco.codegen.gizmo import generate_funcjac_code
from jaco.declarations import Parameter, Species, Variable, Output, NO_CEILING
from jaco.model import Model, Rule
from jaco.processes import (CollisionalIonization, GasPhaseRecombination, Reaction, ThermalTerm, Transfer)
from jaco.symbols import n_, x_

pdv = ThermalTerm(sp.Symbol("pdv_work"), name="PdV work")
sigma, c_ratio, hnu, Td = sp.Symbol("sigma"), sp.Symbol("c_ratio"), sp.Symbol("hnu"), sp.Symbol("Td")
C_LIGHT, EV = 2.9979e10, 1.602176634e-12
SPECIES = [Species(s) for s in ("H", "H+", "e-")]
PHOTONS = [Species("photon_EUV", "radiation", "ionizing photons per H nucleus", floor=0.0),
           Species("photon_ONIR", "radiation", "optical energy per H nucleus [eV]", floor=0.0),
           Species("dust heat", "energy", "energy the gas gives the dust")]


def photoionization(row_factors=None):
    return Reaction("H + photon_EUV -> H+ + e- + photon_ONIR", C_LIGHT * sigma, 3e-12, name="photoionization",
                    row_factors=row_factors if row_factors is not None else
                    {"photon_EUV": c_ratio, "photon_ONIR": c_ratio * hnu}, bibliography=["test"])


def dust_processes():
    """Gas-dust collisions into the dust reservoir, dust emission out of it into the optical band, and dust
    absorption of the optical band"""
    gas_dust = ThermalTerm(1e-33 * n_("Htot") ** 2 * sp.sqrt(sp.Symbol("T")) * (Td - sp.Symbol("T")),
                           name="gas-dust", reservoir="dust heat")
    emission = Transfer(1e-30 * n_("Htot") * Td ** 6, {"dust heat": -1, "photon_ONIR": c_ratio / EV}, name="emission")
    absorption = Transfer(C_LIGHT * 1e-21 * n_("Htot") * n_("photon_ONIR") * EV,
                          {"photon_ONIR": -c_ratio / EV, "dust heat": 1}, name="absorption")
    return [gas_dust, emission, absorption]


def model(time_dependent=("T", "H+", "photon_EUV", "photon_ONIR"), **kw):
    procs = [CollisionalIonization("H"), GasPhaseRecombination("H+"), photoionization(), pdv] + dust_processes()
    return Model(procs, solve_vars=["u", "T", "H+", "photon_EUV", "photon_ONIR", "Td"], time_dependent=time_dependent,
                 species=SPECIES + PHOTONS, parameters=[Parameter(p) for p in ("sigma", "c_ratio", "hnu")],
                 variables=[Variable("Td", "dust heat", floor=2.73, ceiling=1e4, units="K", doc="dust temperature")],
                 **kw)


def test_row_factors_scale_rows():
    """A band transported at c_tilde loses photons at c_tilde/c of the matter rate; a band counted in energy gains the
    energy per event"""
    p = photoionization().with_radiation({"photon_EUV", "photon_ONIR"})  # as a model declaring them radiation
    rate = C_LIGHT * sigma * n_("H") * n_("photon_EUV")
    rows = {k: e.rhs for k, e in p.network.items()}
    assert sp.simplify(rows["H+"] - rate) == 0 and sp.simplify(rows["H"] + rate) == 0
    assert sp.simplify(rows["photon_EUV"] + c_ratio * rate) == 0
    assert sp.simplify(rows["photon_ONIR"] - c_ratio * hnu * rate) == 0
    assert sp.simplify(p.heat - 3e-12 * rate) == 0  # the heat follows the matter rate
    assert sp.simplify(photoionization().network["photon_EUV"].rhs + sp.Symbol("C_2") * c_ratio * rate) == 0
    with pytest.raises(ValueError, match="not in the equation"):
        Reaction("H -> H+ + e-", 1.0, row_factors={"photon_EUV": 0.5}, bibliography=["t"])


def test_energy_weighted_photon_stoichiometry():
    """Non-integer photon consumption as a factor on the photon row: e.g. H2 photodissociation taking N photons of
    mean energy E from a band counted in eV"""
    N, E = 6.5, 12.4
    p = Reaction("H_2 + photon_FUV -> 2H", rate=1e-11 * n_("H_2"), row_factors={"photon_FUV": c_ratio * N * E},
                 bibliography=["t"])
    rows = {k: e.rhs for k, e in p.network.items()}
    assert sp.simplify(rows["photon_FUV"] + c_ratio * N * E * 1e-11 * n_("H_2")) == 0
    assert sp.simplify(rows["H"] - 2e-11 * n_("H_2")) == 0


def test_transfer_and_with_rows():
    t = Transfer(sp.Symbol("P"), {"photon_ONIR": -2.0, "dust heat": 1}, name="t")
    rows = {k: e.rhs for k, e in t.network.items()}
    assert rows == {"heat": 0, "photon_ONIR": -2.0 * sp.Symbol("P"), "dust heat": sp.Symbol("P")}
    line = ThermalTerm(-sp.Symbol("L"), name="line")
    routed = line.with_rows({"photon_ONIR": -line.heat / EV})
    assert routed.network["photon_ONIR"].rhs == sp.Symbol("L") / EV and routed.heat == -sp.Symbol("L")
    assert "photon_ONIR" not in line.network  # the original is unchanged


def test_rule_only():
    seen = []

    def tag(p):
        seen.append(p.name)
        return p
    m = model(rules=[Rule("tag", tag, only={"photoionization", "gas-dust"})])
    m.network
    assert sorted(seen) == ["gas-dust", "photoionization"]
    with pytest.raises(ValueError, match="applies to processes not in the model"):
        model(rules=[Rule("tag", tag, only={"nope"})]).network


def test_variable_declarations():
    with pytest.raises(ValueError, match="not solve variables"):
        Model([pdv], solve_vars=["u", "T"], species=SPECIES, variables=[Variable("Td", "dust heat")])
    with pytest.raises(ValueError, match="steady state"):
        model(time_dependent=("T", "Td"))
    with pytest.raises(ValueError, match="kind"):
        Variable("Td", "dust heat", kind="pressure")
    m = Model([pdv, ThermalTerm(sp.Symbol("q"), name="q")], solve_vars=["u", "T", "Td"], time_dependent=("T",),
              species=SPECIES, variables=[Variable("Td", "dust heat")], parameters=[Parameter("q")])
    with pytest.raises(ValueError, match="no process has the row"):
        m.network


def test_output_reads_a_row():
    """heat_of keys (process, row) take that row of the process as assembled"""
    o = Output("onir_from_dust", heat_of={("emission", "photon_ONIR"): 1, "photoionization": 2})
    m = model(outputs=[o])
    resolved = {r.name: r.expr for r in m.network.outputs}["onir_from_dust"]
    procs = {p.name: p for p in m.subprocesses}
    expected = procs["emission"].network["photon_ONIR"].rhs + 2 * procs["photoionization"].heat
    assert sp.simplify(resolved - expected) == 0


def test_codegen_metadata(tmp_path):
    """Bands are time-dependent solve variables with floor 0 and no ceiling; the dust temperature is a variable of its
    own kind, bounded as declared, whose equation is the dust reservoir's row (not an output any more)"""
    result = generate_funcjac_code(model(), output_dir=str(tmp_path), source_ext=".cc")
    h = (tmp_path / "microphysics_func_jac.h").read_text()
    assert result["var_names"] == ["u", "T", "x_Hplus", "x_photon_EUV", "x_photon_ONIR", "Td"]
    assert "#define JACO_VAR_TIME_DEPENDENT_INIT {0, 1, 1, 1, 1, 0}\n" in h
    assert "#define JACO_VAR_FLOOR_INIT {0.0, 0.0, 1e-20, 0.0, 0.0, 2.73}\n" in h
    assert f"#define JACO_VAR_CEILING_INIT {{{NO_CEILING!r}, {NO_CEILING!r}, 1.0, {NO_CEILING!r}, {NO_CEILING!r}, 10000.0}}\n" in h
    assert ("#define JACO_VAR_KIND_INIT {JACO_KIND_ENERGY, JACO_KIND_GAS_TEMPERATURE, JACO_KIND_SPECIES, "
            "JACO_KIND_RADIATION, JACO_KIND_RADIATION, JACO_KIND_TEMPERATURE}\n") in h
    assert "#define JACO_VAR_CHARGE_INIT {0, 0, 1, 0, 0, 0}\n" in h
    assert "{IDX_x_photon_EUV, PARAM_x_photon_EUV_initial}" in h
    assert "JACO_HAS_PARAM_Td" not in h and "#define JACO_HAS_VAR_Td\n" in h
    assert "dust_heat" not in result["output_names"]
    src = (tmp_path / "microphysics_func_jac.cc").read_text()
    assert "rhs->Td = " in src and "rhs->x_photon_EUV = " in src


def test_generated_rows_conserve_photons(tmp_path):
    """Compile the generated code and evaluate it: at any state the EUV row is -c_ratio times the H+ row's
    photoionization part, the ONIR row gets hnu per event, and the Td row is the dust's energy balance"""
    generate_funcjac_code(model(), output_dir=str(tmp_path), source_ext=".cc")
    (tmp_path / "main.cc").write_text(r'''
#include <stdio.h>
#include <math.h>
#include "microphysics_func_jac.h"
int main() {
    SolveVars v = {}, F; Params p = {}; double J[N_VARS][N_VARS];
    p.n_Htot = 100; p.sigma = 6e-18; p.c_ratio = 1e-3; p.hnu = 20; p.Delta_t = 1e10;
    v.T = 1e4; v.x_Hplus = 0.3; v.x_photon_EUV = 0.02; v.x_photon_ONIR = 0.5; v.Td = 30;
    p.x_Hplus_initial = v.x_Hplus; p.x_photon_EUV_initial = v.x_photon_EUV; p.x_photon_ONIR_initial = v.x_photon_ONIR;
    microphysics_func_jac(&v, &p, &F, J);
    printf("%.17g %.17g %.17g %.17g\n", F.x_photon_EUV, F.x_photon_ONIR, F.Td, J[IDX_Td][IDX_Td]);
    return 0;
}
''')
    srcs = ["main.cc", "microphysics_func_jac.cc"]
    subprocess.run(["g++", "-O1", "-I.", *srcs, "-o", "main"], cwd=tmp_path, check=True)
    out = subprocess.run(["./main"], cwd=tmp_path, check=True, capture_output=True, text=True).stdout.split()
    F_euv, F_onir, F_td, J_tdtd = map(float, out)
    nH, ion = 100.0, C_LIGHT * 6e-18 * 100 * 0.7 * 100 * 0.02  # n_H0 n_gamma c sigma
    assert F_euv == pytest.approx(-1e-3 * ion, rel=1e-12)
    absorbed = C_LIGHT * 1e-21 * nH * nH * 0.5 * EV
    emitted = 1e-30 * nH * 30.0**6
    assert F_onir == pytest.approx(1e-3 * 20 * ion + 1e-3 * (emitted - absorbed) / EV, rel=1e-10)
    assert F_td == pytest.approx(1e-33 * nH**2 * 100 * (1e4 - 30) - emitted + absorbed, rel=1e-10)
    assert J_tdtd == pytest.approx(-1e-33 * nH**2 * 100 - 6e-30 * nH * 30.0**5, rel=1e-10)
