"""Solver metadata in the generated header: per-variable time dependence, bounds, scale, charge and start-of-step
parameter, and the budgets of the eliminated abundances, in the layout cooling/jaco_solver.cc reads."""

import subprocess

import pytest
import sympy as sp

from jaco.codegen.gizmo import generate_funcjac_code
from jaco.declarations import Parameter, Species, NO_CEILING
from jaco.model import Model
from jaco.processes import CollisionalIonization, GasPhaseRecombination, Reaction, ThermalTerm
from jaco.symbols import x_

pdv = ThermalTerm(sp.Symbol("pdv_work"), name="PdV work")
k = sp.Symbol("k")
SPECIES = [Species(s) for s in ("H", "H+", "H_2", "He", "He+", "e-")]
PROCESSES = [CollisionalIonization("H"), GasPhaseRecombination("H+"), CollisionalIonization("He"),
             GasPhaseRecombination("He+"), Reaction("2H -> H_2", k, name="H2 formation", bibliography=["t"]),
             Reaction("H_2 -> 2H", k, name="H2 dissociation", bibliography=["t"]), pdv]


def model(time_dependent=("T", "H_2"), species=SPECIES, processes=PROCESSES, **kw):
    return Model(processes, solve_vars=["u", "T", "H+", "He+", "H_2"], time_dependent=time_dependent,
                 species=species, parameters=[Parameter("k")], **kw)


def header_of(m, path):
    result = generate_funcjac_code(m, output_dir=str(path))
    return result, (path / "microphysics_func_jac.h").read_text()


def test_metadata_macros(tmp_path):
    result, h = header_of(model(), tmp_path)
    assert result["var_names"] == ["u", "T", "x_Hplus", "x_Heplus", "x_H_2"]
    assert "#define JACO_HAS_SOLVER_METADATA\n" in h
    assert "#define JACO_VAR_TIME_DEPENDENT_T\n" in h and "#define JACO_VAR_TIME_DEPENDENT_x_H_2\n" in h
    assert "#define JACO_VAR_TIME_DEPENDENT_x_Hplus" not in h
    assert "#define JACO_VAR_TIME_DEPENDENT_INIT {0, 1, 0, 0, 1}\n" in h
    assert "#define JACO_VAR_FLOOR_INIT {0.0, 0.0, 1e-20, 1e-20, 1e-20}\n" in h
    assert "#define JACO_VAR_CEILING_INIT {1e+300, 1e+300, 1.0, 1.0, 1.0}\n" in h
    assert "#define JACO_VAR_CHARGE_INIT {0, 0, 1, 1, 0}\n" in h
    assert "#define JACO_VAR_INITIAL_PARAM_INIT {-1, PARAM_u_initial, -1, -1, PARAM_x_H_2_initial}\n" in h
    assert "#define JACO_N_TD_SPECIES 1\n" in h and "{IDX_x_H_2, PARAM_x_H_2_initial}, \\\n" in h
    assert result["budgets"] == ["x_H", "x_He"] and result["unbounded"] == []
    assert "{1.0, -1, 2, {IDX_x_Hplus, IDX_x_H_2}, {1.0, 2.0}} /* x_H */" in h
    assert "{0.0, PARAM_y, 1, {IDX_x_Heplus}, {1.0}} /* x_He */" in h


def test_time_dependent_set_follows_the_model(tmp_path):
    _, h = header_of(model(time_dependent=("T", "H+", "H_2")), tmp_path)
    assert "#define JACO_VAR_TIME_DEPENDENT_INIT {0, 1, 1, 0, 1}\n" in h
    assert "#define JACO_N_TD_SPECIES 2\n" in h
    assert "{IDX_x_Hplus, PARAM_x_Hplus_initial}, \\\n    {IDX_x_H_2, PARAM_x_H_2_initial}, \\\n" in h


def test_declared_bounds_and_scale(tmp_path):
    species = [s if s.name != "H+" else Species("H+", floor=1e-30, scale=1e-4) for s in SPECIES]
    _, h = header_of(model(species=species), tmp_path)
    assert "#define JACO_VAR_FLOOR_INIT {0.0, 0.0, 1e-30, 1e-20, 1e-20}\n" in h
    assert "#define JACO_VAR_SCALE_INIT {1.0, 1.0, 0.0001, 1.0, 1.0}\n" in h
    assert Species("photon_EUV", "radiation").abundance_ceiling == NO_CEILING
    with pytest.raises(ValueError, match="floor < ceiling"):
        Species("H", floor=1.0, ceiling=0.5)


def test_nonaffine_elimination_has_no_budget(tmp_path):
    """C+ fixed to a nonlinear function of x_H+ leaves x_C = x_C,tot - x_C+ unbounded by the solver"""
    species = SPECIES + [Species("C"), Species("C+")]
    m = model(species=species, fixed={"C+": sp.Symbol("x_C,tot") * x_("H+") ** 2},
              processes=PROCESSES + [ThermalTerm(-1e-30 * sp.Symbol("n_C") * sp.Symbol("n_Htot"), name="C cooling")])
    result, h = header_of(m, tmp_path)
    assert result["budgets"] == ["x_H", "x_He"] and result["unbounded"] == ["x_C"]
    assert "(no solver budget): x_C */" in h


def test_initializers_fill_the_solver_tables(tmp_path):
    """The macros initialize the tables in the layout jaco_solver.cc declares them with"""
    header_of(model(time_dependent=("T", "H+", "H_2")), tmp_path)
    (tmp_path / "check.c").write_text(r'''
#include <stdio.h>
#include <math.h>
#include "microphysics_func_jac.h"
struct Budget { double total; int total_param; int nterm; int k[JACO_BUDGET_MAX_TERMS]; double w[JACO_BUDGET_MAX_TERMS]; };
struct TD { int k, param; };
static const struct Budget budgets[] = JACO_BUDGETS_INIT;
static const struct TD td[] = JACO_TD_SPECIES_INIT;
static const int tdv[N_VARS] = JACO_VAR_TIME_DEPENDENT_INIT;
static const int init[N_VARS] = JACO_VAR_INITIAL_PARAM_INIT;
static const double floor_[N_VARS] = JACO_VAR_FLOOR_INIT, ceil_[N_VARS] = JACO_VAR_CEILING_INIT;
int main(void) {
    printf("%d %g %d %d %d %g %d %d\n", JACO_N_BUDGETS, budgets[0].total, budgets[1].total_param == PARAM_y,
           budgets[0].k[1] == IDX_x_H_2, td[0].param == PARAM_x_Hplus_initial, budgets[0].w[1],
           tdv[IDX_x_Hplus] + tdv[IDX_x_Heplus], init[IDX_T] == PARAM_u_initial);
    printf("%g %g\n", floor_[IDX_x_H_2], ceil_[IDX_T]);
    return 0;
}
''')
    exe = tmp_path / "check"
    subprocess.run(["gcc", "-Wall", "-Werror", f"-I{tmp_path}", "-o", str(exe), str(tmp_path / "check.c")], check=True)
    out = subprocess.run([str(exe)], check=True, capture_output=True, text=True).stdout.split("\n")
    assert out[0] == "2 1 1 1 1 2 1 1" and out[1] == "1e-20 1e+300"
