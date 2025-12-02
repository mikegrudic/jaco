"""Routines for generating a cooling solver function for GIZMO"""

from jaco.models.wind_comparison import cooling, heating
from sympy.utilities.codegen import C99CodeGen
from sympy.codegen.ast import Assignment
import sympy as sp


def generate_funcjac_code(system, cse=True):
    """Generates a .h and .c pair of C source files specifying the function funcjac"""

    system = cooling + heating
    solve_vars = ["u", "T"]
    time_dependent = ["T"]
    func, jac, indices = system.solver_functions(solve_vars, time_dependent, return_jac=True)
    funcjac = sp.Matrix(func + sp.flatten(jac))
    X = sp.MatrixSymbol("X", len(indices), 1)

    index_defs = []
    for var in indices:
        funcjac = funcjac.subs(
            var, X[indices[var]]
        )  # .subs("T",sp.MatrixSymbol("X",2,1)[1]).subs("u",sp.MatrixSymbol("X",2,1)[0])
        index_defs.append(f"#define IDX_{var} {indices[var]}")
    paramsvars = funcjac.free_symbols.copy()
    paramsvars.remove(X)
    P = sp.MatrixSymbol("params", len(paramsvars), 1)
    assignments = []
    for i, p in enumerate(paramsvars):
        funcjac = funcjac.subs(p, P[i])
        assignments.append(Assignment(P[i], p))
        index_defs.append(f"#define IDX_{p} {i}")
    cg = C99CodeGen(cse=cse)
    routine = cg.routine("microphysics_func_jac", funcjac)  # cg.routine("func",sp.Matrix(func + sp.flatten(jac)))
    cg.write([routine], "microphysics_func_jac", to_files=True)

    with open("indices.h", "w") as F:
        for i in index_defs:
            F.write(i + "\n")
        F.write(f"#define NUM_VARS {len(indices)}")

    with open("assignments.h", "w") as F:
        F.write(f"double params[{len(paramsvars)}];\n")
        F.write(sp.ccode(Assignment(P, sp.Matrix(list(paramsvars)))) + "\n")
        F.write(f"double X[{len(indices)}];\n")
        F.write(sp.ccode(Assignment(X, sp.Matrix(list(indices)))))


if __name__ == "__main__":
    generate_funcjac_code(cooling + heating, cse=True)
