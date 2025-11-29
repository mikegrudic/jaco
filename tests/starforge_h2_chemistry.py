from jaco.models.starforge.h2_chemistry import h2_chemistry, radiative_association
from jaco.eos import EOS
import sympy as sp

# for p in h2_chemistry.subprocesses:
#    print(p.equation, p.rate, p.bibliography)
# print(
#    [str(s) for s in radiative_association.radiative_association("H").network["K"].lhs.atoms(sp.Function)]
# )  # PROBLEM: not finding chemical symbol for LHS of this equation
# print(h2_chemistry.network.symbols)
print(EOS(h2_chemistry.network.chemical_species).density)
