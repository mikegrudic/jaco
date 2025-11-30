from jaco.models.starforge.h2_chemistry import H2_chemistry
# from jaco.eos import EOS
# import sympy as sp

# for p in h2_chemistry.subprocesses:
#    print(p.equation, p.rate, p.bibliography)
# print(
#    [str(s) for s in radiative_association.radiative_association("H").network["K"].lhs.atoms(sp.Function)]
# )  # PROBLEM: not finding chemical symbol for LHS of this equation
# print(h2_chemistry.network.symbols)
# print(EOS(h2_chemistry.network.chemical_species).internal_energy)

from jaco.models.starforge.starforge import make_model
import sympy as sp
from jaco.symbols import n_

model = make_model()
collected = sp.collect(model.network["H-"].rhs, n_("H-"))
x_Hminus = collected.coeff(n_("H-"), 0) / collected.coeff(n_("H-"), 1)
print(x_Hminus)
# print(H2_chemistry.heat)
