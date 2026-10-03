"""Model-prescribed free electrons (fixed_electrons) and intermediates: charge neutrality uses them instead of the
fixed species' charges, and the generated Jacobian, built by the chain rule through the intermediates, matches
finite differences of the generated RHS, which matches the fully substituted equations."""

import numpy as np
import sympy as sp
from jaco.processes import GasPhaseRecombination, Ionization
from jaco.symbols import x_, n_

T, n_Htot = sp.Symbol("T"), sp.Symbol("n_Htot")
w1, w2 = sp.symbols("wA wB")


def toy_network():
    system = GasPhaseRecombination("H+") + Ionization("H", rate=1e-16 * n_Htot * x_("H") * T**0.1)
    net = system.network
    net.fixed_species = {"C+": sp.Float(3e-4)}
    net.intermediates = [(w1, sp.sqrt(T) * x_("H+") + n_Htot * 1e-6), (w2, sp.exp(-w1) * w1 + 1e-5 * T * x_("H+"))]
    net.fixed_electrons = 1e-3 * w2
    return net


def test_charge_neutrality_uses_fixed_electrons():
    net = toy_network()
    red = net.reduced({"T", "n_Htot"}, [])
    x_e = dict(red.substitutions)[x_("e-")]
    assert x_("C+") not in x_e.free_symbols and sp.simplify(x_e - (x_("H+") + 1e-3 * w2)) == 0


def test_generated_jacobian_through_intermediates():
    net = toy_network()
    result = net.generate_code(["H+"], [], language="python", minimal=False)
    scope = {}
    exec(result["code"], scope)
    f = scope["microphysics_func_jac"]
    params = {"T": 300.0, "n_Htot": 50.0, "C_2": 1.3}
    pvec = [params[p] for p in result["param_names"]]

    def rhs_jac(x):
        rhs, jac = [0.0], [[0.0]]
        f([x], pvec, rhs, jac)
        return rhs[0], jac[0][0]

    # fully substituted reference
    inter = dict(net.intermediates)
    ref = toy_network().reduced({"T", "n_Htot"}, [])["H+"].rhs.subs(w2, inter[w2]).subs(w1, inter[w1])
    ref = sp.lambdify((x_("H+"), T, n_Htot), ref.subs(sp.Symbol("C_2"), 1.3))
    for x in (1e-6, 1e-3, 0.3):
        r, j = rhs_jac(x)
        h = 1e-6 * x
        assert np.isclose(r, ref(x, 300.0, 50.0), rtol=1e-12)
        assert np.isclose(j, (rhs_jac(x + h)[0] - rhs_jac(x - h)[0]) / (2 * h), rtol=1e-6)
