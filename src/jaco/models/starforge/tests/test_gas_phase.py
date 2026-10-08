"""The gas-phase fractions of gas_phase.py against Jenkins (2009, ApJ 700, 1299) Table 4, and their use by both the
cooling and the electron budget."""

import numpy as np
import pytest
import sympy as sp
from jaco.symbols import x_
from ..gas_phase import GAS_PHASE_FRACTION, F_STAR, JENKINS_2009_TABLE4, jenkins_gas_fraction, x_gas
from .. import ionization_balance as ib
from ..starforge import carbon_abundances
from ..symbols import T, n_Htot, G_0

# Table 4 columns (6) and (7): [X_gas/H] at F* = 0 and F* = 1
JENKINS_2009_DEPLETIONS = {"C": (-0.112, -0.213), "O": (-0.010, -0.236), "Mg": (-0.270, -1.267),
                           "Si": (-0.223, -1.359), "Fe": (-0.951, -2.236)}


@pytest.mark.parametrize("element", list(JENKINS_2009_DEPLETIONS))
def test_eq10_reproduces_table4(element):
    """Eq. 10 from the tabulated A_X, B_X, z_X gives the tabulated depletions at both ends of F*, to their rounding"""
    assert element in JENKINS_2009_TABLE4
    for F_star, expected in zip((0, 1), JENKINS_2009_DEPLETIONS[element]):
        assert np.log10(jenkins_gas_fraction(element, F_star)) == pytest.approx(expected, abs=1e-3)


def test_si_and_fe_at_the_cold_neutral_medium():
    assert F_STAR == 1
    assert GAS_PHASE_FRACTION["Si"] == pytest.approx(10 ** -1.359, rel=3e-3)
    assert GAS_PHASE_FRACTION["Fe"] == pytest.approx(10 ** -2.236, rel=3e-3)


def test_electrons_and_cooling_share_the_gas():
    """The electron budget's C and Mg are the gas-phase ones, and the cooling's C+ and CO partition the same carbon"""
    assert ib.x_C_gas == x_gas("C") and ib.x_Mg == x_gas("Mg")
    x_C_tot = sp.Symbol("x_C,tot")
    vals = {x_C_tot: 2.7e-4, T: 50.0, G_0: 1.0, x_("H_2"): 0.5}
    gas = float(x_gas("C").subs(vals))
    assert gas == pytest.approx(GAS_PHASE_FRACTION["C"] * 2.7e-4, rel=1e-12)
    for n in (1e-2, 1e3, 1e6):  # all C+ to all CO
        ab = carbon_abundances()
        total = float((ab["C+"] + ab["CO"]).subs(vals).subs(n_Htot, n))
        assert total == pytest.approx(gas, rel=1e-3)
