"""Fine-structure cooling (fine_structure.py): the closed-form level populations against a direct solve of the
statistical equilibrium, the low-density and LTE limits, detailed balance with the CMB, the data as read, and critical
densities and cooling curves against the literature."""

import warnings
import numpy as np
import pytest
import sympy as sp
from jaco.symbols import n_, boltzmann_cgs
from ..fine_structure import ATOMS, transition_rates, level_populations, line_power, FineStructureCooling, C2
from ..symbols import T, z

COOLANTS = ("C+", "C", "O", "Si+", "Fe+")
ARGS = (T, n_("H"), n_("H_2"), n_("e-"), n_("H+"), z, C2)
NO_CMB = -0.9  # z: T_cmb = 0.27 K, no CMB excitation of these lines


def numeric(expr):
    f = sp.lambdify(ARGS, expr, "numpy")

    def call(*state):
        with warnings.catch_warnings(), np.errstate(over="ignore"):  # exp(E / T_cmb) overflows to nbar = 0
            warnings.simplefilter("ignore")
            return float(f(*state))
    return call


def state(T_gas, n_H, n_H2=0.0, x_e=0.0, x_Hp=0.0, redshift=0.0, clumping=1.0):
    return (T_gas, n_H, n_H2, x_e * n_H, x_Hp * n_H, redshift, clumping)


def power(coolant, s):
    """Net line power per atom (erg s^-1)"""
    return numeric(sum(line_power(ATOMS[coolant]).values()))(*s)


def rate_matrix(atom, s):
    R = transition_rates(atom)
    return np.array([[numeric(sp.sympify(R[i][j]))(*s) for j in range(len(R))] for i in range(len(R))])


STATES = [state(10.0, 1e-2, 0, 1e-3), state(80.0, 3e2, 1e2, 1e-4, 1e-5, 0.0, 3.0), state(3e3, 1e1, 0, 1e-2, 1e-2),
          state(8e3, 1e6, 1e5, 1e-3, 1e-3, 5.0), state(2.0, 1e2, 1e2, 1e-4, 0, 0.0)]


@pytest.mark.parametrize("coolant", COOLANTS)
@pytest.mark.parametrize("s", STATES)
def test_closed_form_solves_statistical_equilibrium(coolant, s):
    """The spanning-tree populations solve d n_i/dt = sum_j (n_j R_ji - n_i R_ij) = 0, and the net line power equals
    the net energy the collisions take from the gas"""
    atom = ATOMS[coolant]
    R = rate_matrix(atom, s)
    M = R.T - np.diag(R.sum(axis=1))
    M[-1] = 1
    rhs = np.zeros(len(R))
    rhs[-1] = 1
    direct = np.linalg.solve(M, rhs)
    closed = [numeric(f)(*s) for f in level_populations(transition_rates(atom))]
    assert closed == pytest.approx(direct, rel=1e-9, abs=1e-14)  # elimination loses populations below ~1e-15

    collisional = 0.0
    for (collider, u, l), q in atom.rates.items():
        q_ul = numeric(sp.sympify(q))(*s) * s[6] * {"H": s[1], "o-H2": 0.75 * s[2], "p-H2": 0.25 * s[2], "e-": s[3],
                                                     "H+": s[4]}[collider]
        q_lu = q_ul * atom.g[u] / atom.g[l] * np.exp(-(atom.E[u] - atom.E[l]) / s[0])
        collisional += boltzmann_cgs * (atom.E[u] - atom.E[l]) * (closed[l] * q_lu - closed[u] * q_ul)
    assert power(coolant, s) == pytest.approx(collisional, rel=1e-8, abs=1e-40)


@pytest.mark.parametrize("coolant", COOLANTS)
def test_low_density_limit(coolant):
    """Far below the critical densities every collisional excitation from the ground level is radiated: power per atom
    = sum_u n_H q_0u E_u, linear in the clumping factor"""
    atom, T_gas, n_H = ATOMS[coolant], 300.0, 1e-4
    expected = 0.0
    for u in range(1, len(atom.g)):
        q_u0 = numeric(sp.sympify(atom.rates[("H", u, 0)]))(*state(T_gas, n_H))
        expected += n_H * q_u0 * atom.g[u] / atom.g[0] * np.exp(-atom.E[u] / T_gas) * boltzmann_cgs * atom.E[u]
    for clumping in (1.0, 4.0):
        thin = power(coolant, state(T_gas, n_H, redshift=NO_CMB, clumping=clumping))
        assert thin == pytest.approx(clumping * expected, rel=1e-4, abs=0)


@pytest.mark.parametrize("coolant", COOLANTS)
def test_LTE_limit(coolant):
    """Far above them the levels are Boltzmann-populated and the power is sum A E n_u, whatever the clumping"""
    atom, T_gas = ATOMS[coolant], 500.0
    boltzmann = np.array(atom.g) * np.exp(-np.array(atom.E) / T_gas)
    boltzmann /= boltzmann.sum()
    expected = sum(line.A * line.photon_energy * boltzmann[line.upper] for line in atom.lines)
    for clumping in (1.0, 4.0):
        lte = power(coolant, state(T_gas, 1e14, redshift=NO_CMB, clumping=clumping))
        assert lte == pytest.approx(expected, rel=1e-4, abs=0)


@pytest.mark.parametrize("coolant", COOLANTS)
def test_cmb_detailed_balance(coolant):
    """No net power at T = T_cmb, at any density; heating below it"""
    for n_H in (1e-3, 1e2, 1e8):
        lte = power(coolant, state(2.0 * 2.73, n_H, redshift=1.0))  # T = T_cmb(z = 1)
        hot = power(coolant, state(4.0 * 2.73, n_H, redshift=1.0))
        assert abs(lte) < 1e-9 * hot
        assert power(coolant, state(1.5 * 2.73, n_H, redshift=1.0)) < 0


def test_level_energies():
    """Levels from Table 5's transitions to the ground level; its 2->1 energies agree with their differences to the
    table's two significant figures"""
    for coolant, E21 in (("C", 39), ("O", 98)):
        E = ATOMS[coolant].E
        assert E[2] - E[1] == pytest.approx(E21, abs=2)
    assert ATOMS["O"].g == (5, 3, 1) and ATOMS["C"].g == (1, 3, 5) and ATOMS["Fe+"].g == (10, 8)


def test_piecewise_fits_are_continuous():
    """Every piecewise rate joins its neighbours to within 15% (C + H+ 2->1 at 5000 K, as published); the three rows
    data/README.md corrects would jump by e^890 and 250^0.07 = 1.47"""
    for atom in ATOMS.values():
        for key, q in atom.rates.items():
            if isinstance(q, sp.Piecewise):
                f = sp.lambdify(T, q)
                for _, cond in q.args[:-1]:
                    edge = float(cond.rhs)
                    assert f(edge * (1 + 1e-9)) == pytest.approx(f(edge), rel=0.15, abs=0), (atom.name, key, edge)


# Critical densities n_crit,u = sum_l A_ul / sum_l q_ul(H) against the literature: (coolant, level, T, value, source).
# Santoro & Shull 2006 (ApJ 643, 26) Table 1 use Hollenbach & McKee 1989 rates and A = 2.4e-6 (C+) and 2.1e-4 s^-1
# (Si+); Glover & Jappsen take 2.3e-6 and 2.2e-4, and Roueff's (1990) Si+ + H rate, 1.31x lower at 200 K. Draine 2011
# (Table 17.1) uses newer C + H rates (Abrahamsson et al. 2007); Kaufman et al. 1999 (Table 2) quote one figure.
CRITICAL_DENSITIES = [("C+", 1, 200.0, 2.86e3, "SS06"), ("O", 1, 200.0, 6.08e5, "SS06"),
                      ("Si+", 1, 200.0, 2.75e5, "SS06"), ("Fe+", 1, 200.0, 2.24e6, "SS06"),
                      ("C", 1, 100.0, 620.0, "Draine 2011"), ("C", 2, 100.0, 720.0, "Draine 2011"),
                      ("O", 2, 100.0, 1e5, "Kaufman+99")]


@pytest.mark.parametrize("coolant,u,T_gas,value,source", CRITICAL_DENSITIES)
def test_critical_densities(coolant, u, T_gas, value, source):
    atom = ATOMS[coolant]
    A = sum(line.A for line in atom.lines if line.upper == u)
    q = sum(float(sp.sympify(atom.rates[("H", u, l)]).subs(T, T_gas)) for l in range(u))
    assert abs(np.log10(A / q / value)) < np.log10(1.5)


# Cooling curves read off published figures by pixel position (to ~0.01 dex): log10 Lambda (erg cm^-3 s^-1) for 1e-6
# coolant atoms per cm^3 with n_e / n_H = 1e-4 and no CMB. Grassi et al. 2014 (KROME) Fig. 3, n_H = 1: C I from
# Glover & Jappsen's data, as here. Maio et al. 2007 Fig. 3: C+, O, Si+ and Fe+ from Hollenbach & McKee's data; its
# curves sit 0.125 dex below those rates at their stated n_H = 1 for every species and are matched at n_H = 0.75.
# Their Si+ + H and e- rates differ from Glover & Jappsen's and are rescaled by the published rate ratio; their Fe+ has
# five levels, so it is compared at 100 K only, where the a6D 5/2 level adds 1%.
KROME_CI = {12.0: -30.75, 20.0: -30.31, 100.0: -29.47, 1000.0: -29.05, 3000.0: -28.93}
MAIO = {"C+": {100.0: -29.198, 300.0: -28.907, 1000.0: -28.783},
        "O": {100.0: -30.831, 300.0: -29.814, 1000.0: -29.205, 3000.0: -28.816},
        "Si+": {100.0: -29.875, 300.0: -28.732, 1000.0: -28.377, 3000.0: -28.301},
        "Fe+": {100.0: -30.681}}


def si_plus_rate_ratio(T_gas, x_e=1e-4):
    """Glover & Jappsen's Si+ de-excitation rate per H over Hollenbach & McKee's (as Maio et al. tabulate it)"""
    T2 = T_gas / 100
    return (4.95e-10 * T2**0.24 + x_e * 1.2e-6 * T2**-0.5) / (8.0e-10 * T2**-0.07 + x_e * 1.7e-6 * T2**-0.5)


def test_cooling_curve_C_KROME():
    for T_gas, expected in KROME_CI.items():
        Lambda = 1e-6 * power("C", state(T_gas, 1.0, x_e=1e-4, redshift=NO_CMB))
        assert np.log10(Lambda) == pytest.approx(expected, abs=0.03)


@pytest.mark.parametrize("coolant", list(MAIO))
def test_cooling_curves_Maio(coolant):
    for T_gas, expected in MAIO[coolant].items():
        if coolant == "Si+":
            expected += np.log10(si_plus_rate_ratio(T_gas))
        Lambda = 1e-6 * power(coolant, state(T_gas, 0.75, x_e=1e-4, redshift=NO_CMB))
        assert np.log10(Lambda) == pytest.approx(expected, abs=0.03)


def test_process_metadata():
    """Each process names its lines, with tabulated wavelengths and the photon energies its heat uses"""
    oi = FineStructureCooling("O", sp.Symbol("x"))
    assert oi.name == "[OI] fine-structure cooling"
    by_label = {line.label: line for line in oi.lines}
    assert by_label["[OI] 63 um"].wavelength_um == 63.1
    assert by_label["[OI] 63 um"].photon_energy == pytest.approx(230 * boltzmann_cgs, rel=1e-12, abs=0)
    assert {line.label for line in FineStructureCooling("C+", sp.Symbol("x")).lines} == {"[CII] 158 um"}
    assert {line.label for line in FineStructureCooling("Fe+", sp.Symbol("x")).lines} == {"[FeII] 26 um"}
