"""Fine-structure line cooling of C+, C, O, Si+ and Fe+ (STARFORGE): optically thin lines from levels in statistical
equilibrium.

Each coolant's lowest two or three levels are populated by collisions (with H, ortho- and para-H2, e- and H+ where the
data have rates), spontaneous decay, and absorption and stimulated emission of the CMB at T_cmb = 2.73 (1 + z) K. The
populations are the closed-form solution of the 2- or 3-level statistical equilibrium, so they stay differentiable for
Newton. The heat is the net power of the lines, the sum over u > l of
E_ul A_ul [n_u (1 + nbar_ul) - n_l (g_u / g_l) nbar_ul], nbar the CMB photon occupation number: it vanishes at
T = T_cmb, in place of the cmb_bath_factor approximation of the other low-temperature coolants. The collisional rates
carry the clumping factor C_2, so the emission is clumped below the critical densities and not above them, as the H2
cooling is.

Data (data/README.md): Glover & Jappsen 2007 (GJ07) Tables 5 and 6 for C, C+, O and Si+; Hollenbach & McKee 1989 for
Fe+. The level energies are GJ07's energies of the transitions to the ground level, so the 2->1 energies are their
differences (O: 100 K, printed as 98 K) and detailed balance holds exactly. Ortho- and para-H2 are in the ratio 3:1, as
the model's H2 cooling rates (Glover & Abel 2008) assume.

Levels: C+ and Si+ have a single excited level below 5 eV; C and O are their ground 3P triplets (the 1D2 terms lie at
14,700 and 22,800 K). Fe+ is its 26 um line alone: Fe is the most depleted of the five (gas_phase.py), so that even with
its a6D 5/2..1/2 levels, which would at most double it below 1e4 K, Fe+ cooling stays below 1% of O's in neutral gas.

Abundances (gas_phase.py): C+ and neutral C are the model's gas-phase carbon partition. O, Si+ and Fe+ are in neutral
gas only: O at its gas-phase abundance outside CO, ionized where H is (charge exchange with H+); Si and Fe singly
ionized (ionization potentials below 13.6 eV), doubly where H is ionized.

The terms carry lowtemp_truncation: above ~3e4 K the tabulated metal-line cooling holds these elements, and the C+
abundance (Tielens' interpolation) has no collisional ionization.

TODO: escape-probability (LVG) trapping of [CII] 158 um and [OI] 63 um, whose optical depths reach unity in dense
PDR-like gas; the lines are optically thin here.
"""

from dataclasses import dataclass
from importlib.resources import files
import numpy as np
import sympy as sp
from jaco.processes import ThermalTerm
from jaco.symbols import n_, x_, boltzmann_cgs
from .symbols import T, T_cmb, n_Htot, lowtemp_truncation
from .gas_phase import x_gas

C2 = sp.Symbol("C_2")
H2_ORTHO_FRACTION = sp.Rational(3, 4)  # ortho:para = 3:1
# The quartic-in-ln T fits (form lnpoly, C + e-) run away above this (x2 by 3e4 K, x1e3 by 1e5 K); held at its value
LNPOLY_T_MAX = 2.0e4
DATA = files("jaco.models.starforge").joinpath("data")
ATOMIC_DATA = ("glover_jappsen_2007_table5.txt", "hollenbach_mckee_1989_fe_ii_atomic.txt")
RATE_DATA = ("glover_jappsen_2007_table6.txt", "hollenbach_mckee_1989_fe_ii_rates.txt")
GJ07, HM89 = "2007ApJ...666....1G", "1989ApJ...342..306H"
BIBLIOGRAPHY = {"C+": [GJ07], "C": [GJ07], "O": [GJ07], "Si+": [GJ07],
                "Fe+": [HM89, "2007MNRAS.379..963M", "2014MNRAS.439.2386G"]}
SPECTROSCOPIC = {"C+": "[CII]", "C": "[CI]", "O": "[OI]", "Si": "[SiI]", "Si+": "[SiII]", "Fe+": "[FeII]"}
COLLIDER_DENSITY = {"H": n_("H"), "o-H2": H2_ORTHO_FRACTION * n_("H_2"), "p-H2": (1 - H2_ORTHO_FRACTION) * n_("H_2"),
                    "e-": n_("e-"), "H+": n_("H+")}


@dataclass(frozen=True)
class Line:
    """A radiative transition u -> l: its label, wavelength (um, as tabulated) and photon energy (erg)"""
    label: str
    upper: int
    lower: int
    wavelength_um: float
    photon_energy: float
    A: float


@dataclass(frozen=True)
class Atom:
    """Levels (statistical weight, energy in K above the ground level), radiative lines and collisional de-excitation
    rate coefficients {(collider, u, l): expression in T}"""
    name: str
    g: tuple
    E: tuple
    lines: tuple
    rates: dict


def _records(filename):
    for row in DATA.joinpath(filename).read_text().splitlines():
        if row.strip() and not row.lstrip().startswith("#"):
            yield row.split()


def _transition(s):
    u, l = s.split("->")
    return int(u), int(l)


def _fit(form, c):
    """Rate coefficient (cm^3 s^-1) of a fit form of the data files (header of glover_jappsen_2007_table6.txt)"""
    match form:
        case "const":
            return c[0]
        case "lin":
            return c[0] + c[1] * T
        case "pow":
            return c[0] * T ** c[1]
        case "pow2":
            return c[0] * (T / 100) ** c[1]
        case "powT0":
            return c[0] * (T / c[2]) ** c[1]
        case "exp1":
            return c[0] + c[1] * sp.exp(-T / c[2])
        case "exp2":
            return c[0] + c[1] * sp.exp(-T / c[3]) + c[2] * sp.exp(-2 * T / c[3])
        case "quadpow":
            return (c[0] + c[1] * T + c[2] * T**2) * T ** c[3]
        case "lnpoly":
            T_fit = sp.Min(T, LNPOLY_T_MAX)
            lnT = sp.log(T_fit)
            # Horner form: written as a sum, sympy turns exp(c_k ln T) into T^c_k, which overflows (T^228)
            poly = c[-1]
            for ck in reversed(c[1:-1]):
                poly = ck + lnT * poly
            return c[0] / sp.sqrt(T_fit) * sp.exp(poly)
    raise ValueError(f"unknown fit form {form}")


def _upper_edge(T_range):
    """Upper temperature edge of a T-range of the data files (inf for the last piece)"""
    if T_range == "all" or T_range.startswith("T>"):
        return np.inf
    return float(T_range.split("<=")[-1])


def _read_rates():
    pieces = {}
    for filename in RATE_DATA:
        for coolant, collider, transition, T_range, _refs, form, *coeffs in _records(filename):
            key = (coolant, collider, *_transition(transition))
            pieces.setdefault(key, []).append((_upper_edge(T_range), _fit(form, [float(c) for c in coeffs])))
    rates = {}
    for key, fits in pieces.items():
        fits.sort(key=lambda p: p[0])
        rates[key] = fits[0][1] if len(fits) == 1 else sp.Piecewise(
            *[(q, T <= edge) for edge, q in fits[:-1]], (fits[-1][1], True))
    return rates


def _read_atoms():
    rows = {}
    for filename in ATOMIC_DATA:
        for coolant, transition, g_u, g_l, wavelength, E_K, A in _records(filename):
            rows.setdefault(coolant, {})[_transition(transition)] = (int(g_u), int(g_l), float(wavelength), float(E_K),
                                                                      float(A))
    rates = _read_rates()
    atoms = {}
    for coolant, transitions in rows.items():
        n_levels = 1 + max(u for u, _ in transitions)
        g = [transitions[(1, 0)][1]] + [transitions[(u, 0)][0] for u in range(1, n_levels)]
        E = [0.0] + [transitions[(u, 0)][3] for u in range(1, n_levels)]
        lines = tuple(Line(f"{SPECTROSCOPIC[coolant]} {wl:.0f} um", u, l, wl, boltzmann_cgs * (E[u] - E[l]), A)
                      for (u, l), (_, _, wl, _, A) in sorted(transitions.items()))
        atom_rates = {(c, u, l): q for (cl, c, u, l), q in rates.items() if cl == coolant}
        atoms[coolant] = Atom(coolant, tuple(g), tuple(E), lines, atom_rates)
    return atoms


ATOMS = _read_atoms()


def cmb_occupation(E_K):
    """CMB photon occupation number at a transition of energy E_K (K)"""
    return 1 / (sp.exp(E_K / T_cmb) - 1)


def transition_rates(atom, clumping=C2):
    """R[i][j]: rate (s^-1) per atom in level i of its transitions to level j"""
    n = len(atom.g)
    R = [[sp.S.Zero] * n for _ in range(n)]
    for (collider, u, l), q in atom.rates.items():
        C_ul = clumping * COLLIDER_DENSITY[collider] * q
        R[u][l] += C_ul
        R[l][u] += C_ul * sp.Rational(atom.g[u], atom.g[l]) * sp.exp(-(atom.E[u] - atom.E[l]) / T)
    for line in atom.lines:
        u, l, nbar = line.upper, line.lower, cmb_occupation(atom.E[line.upper] - atom.E[line.lower])
        R[u][l] += line.A * (1 + nbar)
        R[l][u] += sp.Rational(atom.g[u], atom.g[l]) * line.A * nbar
    return R


def level_populations(R):
    """Fractions in each level in statistical equilibrium of a 2- or 3-level system with transition rates R[i][j]: the
    weight of level k is the sum over the spanning trees directed into k of the products of their rates"""
    if len(R) == 2:
        w = [R[1][0], R[0][1]]
    elif len(R) == 3:
        w = [R[1][0] * R[2][0] + R[1][0] * R[2][1] + R[1][2] * R[2][0],
             R[0][1] * R[2][1] + R[0][1] * R[2][0] + R[0][2] * R[2][1],
             R[0][2] * R[1][2] + R[0][2] * R[1][0] + R[0][1] * R[1][2]]
    else:
        raise NotImplementedError("closed-form populations for 2 or 3 levels only")
    total = sum(w)
    return [wk / total for wk in w]


def line_power(atom, clumping=C2):
    """{line label: net power emitted in the line per atom (erg s^-1)}"""
    f = level_populations(transition_rates(atom, clumping))
    power = {}
    for line in atom.lines:
        u, l = line.upper, line.lower
        nbar = cmb_occupation(atom.E[u] - atom.E[l])
        g_ratio = sp.Rational(atom.g[u], atom.g[l])
        power[line.label] = line.photon_energy * line.A * (f[u] * (1 + nbar) - f[l] * g_ratio * nbar)
    return power


class FineStructureCooling(ThermalTerm):
    """Net fine-structure line cooling of one coolant at the given abundance per H nucleus; lines: its transitions
    with their wavelengths and photon energies, for a spectral declaration to route"""

    def __init__(self, coolant, abundance, clumping=C2):
        self.atom = ATOMS[coolant]
        self.abundance = abundance
        self.collisional_clumping = clumping
        self.lines = self.atom.lines
        heat = -n_Htot * abundance * sum(line_power(self.atom, clumping).values()) * lowtemp_truncation
        super().__init__(heat, name=f"{SPECTROSCOPIC[coolant]} fine-structure cooling",
                         bibliography=BIBLIOGRAPHY[coolant])

    def unclumped(self):
        return FineStructureCooling(self.atom.name, self.abundance, clumping=1)


def neutral_gas(x):
    """x in neutral gas only: times the fraction of H nuclei not ionized"""
    return x * (1 - x_("H+"))


def coolant_abundances():
    """{coolant: abundance per H nucleus}"""
    x_C_neutral = x_gas("C") - x_("C+") - x_("CO")
    return {
        "C+": x_("C+"),
        "C": x_C_neutral,
        "O": neutral_gas(x_gas("O") - x_("CO")),
        "Si+": neutral_gas(x_gas("Si")),
        "Fe+": neutral_gas(x_gas("Fe")),
    }


def fine_structure_cooling():
    """The fine-structure cooling processes of every coolant"""
    return [FineStructureCooling(c, x) for c, x in coolant_abundances().items()]
