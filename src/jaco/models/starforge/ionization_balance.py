"""Ionization of the neutral gas by cosmic rays and FUV light, and the free electrons it leaves (STARFORGE).

Cosmic rays ionize atomic H at zeta_H, the attenuated rate GIZMO uses for its CR heating and heavy ions, into H+,
which recombines radiatively (network), on grains (WD01) and by charge transfer to Mg. Ionized H2 becomes molecular
ions (H3+, HCO+, ...) that recombine dissociatively (Fromang+02) or pass their charge to Mg. Mg+, also photoionized,
recombines slowly: radiatively and on grains. C+ follows GIZMO's photo/CR ionization balance (Gong+17). Mg stands
for all low-ionization-potential metals and is undepleted, as in GIZMO's heavy-ion cap. The C+, Mg+ and molecular-ion
abundances depend on n_e and through it on each other; their sum y solves y = F(y) by Newton steps in ln y (F falls
with y, so the root is unique). Thermal K and O+ (locked to H+ by charge exchange) are as in GIZMO. Every two-body
rate carries the clumping factor C_2.
"""

import sympy as sp
from jaco.processes import Reaction
from jaco.symbols import n_
from .symbols import T, n_Htot, G_0, x_, cosmicray_ionization_rate_H as zeta
from .grain_assisted_recombination import alpha_grain, GrainAssistedRecombination
from .cosmic_ray_ionization import cosmic_ray_ionization
from .metal_electrons import f_Cplus, x_C_gas, x_e_ions, alkali_electrons, Oplus_electrons

C2 = sp.Symbol("C_2")
# Gas-phase fraction of Mg, held constant: undepleted, as GIZMO's heavy-ion cap assumes.
# TODO: depletion from Jenkins 2009 (2009ApJ...700.1299J), Eq. 10 with the Mg row of Table 4 (0.54 at F* = 0, 0.054
# at F* = 1). Its F*-<n(H)> fit (Sec. 10.2) is in sight-line mean densities: do not feed it local cell densities.
F_GAS_MG = 1
x_Mg = F_GAS_MG * x_("Mg")  # gas-phase Mg per H; the model's x_Mg parameter is the total
K_CT_HPLUS_MG = 1.1e-9  # H+ + Mg -> Mg+ + H (UMIST 2022; Prasad & Huntress 1980)
K_CT_MOLION_MG = 1.0e-9  # H3+ + Mg -> Mg+ + H2 + H (UMIST 2022; Prasad & Huntress 1980)
GAMMA_MG_DRAINE = 6.59e-11  # Mg photoionization in the Draine (1978) field, 1.7 Habing (Heays+17)
alpha_rr_Mgplus = 2.78e-12 * (T / 300) ** -0.68  # Mg+ radiative recombination (UMIST 2022; Pequignot & Aldrovandi 1986)
beta_molion = 3e-6 / sp.sqrt(T)  # dissociative recombination of molecular ions (Fromang+02)
N_NEWTON = 3
BIBLIOGRAPHY = ["2024A&A...682A.109M", "1980ApJS...43....1P", "2017A&A...602A.105H", "1986A&A...161..169P",
                "2002MNRAS.329...18F", "2001ApJ...563..842W", "2017ApJ...843...38G"]


def ion_abundances(x_e):
    """(x_C+, x_Mg+, x_mol+) per H nucleus at electron abundance x_e. The neutral Mg that molecular ions charge-transfer
    to is the Mg left neutral by the other channels, so that their balance stays explicit."""
    n = n_Htot
    xC = x_C_gas * f_Cplus(x_e, n, C2)
    L_Mg = C2 * n * (alpha_rr_Mgplus * x_e + alpha_grain("Mg+", x_e))
    P_Mg0 = GAMMA_MG_DRAINE * G_0 / 1.7 + C2 * K_CT_HPLUS_MG * n * x_("H+")
    x_Mg0 = x_Mg * L_Mg / (L_Mg + P_Mg0)
    x_mol = 2 * zeta * x_("H_2") / (C2 * n * (beta_molion * x_e + K_CT_MOLION_MG * x_Mg0))
    P_Mg = P_Mg0 + C2 * K_CT_MOLION_MG * n * x_mol
    return xC, x_Mg * P_Mg / (P_Mg + L_Mg), x_mol


def solved_electrons(n_steps=N_NEWTON):
    """(intermediates, fixed_electrons, Mg+ symbol): free electrons per H nucleus beyond the solved H and He ions.
    Newton in s = ln y starts from all C and Mg ionized plus the molecular ions with no other electrons."""
    e_other, e_ions = sp.symbols("eIonOther eIonMetal")
    x_e_other = x_e_ions + e_other
    xe = sp.Dummy("x_e")
    lnF = sp.log(sum(ion_abundances(xe)))
    dlnF = sp.diff(lnF, xe)
    inter = [(e_other, alkali_electrons() + Oplus_electrons())]
    s = sp.Symbol("IbS0")
    # (the 1e-40 keeps the derivative of the start finite at x_H2 = 0)
    inter.append((s, sp.log(x_C_gas + x_Mg + sp.sqrt(2 * zeta * x_("H_2") / (C2 * n_Htot * beta_molion) + 1e-40))))
    for k in range(n_steps):
        g, h, s_next = sp.symbols(f"IbG{k} IbH{k} IbS{k + 1}")
        y = sp.exp(s)
        inter += [(g, lnF.subs(xe, x_e_other + y)), (h, 1 - y * dlnF.subs(xe, x_e_other + y)), (s_next, s + (g - s) / h)]
        s = s_next
    x_Mgplus = sp.Symbol("xIbMgplus")
    inter += [(e_ions, sp.exp(s)), (x_Mgplus, ion_abundances(x_e_other + e_ions)[1])]
    return inter, e_other + e_ions, x_Mgplus


def ChargeTransferToMetals(x_Mgplus):
    """H+ + Mg -> H + Mg+ at k n_H+ n_Mg0 C_2; Mg is not a network species, and the Mg+ it makes is in the electron
    balance"""
    return Reaction("H+ -> H", rate=K_CT_HPLUS_MG * n_("H+") * n_Htot * (x_Mg - x_Mgplus), clumping=C2,
                    name="Charge transfer of H+ to Mg", bibliography=["1980ApJS...43....1P"])


def ionization_processes():
    """(processes, intermediates, fixed_electrons) of the explicit H ionization by cosmic rays and its sinks"""
    inter, electrons, x_Mgplus = solved_electrons()
    processes = [cosmic_ray_ionization("H"), GrainAssistedRecombination("H+"), ChargeTransferToMetals(x_Mgplus)]
    return processes, inter, electrons
