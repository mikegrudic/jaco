"""GIZMO's optically-thick limit on the cooling (cooling.cc 1438-1483, COOL_LOW_TEMPERATURES without
COOL_LOWTEMP_THIN_ONLY), for models without a transported IR band.

GIZMO treats a cell as a slab at temperature T with the Sobolev/TreeCol column Sigma to infinity, radiating from a
photosphere (Rafikov 2007, eta = 1): per unit volume at most

    L = sigma_SB T^4 A_eff nHcgs / (1 + tau),   A_eff = 2.3 m_p / Sigma,   tau = kappa Sigma,

with nHcgs = 0.76 rho / m_p, and caps the net heating or cooling of CoolingRate, |Q| <= L. Q is every term CoolingRate
sums, the hydrodynamic PdV work excluded. Here Q is the heat of every process of the model but PdV work and those named
exempt: the H2/H- reactions, whose chemical heat GIZMO's CoolingRate does not have (capping it would let the gas
dissociate H2 without paying the binding energy). The energy H2 formation radiates is a separate cooling term
(h2_chemistry/chemical_heat.py), capped. Collisional ionization is capped, as GIZMO caps it.

The limit is a rule that scales the energy rows (heat, and the dust heat of gas-dust collisions) of every capped
process by f, so their sum Q becomes Q f. Q and f are intermediates (Q_thick, f_thick), evaluated once per call and
differentiated by the chain rule, so neither the sum nor the cap is expanded into every term.

Deviations from GIZMO:
- the hard cap min(|Q|, L) is the smooth Q f, f = (1 + (Q/L)^P)^(-1/P), P = 16, so that Newton sees a differentiable
  rate: 4.2% under the hard cap at |Q| = L, under 0.1% beyond a factor 1.5 either side (smooth_cap_factor);
- kappa is a parameter, GIZMO's kappa_eff (cooling.cc 1462-1479) evaluated by the host at the cell's temperature at the
  start of the step; GIZMO evaluates it at the temperature of each iterate. L's T^4 is implicit here as in GIZMO.
"""

import sympy as sp
from jaco.declarations import Parameter
from jaco.model import Rule
from jaco.process import Process
from .symbols import T, NH, nH_gizmo_cooling

SIGMA_SB = 5.67e-5  # as GIZMO writes the slab emission
PROTONMASS = 1.67262178e-24  # GIZMO's PROTONMASS_CGS, with which the host converts Sigma to N_H
SMOOTHNESS = 16  # P: the exponent of the smooth minimum
LIMIT_ID = "GIZMO's optically-thick cooling limit"
ENERGY_ROWS = ("heat", "dust heat")  # the gas's energy and the dust reservoir gas-dust collisions exchange it with

kappa_thick = sp.Symbol("kappa_thick")
Q_THICK = sp.Symbol("Q_thick")  # the capped processes' summed heat, uncapped (an intermediate)
F_THICK = sp.Symbol("f_thick")  # the factor the limit puts on it (an intermediate)
KAPPA_THICK = Parameter("kappa_thick", "cm^2 g^-1", 0.1,
                        "opacity of GIZMO's optically-thick cooling limit at the cell's temperature (cooling.cc 1462-1479)")


def thick_limit():
    """L, GIZMO's Lambda_Thick_BlackBody times nHcgs^2 [erg cm^-3 s^-1]: sigma T^4 A_eff nHcgs / (1 + tau), with
    Sigma = N_H m_p, so A_eff nHcgs = 2.3 nHcgs / N_H"""
    return SIGMA_SB * T**4 * 2.3 * nH_gizmo_cooling / (NH * (1 + kappa_thick * NH * PROTONMASS))


def smooth_cap_factor(Q, L, p=SMOOTHNESS):
    """f = (1 + r^p)^(-1/p), r = |Q| / L, so that Q f is a smooth min(|Q|, L) sign(Q). Written as
    max(1, r)^-1 (1 + min(r, 1/r)^p)^(-1/p), the same function, so that r^p cannot overflow. Q should be a symbol:
    Max and Min of the summed heat itself make sympy query its sign, which does not finish."""
    s = (Q / L) ** 2
    return 1 / sp.sqrt(sp.Max(1, s)) * (1 + sp.Min(s, 1 / s) ** (p // 2)) ** sp.Rational(-1, p)


def capped_ids(model, exempt=()):
    """Ids of the processes the limit caps: every process with a heat row but PdV work and those in exempt"""
    skip = set(exempt) | {"PdV work"}
    return [i for i, p in model.processes.items() if i not in skip and "heat" in p.network]


def _scale_energy_rows(p):
    rows = {k: e.rhs for k, e in p.network.items()}
    return Process(p.name, p.bibliography, {k: (v * F_THICK if k in ENERGY_ROWS else v) for k, v in rows.items()})


def with_optically_thick_cap(model, exempt=()):
    """model with GIZMO's optically-thick limit on the summed heat of capped_ids(model, exempt), as the model assembles
    them (rules applied). Apply it last: processes added afterwards are scaled without counting towards Q."""
    unknown = set(exempt) - set(model.processes)
    if unknown:
        raise ValueError(f"the optically-thick limit exempts processes not in the model: {sorted(unknown)}")
    ids = capped_ids(model, exempt)
    Q = sp.Add(*[p.heat for p in model.subprocesses if p.name in set(ids)])
    rule = Rule(LIMIT_ID, _scale_energy_rows, only=ids)
    return model.evolve(rules=(*model.rules, rule), parameters=(*model.parameters, KAPPA_THICK),
                        intermediates=(*model.intermediates, (Q_THICK, Q),
                                       (F_THICK, smooth_cap_factor(Q_THICK, thick_limit()))))
