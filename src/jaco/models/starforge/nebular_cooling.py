"""Nebular forbidden-line cooling of photoionized gas (O+, O++, N+, S+, Ne+).

Kim, Gong, Kim & Ostriker 2023 (ApJS 264, 10), Eq. 47, as implemented in GIZMO's CoolingRate()
under RT_CHEM_PHOTOION && METALS. The tabulated metal-line cooling assumes collisional ionization and
under-predicts forbidden-line cooling in photoionized gas; this term supplies it.

Volumetric rate: f_neb * Z_d * C(T, n_e) * n_e * n_H+, with GIZMO's metallicity factor
Metallicity[0]/SolarAbundances[0] carried by Z_d, which GIZMO's interface sets to exactly that ratio.

f_neb is a runtime switch (1 = on, 0 = off) so one generated model serves builds with and without
explicit photoionization; GIZMO applies this term only under RT_CHEM_PHOTOION.
"""

import sympy as sp
from jaco.math import logistic
from jaco.symbols import n_
from jaco.processes import NBodyProcess
from .symbols import T, Z_dust

f_neb = sp.Symbol("f_neb")

# GIZMO evaluates the fit only for 2e3 < T < 5e4 K. The fit is clamped to that range and the lower
# edge replaced by a logistic step 0.01 dex wide; the upper edge is already ~1e-6 of the peak after
# the CIE taper, so the clamp alone suffices there.
T_NEB_MIN, T_NEB_MAX = 2.0e3, 5.0e4


def nebular_cooling_coefficient(T=T, n_e=n_("e-")):
    """Kim+23 Eq. 47 cooling coefficient in erg cm^3 s^-1 (per n_e n_H+), at solar metallicity."""
    T4 = sp.Min(sp.Max(T, T_NEB_MIN), T_NEB_MAX) / 1.0e4
    lnT4 = sp.log(T4)
    log10_fneb = 0.692 + lnT4 * (-0.586 + lnT4 * (0.816 + lnT4 * (-0.505 + lnT4 * (0.118 + lnT4 * (0.00766 - 0.00508 * lnT4)))))
    # collisional de-excitation; the 1e-20 cm^-3 floor only guards the derivative at n_e = 0
    deexcitation = 1 + 0.12 * ((n_e + 1.0e-20) / 100) ** (0.38 - 0.12 * lnT4)
    # (1 - S): sigmoid handoff to the CIE metal tables over 2e4-3.5e4 K
    cie_taper = logistic(-10 * (T - 2.75e4) / 1.5e4)
    lower_edge = logistic((sp.log(T, 10) - sp.log(T_NEB_MIN, 10)) / 0.01)
    return 3.68e-23 * sp.exp(-3.86 / T4) / sp.sqrt(T4) * 10**log10_fneb / deexcitation * cie_taper * lower_edge


nebular_cooling = NBodyProcess(
    {"e-", "H+"},
    heat_rate_coefficient=-f_neb * Z_dust * nebular_cooling_coefficient(),
    name="Nebular forbidden-line cooling",
    bibliography=["2023ApJS..264...10K"],
)
