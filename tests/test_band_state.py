"""State-dependent opacities in the band layer: a separable scale A(state), a cross-section tabulated in one state
variable (tables in s, or in (T_rad, s) on a tracked band), jumps in the state reported rather than smoothed, and dust
opacities from Draine's tables with their band means against GIZMO's constants."""

import warnings

import numpy as np
import pytest
import sympy as sp

from jaco.bands import (Absorber, Band, BandSet, DiscontinuityWarning, GIZMO_STARFORGE, Projector,
                        STARFORGE_RT, Spectral, Table2D, TTable, band_means, draine_mw31, read_draine)
from jaco.bands.dust import DATA_DIR, DRAINE_MW31_FILE, HC_EV_MICRON
from jaco.bands.projection import ppl_absorption

Td, Zd, x_e = sp.symbols("Td Z_d x_e")
KAPPA0 = Spectral(lambda E: 50.0 * E**-1.5 + 3.0)
S_GRID = tuple(np.logspace(1.0, 3.5, 51))
T_GRID = tuple(np.logspace(0.0, 4.0, 41))
BANDS = BandSet([Band("NUV", 3.444, 8.0), Band("IR", 0.001, 0.4133, shape="tracked", temperature="T_rad")])


def test_separable_scale():
    """sigma_N, chi_E and chi_S are the projection times the scale (scattering_scale for chi_S); chi_F is their sum"""
    survival = sp.exp(-Td / 1500)
    plain = Projector(BANDS, [Absorber("dust", KAPPA0, scattering=2.0)]).absorption("dust", "NUV")
    scaled = Projector(BANDS, [Absorber("dust", KAPPA0, scattering=2.0, scale=Zd * survival,
                                        scattering_scale=x_e)]).absorption("dust", "NUV")
    assert plain.chi_F == pytest.approx(plain.chi_E + plain.chi_S, rel=1e-12) and plain.chi_S == pytest.approx(2.0)
    assert scaled.chi_E == plain.chi_E  # the numbers are the plain projection
    assert scaled.expr("chi_E") == sp.Float(plain.chi_E) * Zd * survival
    assert scaled.expr("sigma_N") == sp.Float(plain.sigma_N) * Zd * survival
    assert scaled.expr("chi_S") == sp.Float(plain.chi_S) * x_e
    assert sp.simplify(scaled.expr("chi_F") - (sp.Float(plain.chi_E) * Zd * survival + 2.0 * x_e)) == 0
    assert scaled.expr("E_abs") == sp.Float(plain.E_abs)  # energies per photon are not scaled
    same = Projector(BANDS, [Absorber("dust", KAPPA0, scattering=2.0, scale=Zd)]).absorption("dust", "NUV")
    assert same.expr("chi_F") == sp.Float(plain.chi_F) * Zd
    assert plain.expr("chi_E") == sp.Float(plain.chi_E)  # no scale, no factor


def growth(s):
    return 1.0 + (s / 300.0) ** 2


def test_state_tables_ppl_and_tracked():
    """kappa(E; s) = g(s) kappa0(E): the tables equal the fixed-state projections at the nodes and g(s) times the
    plain projection between them to the interpolation's accuracy; on a tracked band, a Table2D in (T_rad, s) whose
    Kirchhoff value at (T_d, T_d) evaluates from its expression"""
    a = Absorber("dust", lambda E, s: growth(s) * KAPPA0(E), state="Td", state_grid=S_GRID)
    proj = Projector(BANDS, [a], T_grid=T_GRID)
    plain = Projector(BANDS, [Absorber("dust", KAPPA0)], T_grid=T_GRID)
    nuv, ir = proj.absorption("dust", "NUV"), proj.absorption("dust", "IR")
    assert isinstance(nuv.chi_E, TTable) and isinstance(ir.chi_E, Table2D) and not nuv.discontinuities
    ref = plain.absorption("dust", "NUV")
    assert np.allclose(nuv.chi_E.values, growth(np.array(S_GRID)) * ref.chi_E, rtol=1e-12)
    assert np.allclose(nuv.E_abs.values, ref.E_abs, rtol=1e-12)  # a separable state leaves energies alone
    s = 137.0
    assert nuv.chi_E(s) == pytest.approx(growth(s) * ref.chi_E, rel=2e-3)
    assert float(nuv.expr("chi_E").subs(Td, s).doit()) == pytest.approx(nuv.chi_E(s), rel=1e-10)
    ref_ir = plain.absorption("dust", "IR").chi_E
    for T_rad in (3.0, 30.0, 300.0):
        assert ir.chi_E(T_rad, s) == pytest.approx(growth(s) * ref_ir(T_rad), rel=3e-3)
    kirchhoff = ir.expr("chi_E", Td, Td)
    assert kirchhoff.free_symbols == {Td}
    assert float(kirchhoff.subs(Td, s)) == pytest.approx(float(ir.chi_E(s, s)), rel=1e-10)
    assert ir.expr("chi_E").free_symbols == {sp.Symbol("T_rad"), Td}


def zoned(E, s, width=0.0):
    """A composition switch at 160 K, a step (width 0) or a smoothstep of the given width"""
    if width == 0:
        f = np.where(s < 160.0, 1.0, 0.4)
    else:
        x = np.clip((s - 160.0 + width / 2) / width, 0, 1)
        f = 1 - 0.6 * x * x * (3 - 2 * x)
    return f * KAPPA0(E)


def test_jump_in_state_is_reported():
    a = Absorber("zoned", zoned, state="Td", state_grid=S_GRID)
    with pytest.warns(DiscontinuityWarning, match="x0.4 at Td = 160"):
        c = Projector(BANDS, [a], T_grid=T_GRID).absorption("zoned", "NUV")
    (field, s, factor), = c.discontinuities
    assert field == "chi_E" and s == pytest.approx(160.0, rel=1e-6) and factor == pytest.approx(0.4, rel=1e-6)
    with pytest.raises(ValueError, match="jumps in Td"):
        Projector(BANDS, [a], T_grid=T_GRID, discontinuities="raise").absorption("zoned", "IR")
    smooth = Absorber("smooth", lambda E, s: zoned(E, s, width=10.0), state="Td", state_grid=S_GRID)
    with warnings.catch_warnings():
        warnings.simplefilter("error", DiscontinuityWarning)
        p = Projector(BANDS, [smooth], T_grid=T_GRID)
        assert p.absorption("smooth", "NUV").discontinuities == ()
        assert p.absorption("smooth", "IR").discontinuities == ()


def test_state_declaration_checks():
    with pytest.raises(ValueError, match="log-uniform"):
        Absorber("x", lambda E, s: E, state="s", state_grid=(1.0, 2.0, 5.0))
    with pytest.raises(ValueError, match="positive state_grid"):
        Absorber("x", lambda E, s: E, state="s")
    with pytest.raises(ValueError, match="without a state"):
        Absorber("x", 1.0, state_grid=S_GRID)
    with pytest.raises(ValueError, match="at_state"):
        ppl_absorption(Absorber("x", lambda E, s: E, state="s", state_grid=S_GRID), 1.0, 2.0, -1.0)


# --- dust tables -----------------------------------------------------------------------------------------------------

DRAINE_LIKE = """Extinction, Scattering, and Absorption Properties (format test)

Tabulated quantities:
lambda  = wavelength in vacuo (micron)

2.000E-26 = M_dust per H nucleon (gram/H) for this dust model
1.000E+02 = M_gas/M_dust for this dust model (assuming He/H=0.096)

  lambda    albedo   <cos>  C_ext/H    K_abs   <cos^2>
 (micron)                   (cm^2/H)  (cm^2/g)          comment
1.00000E+01 0.0000  0.0000 2.000E-24 1.000E+02 0.40000
1.00000E+00 0.5000  0.5000 4.000E-22 1.000E+04 0.40000 some comment
1.00000E-01 0.2500  0.6000 8.000E-21 3.000E+05 0.50000
"""


def test_read_draine_format(tmp_path):
    path = tmp_path / "kext_test.all"
    path.write_text(DRAINE_LIKE)
    d = read_draine(path)
    assert np.allclose(d.E, HC_EV_MICRON / np.array([10.0, 1.0, 0.1]))
    assert np.allclose(d.kappa_abs, [1e2, 1e4, 3e5]) and np.allclose(d.kappa_sca, [0.0, 1e4, 1e5])
    assert d.dust_to_gas == pytest.approx(0.01)
    E = np.sqrt(d.E[1] * d.E[2])
    assert d.absorption(E) == pytest.approx(np.sqrt(1e4 * 3e5), rel=1e-12)  # log-log between the rows
    assert d.absorption(1e-6) == pytest.approx(1e2) and d.absorption(1e6) == pytest.approx(3e5)  # clamped
    gas = Projector(BandSet([Band("V", 1.5, 3.0)]), [d.absorber()]).absorption("dust", "V")
    per_dust = Projector(BandSet([Band("V", 1.5, 3.0)]), [d.absorber(per="dust")]).absorption("dust", "V")
    assert gas.chi_E == pytest.approx(0.01 * per_dust.chi_E, rel=1e-12)
    assert gas.chi_S == pytest.approx(0.01 * per_dust.chi_S, rel=1e-12)


GIZMO_RT_KAPPA = {"FUV": 720.0, "NUV": 480.0, "ONIR": 180.0}  # rt_kappa, cm^2/g at solar Z, extinction; half absorbed


@pytest.mark.skipif(not (DATA_DIR / DRAINE_MW31_FILE).is_file(), reason=f"{DRAINE_MW31_FILE} not in jaco/bands/data")
def test_draine_mw31_against_gizmo():
    """Draine's MW R_V=3.1 table: its own consistency (K_abs = (1 - albedo) C_ext/H / M_dust, A_V/N_H of the
    D03 normalization) and its band means per gram of gas against GIZMO's rt_kappa"""
    d = draine_mw31()
    assert d.dust_to_gas == pytest.approx(1 / 165.3, rel=1e-3)
    E_V = HC_EV_MICRON / 0.55
    kappa_ext_V = d.absorption(E_V) + d.scattering(E_V)
    A_V_per_N_H = 1.086 * kappa_ext_V * 1.398e-26  # mag cm^2
    assert A_V_per_N_H == pytest.approx(5.3e-22, rel=0.1)
    for label, bands in (("GIZMO_STARFORGE", GIZMO_STARFORGE), ("STARFORGE_RT", STARFORGE_RT)):
        means = band_means(d, bands)
        print(f"\n{label}: band means per gram of gas [cm^2/g] (chi_E absorption, chi_F extinction)")
        for b, (chi_E, chi_S, chi_F) in means.items():
            if bands[b].shape == "tracked":
                print(f"  {b}: chi_E at T_rad = 10/30/100 K: " + ", ".join(f"{chi_E(T):.3g}" for T in (10, 30, 100)))
                continue
            g = GIZMO_RT_KAPPA.get(b)
            ratios = f"  chi_F/{g:g} = {chi_F / g:.2f}, chi_E/{g / 2:g} = {chi_E / (g / 2):.2f}" if g else ""
            print(f"  {b}: chi_E {chi_E:.4g}, chi_F {chi_F:.4g}, albedo {chi_S / chi_F:.2f}{ratios}")
            if g:
                assert 0.1 < chi_F / g < 10
