"""The band layer (jaco.bands): band sets, projections of absorbers and emitters onto them (closed forms against
quadrature, invariance under splitting a band, energy closure of emission), overrides, the slope fit to GIZMO's
ionizing spectrum, and the band-band coupling check."""

import warnings

import numpy as np
import pytest
import sympy as sp

from jaco.bands import (Absorber, Band, BandSet, Blackbody, Continuum, Line, PowerLaw, Projector, SIGMA_HI, SIGMA_HEI,
                        SIGMA_HEII, STARFORGE_RT, Spectral, TTable, ThermalEmission, band_structure,
                        combined_ionizing_residuals, fit_slope, hardening_mismatch, integrate, ionizing_fits, planck)
from jaco.bands.projection import mean_photon_energy, ppl_absorption, spectrum_absorption
from jaco.bands.spectra import SIGMA_SB, power_integral
from jaco.model import Model
from jaco.processes import Transfer
from jaco.symbols import n_

HI = Absorber("H photoionization", SIGMA_HI, absorber="H", E_th=13.6)
T_SMALL = tuple(np.logspace(0.0, 5.0, 51))


# --- spectra and bands ---------------------------------------------------------------------------------------------

def test_planck_integrates_to_sigma_T4_over_pi():
    T = np.array([3.0, 300.0, 3e4, 3e6])
    assert np.allclose(integrate(planck, 1e-9, 1e7, T) / (SIGMA_SB * T**4 / np.pi), 1, rtol=1e-9)


def test_power_integral_through_minus_one():
    for p in (-1.0, -1.0 + 1e-14, -2.5, 0.0, 3.0):
        quad = integrate(lambda E: E**p, 2.0, 50.0)
        assert power_integral(p, 2.0, 50.0) == pytest.approx(quad, rel=1e-10)


def test_verner_thresholds():
    # Verner et al. 1996 threshold cross-sections: 6.30, 7.42 and 1.58 Mb (hydrogenic Z^-2 scaling for He+)
    assert SIGMA_HI(13.6) == pytest.approx(6.35e-18, rel=0.01)
    assert SIGMA_HEI(24.59) == pytest.approx(7.42e-18, rel=0.01)
    assert SIGMA_HEII(54.42) == pytest.approx(SIGMA_HI(13.6) / 4, rel=0.01)
    assert SIGMA_HI(13.5) == 0 and SIGMA_HEI(24.5) == 0


def test_band_set():
    bands = BandSet([Band("EUV", 13.6, 500.0, unit="photons"), Band("FUV", 8.0, 13.6), Band("IR", 0.001, 0.4133,
                                                                                           shape="tracked")])
    assert bands.index("FUV") == 1 and bands.names == ("EUV", "FUV", "IR")
    assert bands.band_at(10.0).name == "FUV" and bands.band_at(13.6).name == "EUV" and bands.band_at(1.0) is None
    assert bands.gaps(1e-4, 1e3) == [(1e-4, 0.001), (0.4133, 8.0), (500.0, 1e3)]
    split = bands.split("FUV", [11.2], names=["FUV_lo", "FUV_hi"])
    assert split.names == ("EUV", "FUV_lo", "FUV_hi", "IR") and split["FUV_hi"].E_lo == 11.2
    assert [d.kind for d in split.declarations()] == ["radiation"] * 4
    assert split["FUV_lo"].species == "photon_FUV_lo"
    with pytest.raises(ValueError, match="overlap"):
        BandSet([Band("a", 1.0, 3.0), Band("b", 2.0, 4.0)])
    with pytest.raises(ValueError, match="duplicate"):
        BandSet([Band("a", 1.0, 2.0), Band("a", 3.0, 4.0)])
    with pytest.raises(ValueError):
        Band("x", 2.0, 1.0)
    with pytest.raises(ValueError):
        bands.split("FUV", [20.0])


# --- absorption ----------------------------------------------------------------------------------------------------

def quadrature_twin(absorber):
    """The same absorber with its power laws hidden from the closed forms"""
    def hide(f):
        return None if f is None else Spectral(f, E_min=getattr(f, "E_min", 0.0))
    return Absorber(absorber.name, hide(absorber.cross_section), absorber.absorber, absorber.E_th, absorber.heat_to,
                    hide(absorber.heat_yield), absorber.remainder_to, hide(absorber.scattering))


@pytest.mark.parametrize("slope", [-1.0, -3.7, 0.5, 0.0, -2.0])
def test_closed_forms_match_quadrature(slope):
    """Eq. 30 of He, Wibking & Krumholz (2024) and the per-photon energies, against quadrature, with a threshold inside
    the band, a constant yield, and scattering; slope 0 with sigma ~ E^-1 hits the logarithmic special case"""
    absorbers = [
        Absorber("pl", PowerLaw(2e-18, 13.6, -3.0, E_min=13.6), E_th=13.6),
        Absorber("pl inside", PowerLaw(1e-21, 1.0, -1.0, E_min=10.0), E_th=2.0, heat_yield=0.3,
                 remainder_to="dust heat", scattering=PowerLaw(5e-22, 1.0, -0.5)),
    ]
    for a in absorbers:
        for lo, hi in ((13.6, 500.0), (8.0, 13.6), (5.0, 40.0)):
            closed = ppl_absorption(a, lo, hi, slope, "exact")
            quad = ppl_absorption(quadrature_twin(a), lo, hi, slope, "exact")
            if closed is None:
                assert quad is None
                continue
            for k in closed:
                assert closed[k] == pytest.approx(quad[k], rel=1e-9, abs=1e-30), (a.name, lo, hi, k)
            assert closed["E_abs"] == pytest.approx(a.E_th + closed["heat"] + closed["remainder"], rel=1e-12)


def test_edges_representation():
    """In a ppl band, opacity "edges" is the power law through the band-edge values: identical for a power law, and
    for a curved cross-section the closed form of that power law"""
    pl = Absorber("pl", PowerLaw(2e-18, 13.6, -3.0, E_min=13.6), E_th=13.6)
    for k, v in ppl_absorption(pl, 13.6, 100.0, -1.0, "edges").items():
        assert v == pytest.approx(ppl_absorption(pl, 13.6, 100.0, -1.0, "exact")[k], rel=1e-12)
    s = np.log(SIGMA_HI(500.0) / SIGMA_HI(13.6)) / np.log(500 / 13.6)
    edge = Absorber("edge", PowerLaw(float(SIGMA_HI(13.6)), 13.6, s, E_min=13.6), E_th=13.6)
    for k, v in ppl_absorption(HI, 13.6, 500.0, -1.0, "edges").items():
        assert v == pytest.approx(ppl_absorption(edge, 13.6, 500.0, -1.0, "exact")[k], rel=1e-12)


def band_totals(proj, absorber, weights):
    """Absorptions, absorbed power, heat and momentum summed over the bands, with photons N_b and energy U_b in each"""
    rows = np.zeros(4)
    for b in proj.bands.names:
        c = proj.absorption(absorber, b)
        if c is None:
            continue
        N, U = weights[b]
        rows += [c.sigma_N * N, c.chi_E * U, c.heat * c.sigma_N * N, c.chi_F * U]
    return rows


def test_split_invariance():
    """Sub-bands whose spectra are the parent's restricted to them absorb, heat and push exactly as the parent does,
    for any cross-section under opacity "exact" and for a power law under "edges"; a curved cross-section under
    "edges" converges to it"""
    curved = Absorber("curved", Spectral(lambda E: 1e-21 * (E**-1.3 + 0.2 * np.exp(-(E - 5.0) ** 2)), E_min=3.0),
                      E_th=3.0, heat_yield=lambda E: 0.1 + 0.02 * E, remainder_to="dust heat", scattering=3e-22)
    straight = Absorber("straight", PowerLaw(1e-21, 1.0, -1.3, E_min=3.0), E_th=3.0)
    parent = BandSet([Band("P", 1.0, 20.0, slope=-1.7)])

    def totals(bands, absorber, opacity):
        w = {b.name: (power_integral(b.slope - 1, b.E_lo, b.E_hi), power_integral(b.slope, b.E_lo, b.E_hi))
             for b in bands}
        return band_totals(Projector(bands, [absorber], opacity=opacity), absorber.name, w)

    for absorber, opacity in ((curved, "exact"), (straight, "exact"), (straight, "edges")):
        ref = totals(parent, absorber, opacity)
        for edges in ([2.5], [2.5, 7.0, 11.2], list(np.geomspace(1.0, 20.0, 17)[1:-1])):
            assert np.allclose(totals(parent.split("P", edges), absorber, opacity), ref, rtol=1e-10, atol=0)
    ref = totals(parent, curved, "exact")
    errors = [np.max(np.abs(totals(parent.split("P", list(np.geomspace(1.0, 20.0, n + 1)[1:-1])), curved, "edges")
                            / ref - 1)) for n in (2, 8, 32)]
    assert errors[0] > errors[1] > 8 * errors[2]


def test_split_invariance_tracked():
    """The same for a tracked band, whose spectrum is a blackbody at T_rad: at every tabulated T_rad"""
    k = Absorber("dust", Spectral(lambda E: 50.0 * E**-1.5 + 3.0))
    parent = BandSet([Band("IR", 0.001, 0.4133, shape="tracked")])
    T = np.asarray(T_SMALL)

    def totals(bands):
        proj = Projector(bands, [k], T_grid=T_SMALL)
        U = {b.name: integrate(planck, b.E_lo, b.E_hi, T) for b in bands}
        N = {b.name: integrate(lambda E, t: planck(E, t) / E, b.E_lo, b.E_hi, T) for b in bands}
        out = 0
        for b in bands.names:
            c = proj.absorption("dust", b)
            out = out + np.array([c.sigma_N.values * N[b], c.chi_E.values * U[b]])
        return out

    ref = totals(parent)
    split = totals(parent.split("IR", [0.01, 0.05]))
    ok = ref[1] > 1e-200
    assert np.allclose(split[:, ok], ref[:, ok], rtol=1e-9, atol=0)


def test_tracked_band_is_planck_mean():
    """A tracked band's energy-mean opacity at T_rad is the band's Planck mean at T = T_rad"""
    kappa = Spectral(lambda E: 50.0 * E**-1.5 + 3.0)
    bands = BandSet([Band("IR", 0.001, 0.4133, shape="tracked")])
    proj = Projector(bands, [Absorber("dust", kappa)], [ThermalEmission("dust emission", kappa)], T_grid=T_SMALL)
    planck_mean = proj.emission("dust emission").planck_mean["IR"].values
    assert np.allclose(proj.absorption("dust", "IR").chi_E.values, planck_mean, rtol=1e-12, atol=0)
    hnu = proj.mean_photon_energy("IR")
    assert isinstance(hnu, TTable) and hnu(1e6) < 0.4133 and hnu(100.0) == pytest.approx(2.70 * 8.617e-3, rel=0.01)


def test_mean_photon_energy_and_no_overlap():
    bands = BandSet([Band("FUV", 8.0, 13.6), Band("EUV", 13.6, 500.0, unit="photons", slope=-3.0)])
    proj = Projector(bands, [HI])
    assert proj.mean_photon_energy("FUV") == pytest.approx(np.log(13.6 / 8.0) / (1 / 8.0 - 1 / 13.6), rel=1e-12)
    assert proj.absorption("H photoionization", "FUV") is None
    assert set(proj.absorptions()) == {("H photoionization", "EUV")}


# --- emission ------------------------------------------------------------------------------------------------------

GIZMO_LIKE = BandSet([Band("IR", 0.001, 0.4133, shape="tracked"), Band("ONIR", 0.4133, 3.444), Band("NUV", 3.444, 8.0),
                      Band("LW", 11.2, 13.6), Band("EUV", 13.6, 500.0, unit="photons")])  # gap 8-11.2 eV


@pytest.mark.parametrize("kirchhoff", ["band", "planck"])
def test_thermal_emission_closure(kirchhoff):
    """Band fractions plus the escaping remainder sum to one, and the band powers 4 pi chi_B,b B_b plus the escape
    make the total 4 pi kappa_planck sigma T^4 / pi. With "planck" the total is the exact Planck mean's, integrated
    independently of the bands; with "band" a ppl band emits at its energy-mean opacity chi_E,b and a tracked band at
    its Planck mean, which is the same thing for a blackbody spectrum"""
    kappa = Spectral(lambda E: 3.0 / E + 0.5)
    proj = Projector(GIZMO_LIKE, [Absorber("dust", kappa)], [ThermalEmission("dust emission", kappa)],
                     kirchhoff=kirchhoff)
    e = proj.emission("dust emission")
    T = e.escape.T
    total = sum(f.values for f in e.fraction.values()) + e.escape.values
    assert np.allclose(total, 1, rtol=1e-10)
    power = sum(e.planck_mean[b].values * e.band_planck[b].values for b in GIZMO_LIKE.names)
    expected = e.kappa_planck.values * SIGMA_SB * T**4 / np.pi * (1 - e.escape.values)
    hot = T > 30  # where every band's Planck integral is resolved in double precision
    assert np.allclose(power[hot], expected[hot], rtol=1e-8, atol=0)
    i = np.searchsorted(T, 5e4)
    assert 0.01 < e.escape.values[i] < 0.2 and e.fraction["IR"](100.0) > 0.99 and e.escape(1.0) > 0.99
    exact = Projector(GIZMO_LIKE, emitters=[ThermalEmission("dust emission", kappa)], kirchhoff="planck")
    exact = exact.emission("dust emission")
    assert np.allclose(e.planck_mean["IR"].values, exact.planck_mean["IR"].values, rtol=1e-12, atol=0)
    if kirchhoff == "band":
        for b in ("ONIR", "NUV", "LW", "EUV"):
            assert np.all(e.planck_mean[b].values == proj.absorption("dust", b).chi_E)
        assert not np.allclose(e.kappa_planck.values, exact.kappa_planck.values, rtol=1e-3)


def test_planck_mean_survives_underflow():
    """For a constant opacity every band's Planck mean is that constant, also where B_E underflows (EUV at 1 K)"""
    e = Projector(GIZMO_LIKE, emitters=[ThermalEmission("grey", 7.0)]).emission("grey")
    for b in GIZMO_LIKE.names:
        assert np.allclose(e.planck_mean[b].values, 7.0, rtol=1e-12)
    assert np.allclose(e.kappa_planck.values, 7.0, rtol=1e-9)


def test_continuum_and_line_closure():
    ff = Continuum("free-free", Spectral(lambda E, T: np.exp(-E / (8.617e-5 * T)) / np.sqrt(T),
                                         temperature_dependent=True), E_range=(1e-6, 1e5))
    flat = Continuum("flat", lambda E: np.ones_like(E), E_range=(1.0, 20.0))
    lines = [Line("Lyman alpha", 10.2), Line("H alpha", 1.89)]
    proj = Projector(GIZMO_LIKE, emitters=[ff, flat, *lines])
    e = proj.emission("free-free")
    assert np.allclose(sum(f.values for f in e.fraction.values()) + e.escape.values, 1, rtol=1e-10)
    e = proj.emission("flat")
    assert e.escape == pytest.approx((11.2 - 8.0) / 19.0, rel=1e-10)
    assert e.fraction["NUV"] == pytest.approx((8.0 - 3.444) / 19.0, rel=1e-10)
    assert e.photons_per_eV["NUV"] == pytest.approx(np.log(8.0 / 3.444) / 19.0, rel=1e-10)
    assert proj.emission("Lyman alpha").escape == 1.0  # in the 8-11.2 eV gap
    halpha = proj.emission("H alpha")
    assert halpha.fraction["ONIR"] == 1.0 and halpha.escape == 0.0 and halpha.photons_per_eV["ONIR"] == 1 / 1.89


# --- overrides and symbolic forms ----------------------------------------------------------------------------------

def test_overrides():
    sigma, eps = sp.symbols("sigma_HI eps_HI")
    bands = BandSet([Band("EUV", 13.6, 500.0, unit="photons"), Band("ONIR", 0.4133, 3.444)])
    proj = Projector(bands, [HI], [Line("Halpha", 1.89)],
                     overrides={("H photoionization", "EUV"): {"sigma_N": sigma, "heat": eps},
                                ("H photoionization", "ONIR"): {"chi_E": 1e-22},
                                ("Halpha", "ONIR"): {"fraction": 0.5}})
    c = proj.absorption("H photoionization", "EUV")
    assert c.sigma_N == sigma and c.heat == eps and c.overridden == {"sigma_N", "heat"}
    assert c.chi_E == pytest.approx(Projector(bands, [HI]).absorption("H photoionization", "EUV").chi_E)
    c = proj.absorption("H photoionization", "ONIR")  # no overlap: the override alone couples them
    assert c.chi_E == 1e-22 and c.sigma_N == 0
    e = proj.emission("Halpha")
    assert e.fraction["ONIR"] == 0.5 and e.escape == 0.5
    with pytest.raises(KeyError, match="unknown coefficients"):
        Projector(bands, [HI], overrides={("H photoionization", "EUV"): {"sigma": 1.0}})
    with pytest.raises(KeyError):
        Projector(bands, [HI], overrides={("nope", "EUV"): {"chi_E": 1.0}})
    with pytest.raises(KeyError):
        Projector(bands, [HI], overrides={("H photoionization", "FUV"): {"chi_E": 1.0}})


def test_symbolic_forms():
    T_rad = sp.Symbol("T_rad")
    bands = BandSet([Band("EUV", 13.6, 500.0, unit="photons"), Band("IR", 0.001, 0.4133, shape="tracked")])
    kappa = Absorber("dust", Spectral(lambda E: 50.0 * E**-1.5 + 3.0))
    proj = Projector(bands, [HI, kappa], T_grid=T_SMALL)
    c = proj.absorption("H photoionization", "EUV")
    assert c.expr("sigma_N") == sp.Float(c.sigma_N)
    c = proj.absorption("dust", "IR")
    e = c.expr("chi_E", T_rad)
    assert float(e.subs(T_rad, 300.0).doit()) == pytest.approx(c.chi_E(300.0), rel=1e-12)


# --- GIZMO's ionizing band -----------------------------------------------------------------------------------------

def gizmo_rt_get_sigma(T_eff=4e4, n=10000):
    """rt_get_sigma (radiation/rt_chem.cc) for its single H-ionizing band: rectangle sums over 13.6-500 eV of a
    blackbody, with the hydrogenic cross-section and GIZMO's constants. Returns (sigma_HI, G_HI, nu_eff)"""
    eV, k_B, hc = 1.60217733e-12, 1.38066e-16, 2.9979e10 * 6.6262e-27
    e = 13.6 + np.arange(n) * (500.0 - 13.6) / (n - 1)
    I_nu = 2.0 * (e * eV) ** 3 / hc**2 / np.expm1(e * eV / (k_B * T_eff))
    with np.errstate(divide="ignore", invalid="ignore"):
        f = np.sqrt(e / 13.6 - 1.0)
        sig = 6.3e-18 * (13.6 / e) ** 4 * np.exp(4 - 4 * np.arctan(f) / f) / (1.0 - np.exp(-2 * np.pi / f))
    sig = np.where(e <= 13.6, 6.3e-18, sig)
    n_photons = I_nu / e
    return (np.sum(sig * n_photons) / np.sum(n_photons), np.sum(sig * (e - 13.6) * n_photons) / np.sum(sig * n_photons),
            np.sum(I_nu) / np.sum(n_photons))


def test_gizmo_ionizing_band_reference():
    """GIZMO prints sigma_HI=3.368e-18 G_HI=2.99 nu_eff=18.6 for T_eff = 4e4 K; quadrature with Verner et al. (1996)
    gives the same to its rectangle rule's accuracy"""
    s, G, nu = gizmo_rt_get_sigma()
    assert (f"{s:.3e}", f"{G:.3g}", f"{nu:.3g}") == ("3.368e-18", "2.99", "18.6")
    ref = spectrum_absorption(HI, 13.6, 500.0, Blackbody(4e4))
    assert ref["sigma_N"] == pytest.approx(s, rel=2e-3)
    assert ref["E_abs"] - 13.6 == pytest.approx(G, rel=1e-2)
    assert mean_photon_energy(13.6, 500.0, Blackbody(4e4)) == pytest.approx(nu, rel=2e-3)


EUV_TARGETS = {"G_HI": ("H photoionization", "excess"), "nu_eff": "mean_photon_energy",
               "sigma_HI": ("H photoionization", "sigma_N")}


@pytest.mark.parametrize("opacity", ["edges", "exact"])
def test_euv_slope_fit(opacity):
    """One slope cannot reproduce GIZMO's 4e4 K ionizing band: matching nu_eff leaves G_HI ~20% low, matching G_HI
    leaves nu_eff 15-22% high; the joint fit splits the difference"""
    band = Band("EUV", 13.6, 500.0, unit="photons")
    bb = Blackbody(4e4)
    report = list(EUV_TARGETS.values())
    fits = {name: fit_slope(band, bb, [t], [HI], opacity, report=report) for name, t in EUV_TARGETS.items()
            if name != "sigma_HI"}
    joint = fit_slope(band, bb, [EUV_TARGETS["G_HI"], EUV_TARGETS["nu_eff"]], [HI], opacity, report=report)
    for name, f in fits.items():
        assert abs(f.rows[EUV_TARGETS[name]][2]) < 1e-9
    assert fits["nu_eff"].slope < joint.slope < fits["G_HI"].slope
    res = {name: joint.rows[t][2] for name, t in EUV_TARGETS.items()}
    assert res["G_HI"] < -0.05 and res["nu_eff"] > 0.05 and abs(res["sigma_HI"]) < 0.1
    print(f"\nopacity {opacity}:\n" + "\n".join(f.table() for f in [*fits.values(), joint]))


def test_starforge_rt_spec():
    """STARFORGE_RT: the ionizing band split at 24.59 eV, each sub-band's slope fitted to GIZMO's 4e4 K blackbody
    (G_HI and nu_eff jointly). Per sub-band G_HI and nu_eff within 1%, sigma_HI within 1%; over 13.6-500 eV, with
    each sub-band holding the blackbody's photons in it, GIZMO's sigma_HI, G_HI and nu_eff within 0.6%"""
    assert STARFORGE_RT.names == ("EUV_H", "EUV_He", "FUV", "NUV", "ONIR", "IR")
    assert [b.unit for b in STARFORGE_RT] == ["photons"] * 2 + ["energy"] * 4
    assert STARFORGE_RT["IR"].shape == "tracked" and STARFORGE_RT["EUV_H"].E_hi == STARFORGE_RT["EUV_He"].E_lo == 24.59
    for f in ionizing_fits():
        assert all(abs(res) < 0.01 for _, _, res in f.rows.values()), f.table()
    ionizing = [STARFORGE_RT["EUV_H"], STARFORGE_RT["EUV_He"]]
    assert all(abs(r) < 0.006 for r in combined_ionizing_residuals(ionizing))
    print("\n" + "\n".join(f.table() for f in ionizing_fits()))
    print("over 13.6-500 eV: sigma_HI, G_HI, nu_eff residuals",
          ", ".join(f"{r:+.2%}" for r in combined_ionizing_residuals(ionizing)))


def test_ionizing_hardening_mismatch():
    """A photons band loses its mean photon energy hnu per absorption while the products receive E_abs < hnu: the
    cross-section takes the soft photons, so the true spectrum hardens and a fixed slope cannot. The mismatch is mostly
    the blackbody's own (absorbed photons are softer than the mean), the fixed slope adds <= 1.3 points; the split
    halves it against one band"""
    rows = [hardening_mismatch(STARFORGE_RT[b]) for b in ("EUV_H", "EUV_He")]
    one = fit_slope(Band("EUV", 13.6, 500.0, unit="photons"), Blackbody(4e4), ["mean_photon_energy"], [HI]).band
    single = hardening_mismatch(one)
    print("\nband      hnu [eV]  E_abs [eV]  (hnu - E_abs)/hnu   blackbody: hnu  E_abs  mismatch")
    for r in rows + [single]:
        print(f"{r.band:8s} {r.hnu:9.3f} {r.E_abs:11.3f} {r.mismatch:14.2%}   {r.hnu_reference:15.3f} "
              f"{r.E_abs_reference:6.3f} {r.mismatch_reference:8.2%}")
    for r in rows:
        assert 0.03 < r.mismatch_reference < r.mismatch < r.mismatch_reference + 0.015
    assert single.mismatch > 2 * max(r.mismatch for r in rows)


# --- band-band structure -------------------------------------------------------------------------------------------

def test_band_structure():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        from jaco.models.starforge import radiation as rt
    euv_donation = rt.photoionization(donation=rt.ONIR)
    dust = [rt.dust_band_absorption(b) for b in (rt.FUV, rt.NUV, rt.ONIR)]
    s = band_structure(dust, rt.BANDS)
    assert s.diagonal and s.metadata() == {"band_jacobian_diagonal": True, "band_couplings": []}
    s = band_structure(dust + [euv_donation], rt.BANDS)
    assert not s.diagonal and s.pairs == [("photon_EUV", "photon_ONIR")]
    # through a derived quantity: a row of the IR band reading G_0, an expression of the FUV band
    G0 = sp.Symbol("G_0")
    proc = Transfer(1e-30 * G0 * n_("H"), {"photon_IR": 1}, name="G_0-driven emission")
    assert band_structure([proc], rt.BANDS).diagonal
    model = Model([proc], derived={"G_0": 1e3 * n_("photon_FUV")})
    assert band_structure(model, rt.BANDS).pairs == [("photon_FUV", "photon_IR")]


@pytest.mark.slow
def test_legacy_rt_model_band_structure():
    """GIZMO's coupling is not band-diagonal: the ionizing band donates to the optical one, the dust-absorbed bands'
    energy is copied into the IR band, and cooling routed into the bands depends on other bands"""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        from jaco.models.starforge_legacy_RT import make_model
        from jaco.models.starforge.radiation import BANDS
        s = band_structure(make_model(), BANDS)
    assert not s.diagonal
    assert {("photon_EUV", "photon_ONIR"), ("photon_FUV", "photon_IR"), ("photon_NUV", "photon_IR"),
            ("photon_ONIR", "photon_IR")} <= set(s.pairs)
