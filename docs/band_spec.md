# Radiation bands as a spec (`jaco.bands`)

Status: core layer, phase P2. Nothing in `jaco.models` uses it yet; the generated code is unchanged.
Decided 2026-10-06 (see "Decisions"): STARFORGE_RT splits the ionizing band at 24.59 eV, in-band opacities are exact
quadrature by default, thermal emission into a ppl band is Kirchhoff per band by default.

Processes declare which photon energies they interact with: a cross-section sigma(E) or opacity kappa(E) for an
absorber, a spectrum j(E) for an emitter. A `BandSet` turns the declarations into per-band coefficients at
code-generation time. Refining a band (splitting FUV at 11.2 eV, EUV at the He thresholds) changes only the band set.

## API

```python
from jaco.bands import Band, BandSet, Absorber, ThermalEmission, Line, Projector, SIGMA_HI, STARFORGE_RT

bands = BandSet([
    Band("EUV", 13.6, 500.0, unit="photons", slope=-3.70),  # fixed-slope spectrum u_E ~ E^slope
    Band("FUV", 8.0, 13.6),                                    # default: energy, slope -1 (E u_E = const)
    Band("IR", 0.001, 0.4133, shape="tracked"),                # dilute blackbody at the per-cell T_rad
])
bands = bands.split("FUV", [11.2], names=["FUV_lo", "LW"])     # refinement: only the spec changes
bands = STARFORGE_RT                                           # EUV_H, EUV_He, FUV, NUV, ONIR, IR

proj = Projector(bands,
                 absorbers=[Absorber("H photoionization", SIGMA_HI, absorber="H", E_th=13.6),
                            Absorber("dust", kappa_dust, heat_to="dust heat", scattering=kappa_scat)],
                 emitters=[ThermalEmission("dust emission", kappa_dust, source="dust heat"), Line("Halpha", 1.89)],
                 overrides={("H photoionization", "EUV_H"): {"sigma_N": sp.Symbol("sigma_HI")}})
c = proj.absorption("H photoionization", "EUV_H") # AbsorptionCoefficients, or None if no overlap
e = proj.emission("dust emission")                 # EmissionFractions
proj.mean_photon_energy("EUV_He")                  # <h nu>_b
```

| Object | What it is |
|---|---|
| `Band(name, E_lo, E_hi, unit, shape, slope, temperature, doc)` | `unit` "photons" or "energy" (what the species `photon_<name>` counts); `shape` "ppl" (fixed slope) or "tracked" (blackbody at `T_rad`, the symbol named `temperature`); `doc` describes the species. |
| `BandSet` | Ordered, non-overlapping, gaps allowed; `index`, `by_species`, `band_at(E)`, `gaps`, `split`, `replace`, `declarations()` (radiation `Species`). |
| `Absorber(name, cross_section, absorber, E_th, heat_to, heat_yield, remainder_to, scattering)` | Per absorbed photon: `E_th` to chemistry, `Y(E)(E - E_th)` to `heat_to`, the rest to `remainder_to`. Scattering enters the flux mean only. `cross_section=None`: every coefficient is an override. |
| `Line(name, E0, source)`, `Continuum(name, j, source, E_range)`, `ThermalEmission(name, kappa, source)`, `RoutedEmission(name, source)` | Emitters; `j(E)` or `j(E, T)`; thermal is `4 pi kappa(E) B_E(T)`; routed: band fractions given as overrides. `source` is the row losing the energy. |
| `PowerLaw`, `Spectral`, `Blackbody`, `planck`, `SIGMA_HI/HEI/HEII` | Spectral functions of E [eV]; Verner et al. (1996) fits. A `PowerLaw` gets closed forms. |
| `Projector(bands, absorbers, emitters, overrides, band_overrides, opacity="exact", kirchhoff="band", T_grid)` | Projections, cached. Ppl bands give floats, tracked bands and emitter temperatures give `TTable`s (log-log in T). `.expr(name, T)`/`symbolic` give sympy (`Float`, a `piecewise_linear` in T, or an override expression with the band's T_rad replaced by T). Band level: `mean_photon_energy`, `compton_temperature` (T_C = ⟨E⟩_u / 4k), each overridable. |
| `fit_slope(band, reference, targets, absorbers, opacity, weights, report)` | Slope of a ppl band matching chosen band means of a reference spectrum. |
| `GIZMO_STARFORGE` (`jaco.bands.specs`) | GIZMO's five STARFORGE bands as they are; the legacy model's band set. |
| `STARFORGE_RT`, `starforge_rt(T_eff, opacity)` (`jaco.bands.specs`) | GIZMO's STARFORGE bands, ionizing band split at 24.59 eV, ionizing slopes fitted to a 4e4 K blackbody; `ionizing_fits`, `hardening_mismatch`, `combined_ionizing_residuals` report the fit. |
| `band_structure(model or processes, bands)` | Which processes make a band's row depend on another band; `.diagonal`, `.metadata()`. |

## Projection

Band b = [a, c], energy spectrum u_E inside it (ppl: E^alpha; tracked: B_E(T_rad)), photon spectrum u_E / E. For an
absorber sigma(E) (zero below max(E_th, its own support)):

| Coefficient | Definition | Use |
|---|---|---|
| `sigma_N` | ∫σ u/E / ∫u/E | absorptions per absorber: c σ_N n_γ |
| `chi_E` | ∫σ u / ∫u | absorbed power per absorber: c χ_E u_b |
| `chi_F` | ∫(σ + σ_s) u / ∫u | momentum (fixed slope: the flux has u's shape) |
| `E_abs` | ∫σ u / ∫σ u/E | mean energy of an absorbed photon; `excess` = E_abs − E_th |
| `heat`, `remainder` | ∫σ Y (E−E_th) u/E / ∫σ u/E, and excess − heat | E_th + heat + remainder = E_abs exactly |

Closed forms when σ and u are power laws: with σ = σ_0 (E/E_0)^s on [lo, c], lo = max(a, E_min), and
I_p(x, y) = ∫_x^y E^p dE (log form at p = −1),
χ_E = σ_0 E_0^−s I_{α+s}(lo, c) / I_α(a, c), σ_N = σ_0 E_0^−s I_{α+s−1}(lo, c) / I_{α−1}(a, c),
E_abs = I_{α+s}(lo, c) / I_{α+s−1}(lo, c). With lo = a this is He, Wibking & Krumholz (2024) Eq. 30. Otherwise
Gauss-Legendre quadrature in ln E (8 nodes per 0.1 in ln E), which matches the closed forms to 1e-9 or better.

`opacity="exact"` (default) projects the declared function. `"edges"` (the paper's PPL) replaces σ and σ_s inside a
ppl band by the power law through their values at the band's edges, then uses the closed forms. Tracked bands always
use the declared function.

Emission: `fraction[b]` = ∫_b j / ∫ j, `photons_per_eV[b]` = ∫_b j/E / ∫ j, and `escape` (energy in no band) is
integrated over the gaps independently, so for a line or continuum Σ fraction + escape = 1 is a check, not an identity.
A thermal emitter's band b emits 4π χ_B,b B_b(T), with `band_planck` B_b(T) and `planck_mean` χ_B,b:

- `kirchhoff="band"` (default, the paper's choice): χ_B,b = χ_E,b, the band's energy-mean opacity of κ under its
  fixed-slope spectrum (same `opacity` model), so emission and absorption in the band balance exactly at
  u_b = 4π B_b / c. The total is the sum over the bands plus the gaps; it depends on the band set.
- `kirchhoff="planck"`: the exact Planck mean ∫_b κB / ∫_b B, the total ∫κB independently of the bands.

A tracked band's spectrum is a blackbody, where the two coincide (its exact Planck mean either way, computed with the
band's lower-edge Wien factor taken out so that it stays finite where B underflows). `kappa_planck` is the emitted
power over 4π∫B, so the emitter's total is 4π κ_planck σT⁴/π under either choice.

Overrides: `{(process, band): {coefficient: value}}`, numbers or sympy expressions, replace projected values; an
override may couple an absorber to a band it does not overlap; an overridden emission fraction sets escape to
1 − Σ fractions. Results record `overridden`.

## Results (tests/test_bands.py, tests/test_band_convergence.py)

GIZMO's ionizing band (rt_get_sigma, T_eff = 4e4 K, 13.6-500 eV). Replicating its rectangle rule and hydrogenic
cross-section reproduces its printout, σ_HI = 3.368e-18 cm², G_HI = 2.99 eV, ν_eff = 18.6 eV. Quadrature with
Verner (1996): 3.370e-18, 3.012, 18.635.

STARFORGE_RT (decided): EUV_H [13.6, 24.59] and EUV_He [24.59, 500] eV, each slope fitted to the blackbody restricted
to the sub-band by least squares on ln G_HI and ln ν_eff (H photoionization only, as in GIZMO without
RT_CHEM_PHOTOION_HE). Residuals, `opacity="exact"`:

| Sub-band | slope | G_HI | ν_eff | σ_HI (not fitted) |
|---|---|---|---|---|
| EUV_H | −2.114 | −0.09% (2.699 vs 2.702 eV) | +0.43% (17.46 vs 17.38 eV) | −0.33% |
| EUV_He | −6.315 | −0.95% (13.99 vs 14.13 eV) | +0.84% (29.22 vs 28.97 eV) | +0.62% |
| 13.6-500 eV, each sub-band holding the blackbody's photons | | −0.11% | +0.50% | −0.30% |

With `opacity="edges"` the slopes are −2.103 and −6.235 and the combined residuals −0.38% (G_HI), +0.57% (ν_eff),
−1.04% (σ_HI).

A single band cannot do this: ν_eff depends on the slope alone, and G_HI needs a slope a full unit shallower.

| Single 13.6-500 eV band fit | opacity | slope | G_HI | ν_eff | σ_HI |
|---|---|---|---|---|---|
| ν_eff | exact | −3.70 | −17.4% | 0 | +8.5% |
| G_HI | exact | −2.73 | 0 | +14.9% | −5.8% |
| G_HI and ν_eff, least squares in log | exact | −3.04 | −6.3% | +8.7% | −0.7% |
| ν_eff | edges | −3.70 | −21.3% | 0 | +3.5% |
| G_HI | edges | −2.48 | 0 | +21.7% | −15.2% |
| G_HI and ν_eff, least squares in log | edges | −2.91 | −8.7% | +11.0% | −7.7% |

Hardening mismatch (`hardening_mismatch`). A photons band loses its mean photon energy ⟨hν⟩_b per absorption; the
products receive the mean absorbed photon energy E_abs = 13.6 eV + G_HI. Since σ ∝ E⁻³ takes the soft photons,
E_abs < ⟨hν⟩_b, and the energy difference is lost from the band without going anywhere. The true spectrum hardens as
it is absorbed; a fixed slope cannot.

| Band | ⟨hν⟩_b | E_abs | (⟨hν⟩ − E_abs)/⟨hν⟩ | same, 4e4 K blackbody in the band |
|---|---|---|---|---|
| EUV_H | 17.46 | 16.30 | 6.6% | 6.2% |
| EUV_He | 29.22 | 27.59 | 5.6% | 4.3% |
| single 13.6-500 eV band, ν_eff fitted | 18.64 | 16.09 | 13.7% | 10.9% |

Most of the mismatch belongs to the blackbody itself (absorbed photons are softer than the band's mean), not to the
fixed slope, which adds 0.4 and 1.3 points. The split halves it against one band. Removing it needs a second moment
per band (photon number and energy, i.e. a free slope); not done, pending a decision.

Band convergence, after the paper's §4.1.2: χ ∝ E^−2 with its Kirchhoff emission, N bands log-spaced over 0.01-10 eV
(plus two closing bands below and above), against 512 bands (256 vs 512 bands differ by < 1e-4).

| N | relaxation, `kirchhoff="band"` (default): max error / at equilibrium | relaxation, `kirchhoff="planck"` | absorbed fraction of a 6000 K blackbody at τ(1 eV) = 1 |
|---|---|---|---|
| 4 | 1.2e-1 / 5e-7 | 5.6e-2 / 3.4e-2 | +43% |
| 8 | 4.8e-2 / 2e-7 | 2.5e-2 / 9.2e-4 | +9.4% |
| 16 | 1.2e-2 / 6e-8 | 6.5e-3 / 3e-6 | +2.4% |
| 32 | | | +0.6% |

The relaxation error is |T − T_ref| / (T0 − T_eq), for gas at 3000 K with C_v = a T0³ and no radiation initially;
the reference is 512 bands with the default emission. All columns converge at second order. Kirchhoff per band
reaches the exact equilibrium at any band count, at twice the transient error. About half the absorbed-fraction error
is the fixed slope. The rest is one mean opacity inside an exponential: weighting χ_E,b by the true blackbody gives
+19% / +4.4% / +1.2% / +0.3%.

The legacy RT model is not band-diagonal (`band_structure`): EUV→ONIR (donation), FUV/NUV/ONIR→IR (the copied dust
absorption), and cooling routed into NUV/IR whose rate reads other bands (f_recNUV, G_0, T_bg, Compton off every band).

## Decisions (2026-10-06)

1. **EUV: split at 24.59 eV.** `STARFORGE_RT` = EUV_H, EUV_He (photons), FUV, NUV, ONIR (energy), IR (tracked). Each
   ionizing slope matches the 4e4 K blackbody's G_HI and ν_eff in its sub-band jointly; all residuals < 1%.
2. **In-band opacity: exact quadrature by default.** The projection runs once at codegen, so quadrature is free, and
   it keeps a split band consistent with its parent for any σ(E). `opacity="edges"` (the paper's PPL) remains an
   option.
3. **Thermal emission into ppl bands: Kirchhoff per band by default** (χ_B,b := χ_E,b), for the exact
   radiation-matter equilibrium at any band count (the dust-IR loop equilibrates in thick cells). The exact Planck
   mean remains as `kirchhoff="planck"`.
4. **Photon versus energy bookkeeping: photons exact, hardening booked.** The per-sub-band mismatch is 6.6% (EUV_H)
   and 5.6% (EUV_He) of ⟨hν⟩_b, against 13.7% for one band (table above). Photons bands conserve photon number
   exactly (ionizations set Strömgren radii). The deferred energy, ⟨hν⟩_b − E_abs per absorption (the energy the true
   spectrum would keep as hardening and deliver downstream), is an explicit term in the energy ledger: an Output, so
   that closed-box ledgers close with it. No second moment. Check: HII_region with the ionizing range split into ~8
   sub-bands as a converged reference for the 2-band ionization-front temperature and radius.

## Open choices

1. **Emitted photons in a photons band:** f_b P / ⟨hν⟩_b (conserves energy) or `photons_per_eV` (conserves photons).
   Both are provided.
2. **Tables.** `TTable.expr` inlines a log-log `piecewise_linear` (default grid 1-1e8 K at 0.025 dex, 321 points).
   For many tables, a 1-D `Table` atom in `jaco_tables.hdf5` would be leaner.
3. **Flux mean.** Under a fixed slope χ_F = χ_E + scattering. There is no Rosseland or diffusion-limit correction
   (the paper's χ_F,diff).

## The legacy model on the band layer

`starforge_legacy_RT` and `starforge_legacy_RT_EUV` take their bands from `GIZMO_STARFORGE` and GIZMO's coefficients
from one `Projector`, `radiation.GIZMO_RT`, whose absorbers and emitters declare no spectra and whose overrides are
GIZMO's constants. The generated code is byte-identical (golden hashes unchanged).

| GIZMO quantity | Now |
|---|---|
| band species, order, docs, edges (13.6, 500 in the IR tail and G_LW), T_rad | `GIZMO_STARFORGE` (`Band.doc`, `Band.temperature`) |
| σ_HI, ε_HI (rt_ion_sigma_HI, rt_ion_G_HI) | override (H photoionization, EUV): `sigma_N`, `heat` = ε_HI / eV (exact Rational eV, so ε_HI stays a bare symbol) |
| hν_EUV (rt_nu_eff_eV) | band override EUV `hnu`; `band_energy_eV` uses it for any photons band |
| DUST_BAND_OPACITY (720/480/180), albedo 1/2 | override (dust, FUV/NUV/ONIR): `chi_F` = κ, `chi_E` = κ/2; `DUST_BAND_OPACITY` is now derived |
| IR dust and gas opacities (rt_kappa_adaptive_IR_band) | override (dust, IR) and (gas, IR): `chi_E` = the existing expressions |
| dust IR emission opacity κ(T_d, T_d) | the (dust, IR) `chi_E` at T_rad = T_d (`expr("chi_E", Td)`: Kirchhoff on the tracked band) |
| all dust emission into IR | `RoutedEmission` "Dust emission into photon_IR", fraction 1 in IR |
| Compton temperatures (2340 hν_EUV, 24400, 12000, 2800 K, T_rad) | band override `T_compton` |
| CoolingRate's routing of cooling into NUV/IR (f_recNUV, the 1e5 K free-free split, NUV→IR above T_rad = 1e4 K) | `RoutedEmission` per routed process in `starforge_legacy_RT.COOLING_ROUTES`, fraction overrides |

Kept in `radiation.py`, as GIZMO's composition laws and scheme: rt_kappa's max(0.02 + 0.35 x_e, κ max(Z floor, Z f_d))
and the PE band's Z floor; the start-of-step dust temperature in the kick (an `xreplace` of the override); the kick
factors; the donation to ONIR; the IR double-count copy; the gas IR absorption's kick share; 2.16e-35 (Compton).

Gaps found, i.e. what the band layer could not express from spectra and what had to be added:

1. **Tracked bands are whole blackbodies in GIZMO.** GIZMO's IR band holds a full blackbody at T_rad; its edges are
   nominal, and its tail above 13.6 eV photoionizes and above 11.2 eV adds to G_LW. The band layer confines a band's
   spectrum to its edges, so a band cannot couple outside them. The tail coupling is hand-written with the EUV band's
   coefficients and GIZMO's blackbody_lum_frac fit.
2. **Products into a band.** An absorber's products go only to heat rows. GIZMO donates the ionizing photon's hν_EUV
   to ONIR (on top of ε_HI to the gas) and copies the dust-absorbed FUV/NUV/ONIR energy into IR a second time. Both
   create energy and stay explicit processes; a conserving model would express reprocessing as absorption into a
   reservoir plus emission.
3. **Non-additive opacity laws and composition scaling.** Absorbers add; GIZMO takes max(gas floor, dust), floors
   Z f_dust in one band, and scales by dust survival f_d(T_d). The band layer holds the κ per unit absorber; the law
   that combines and scales them stays in the model.
4. **State-dependent opacities.** GIZMO's IR opacities are tables in (T_rad, T_d) with composition zones, plus a gas
   opacity in x_e, x_H+, ρ, u. A band layer κ(E) has no matter-state dependence; these are override expressions.
   Added: `AbsorptionCoefficients.expr(name, T)` evaluates a tracked band's coefficient at another temperature
   (Kirchhoff emission at T_d).
5. **Band-level quantities.** Only (process, band) overrides existed. Added `band_overrides` for ⟨hν⟩ (a runtime
   parameter here) and a Compton temperature, which the layer lacked entirely (now projected as ⟨E⟩_u / 4k).
   GIZMO's constants sit below the fixed-slope projections: FUV 24400 vs 30600 K, NUV 12000 vs 15700 K, ONIR 2800 vs
   4100 K; for the tracked IR band, T_rad vs 0.998/0.958/0.71 T_rad at T_rad = 10/100/1000 K (the band's upper edge
   truncates a hot blackbody).
6. **Declarations without spectra.** GIZMO's band constants have no σ(E) behind them. Added `Absorber(cross_section=
   None)` (coefficients only from overrides; those not given are None) and `RoutedEmission` (fractions only).
7. **Units.** The layer's energies are eV, GIZMO's ε_HI is erg; an exact Rational eV keeps ε_HI a bare symbol through
   the round trip. Coefficients that are runtime parameters (σ_HI, ε_HI, hν_EUV, computed by GIZMO from T_eff) enter
   as symbol overrides: the layer has no notion of a coefficient projected at run time.

## Next steps

1. Done: the legacy RT models on `GIZMO_STARFORGE` with overrides (above). Next, `starforge_RT` on `STARFORGE_RT`
   without overrides.
2. Process factories: an absorber becomes one `Reaction`/`Transfer` per overlapping band (`row_factors` rsol, and
   ⟨hν⟩_b or E_abs by the band's unit). An emitter becomes rows f_b × its heat on its source process (an `Output`
   while the band is not solved).
3. Host interface macros from the band set: `JACO_N_BANDS`, `JACO_BAND_<NAME>` (index), `JACO_BAND_<NAME>_E_LO/E_HI`,
   `JACO_BAND_<NAME>_HNU` (⟨hν⟩; a function of T_rad for tracked bands), per-absorber `_KAPPA_E`/`_KAPPA_F` for the
   transport opacity, and `JACO_BAND_JACOBIAN_DIAGONAL` from `band_structure`. GIZMO static-asserts its
   `RT_FREQ_BIN_*` against them.
4. Solver: with a diagonal band block, eliminate the bands by a Schur complement onto the matter variables, at O(N_b)
   per Newton step (the paper's §3).
