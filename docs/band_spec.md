# Radiation bands as a spec (`jaco.bands`)

Status: core layer, phase P2. Nothing in `jaco.models` uses it yet; the generated code is unchanged.

Processes declare which photon energies they interact with: a cross-section sigma(E) or opacity kappa(E) for an
absorber, a spectrum j(E) for an emitter. A `BandSet` turns the declarations into per-band coefficients at
code-generation time. Refining a band (splitting FUV at 11.2 eV, EUV at the He thresholds) changes only the band set.

## API

```python
from jaco.bands import Band, BandSet, Absorber, ThermalEmission, Line, Continuum, Projector, SIGMA_HI, fit_slope

bands = BandSet([
    Band("EUV", 13.6, 500.0, unit="photons", slope=-3.70),  # fixed-slope spectrum u_E ~ E^slope
    Band("FUV", 8.0, 13.6),                                    # default: energy, slope -1 (E u_E = const)
    Band("IR", 0.001, 0.4133, shape="tracked"),                # dilute blackbody at the per-cell T_rad
])
bands = bands.split("FUV", [11.2], names=["FUV_lo", "LW"])     # refinement: only the spec changes

proj = Projector(bands,
                 absorbers=[Absorber("H photoionization", SIGMA_HI, absorber="H", E_th=13.6),
                            Absorber("dust", kappa_dust, heat_to="dust heat", scattering=kappa_scat)],
                 emitters=[ThermalEmission("dust emission", kappa_dust, source="dust heat"), Line("Halpha", 1.89)],
                 overrides={("H photoionization", "EUV"): {"sigma_N": sp.Symbol("sigma_HI")}})
c = proj.absorption("H photoionization", "EUV")   # AbsorptionCoefficients, or None if no overlap
e = proj.emission("dust emission")                 # EmissionFractions
proj.mean_photon_energy("EUV")                     # <h nu>_b
```

| Object | What it is |
|---|---|
| `Band(name, E_lo, E_hi, unit, shape, slope)` | `unit` "photons" or "energy" (what the species `photon_<name>` counts); `shape` "ppl" (fixed slope) or "tracked" (blackbody at T_rad). |
| `BandSet` | Ordered, non-overlapping, gaps allowed; `index`, `band_at(E)`, `gaps`, `split`, `replace`, `declarations()` (radiation `Species`). |
| `Absorber(name, cross_section, absorber, E_th, heat_to, heat_yield, remainder_to, scattering)` | Per absorbed photon: `E_th` to chemistry, `Y(E)(E - E_th)` to `heat_to`, the rest to `remainder_to`. Scattering enters the flux mean only. |
| `Line(name, E0, source)`, `Continuum(name, j, source, E_range)`, `ThermalEmission(name, kappa, source)` | Emitters; `j(E)` or `j(E, T)`; thermal is `4 pi kappa(E) B_E(T)`. `source` is the row losing the energy. |
| `PowerLaw`, `Spectral`, `Blackbody`, `planck`, `SIGMA_HI/HEI/HEII` | Spectral functions of E [eV]; Verner et al. (1996) fits. A `PowerLaw` gets closed forms. |
| `Projector(bands, absorbers, emitters, overrides, opacity, T_grid)` | Projections, cached. Ppl bands give floats, tracked bands and emitter temperatures give `TTable`s (log-log in T). `.expr(...)`/`symbolic` give sympy (`Float`, or a `piecewise_linear` in T). |
| `fit_slope(band, reference, targets, absorbers, opacity, weights, report)` | Slope of a ppl band matching chosen band means of a reference spectrum. |
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

`opacity="edges"` (the default, the paper's PPL) replaces σ and σ_s inside a ppl band by the power law through their
values at the band's edges, then uses the closed forms; `"exact"` projects the declared function. Tracked bands always
use the declared function.

Emission: `fraction[b]` = ∫_b j / ∫ j, `photons_per_eV[b]` = ∫_b j/E / ∫ j, and `escape` (energy in no band) is
integrated over the gaps independently, so Σ fraction + escape = 1 is a check, not an identity. A thermal emitter also
gets per band `planck_mean` χ_B,b(T) = ∫_b κB / ∫_b B (computed with the band's lower-edge Wien factor taken out, so it
stays finite where B underflows), `band_planck` B_b(T), and the total `kappa_planck`. A band's thermal emission is
4π χ_B,b B_b.

Overrides: `{(process, band): {coefficient: value}}`, numbers or sympy expressions, replace projected values; an
override may couple an absorber to a band it does not overlap; an overridden emission fraction sets escape to
1 − Σ fractions. Results record `overridden`.

## Results (tests/test_bands.py, tests/test_band_convergence.py)

GIZMO's ionizing band (rt_get_sigma, T_eff = 4e4 K, 13.6-500 eV). Replicating its rectangle rule and hydrogenic
cross-section reproduces its printout, σ_HI = 3.368e-18 cm², G_HI = 2.99 eV, ν_eff = 18.6 eV. Quadrature with
Verner (1996): 3.370e-18, 3.012, 18.635.

One slope cannot match both G_HI and ν_eff: a power-law spectrum fixes ν_eff by the slope alone, and G_HI
needs a slope a full unit shallower.

| Fit | opacity | slope | G_HI | ν_eff | σ_HI |
|---|---|---|---|---|---|
| ν_eff | exact | −3.70 | −17.4% | 0 | +8.5% |
| G_HI | exact | −2.73 | 0 | +14.9% | −5.8% |
| G_HI and ν_eff, least squares in log | exact | −3.04 | −6.3% | +8.7% | −0.7% |
| ν_eff | edges | −3.70 | −21.3% | 0 | +3.5% |
| G_HI | edges | −2.48 | 0 | +21.7% | −15.2% |
| G_HI and ν_eff, least squares in log | edges | −2.91 | −8.7% | +11.0% | −7.7% |
| two bands cut at 24.59 eV, each fitted to its own ν_eff | edges / exact | −2.27, −6.61 | −2.3% / −1.9% | 0 | +0.1% / +0.7% |

Band convergence, after the paper's §4.1.2: χ ∝ E^−2 with its Kirchhoff emission, N bands log-spaced over 0.01-10 eV
(plus two closing bands below and above), against 512 bands (256 vs 512 bands differ by < 1e-4).

| N | relaxation, emission at Planck mean: max error / at equilibrium | relaxation, emission at χ_E,b (Kirchhoff per band) | absorbed fraction of a 6000 K blackbody at τ(1 eV) = 1 |
|---|---|---|---|
| 4 | 5.6e-2 / 3.4e-2 | 1.2e-1 / 4e-6 | +43% |
| 8 | 2.5e-2 / 9.2e-4 | 4.8e-2 / 3e-6 | +9.4% |
| 16 | 6.5e-3 / 3e-7 | 1.2e-2 / 3e-6 | +2.4% |
| 32 | | | +0.6% |

The relaxation error is |T − T_ref| / (T0 − T_eq), for gas at 3000 K with C_v = a T0³ and no radiation initially.
All columns converge at second order. About half the absorbed-fraction error is the fixed slope. The rest is one mean
opacity inside an exponential: weighting χ_E,b by the true blackbody gives +19% / +4.4% / +1.2% / +0.3%.

The legacy RT model is not band-diagonal (`band_structure`): EUV→ONIR (donation), FUV/NUV/ONIR→IR (the copied dust
absorption), and cooling routed into NUV/IR whose rate reads other bands (f_recNUV, G_0, T_bg, Compton off every band).

## Open choices

1. **EUV.** Options: (a) match ν_eff, so the stellar photon count and energy agree, with G_HI 17-21% low;
   (b) match G_HI; (c) the joint fit; (d) split at 24.59 eV, which matches all three within 2% and is the cut He
   photoionization needs anyway. I recommend (d); the spec cost is one band.
2. **Opacity inside a ppl band: "edges" or "exact".** "edges" follows the decided design and the paper. It is
   split-invariant only for true power laws, and on a wide band it is not small: on the single EUV band at fixed
   slope it lowers σ_N and G by ~5% against "exact". QUOKKA needs edge values because it evaluates opacities at run
   time. Here the projection runs once at codegen, where quadrature costs nothing. Recommend "exact" as the default.
3. **Thermal emission: exact Planck mean, or χ_B,b := χ_E,b.** The exact Planck mean (current) gets the emitted power
   right and halves the transient error, but its radiation-matter equilibrium is wrong at coarse bands (3.4% at 4).
   Kirchhoff per band (the paper) gives the exact equilibrium, but its total cooling then depends on the band count.
   For the dust-IR loop, where T_d and T_rad equilibrate in thick cells, the second may matter more. Both are one
   switch.
4. **Photon versus energy bookkeeping.** A fixed slope cannot harden. A photons band loses ⟨hν⟩_b of energy per
   absorption while the products take E_abs (EUV: 18.6 vs 16.0-16.1 eV). An energy band has the converse error on
   photon number. GIZMO has the same mismatch (G_HI vs ν_eff − 13.6). Accept it, or say which unit each band conserves.
5. **Emitted photons in a photons band:** f_b P / ⟨hν⟩_b (conserves energy) or `photons_per_eV` (conserves photons).
   Both are provided.
6. **Tables.** `TTable.expr` inlines a log-log `piecewise_linear` (default grid 1-1e8 K at 0.025 dex, 321 points).
   For many tables, a 1-D `Table` atom in `jaco_tables.hdf5` would be leaner.
7. **Flux mean.** Under a fixed slope χ_F = χ_E + scattering. There is no Rosseland or diffusion-limit correction
   (the paper's χ_F,diff).

## Next steps

1. Port the legacy RT model onto a 5-band spec (EUV, FUV, NUV, ONIR, IR) with overrides carrying GIZMO's constants:
   `DUST_BAND_OPACITY` per band, and σ_HI, ε_HI, hν_EUV as runtime `Parameter` symbols. Check it against the golden
   hashes.
2. Process factories: an absorber becomes one `Reaction`/`Transfer` per overlapping band (`row_factors` rsol, and
   ⟨hν⟩_b or E_abs by the band's unit). An emitter becomes rows f_b × its heat on its source process (an `Output`
   while the band is not solved).
3. Host interface macros from the band set: `JACO_N_BANDS`, `JACO_BAND_<NAME>` (index), `JACO_BAND_<NAME>_E_LO/E_HI`,
   `JACO_BAND_<NAME>_HNU` (⟨hν⟩; a function of T_rad for tracked bands), per-absorber `_KAPPA_E`/`_KAPPA_F` for the
   transport opacity, and `JACO_BAND_JACOBIAN_DIAGONAL` from `band_structure`. GIZMO static-asserts its
   `RT_FREQ_BIN_*` against them.
4. Solver: with a diagonal band block, eliminate the bands by a Schur complement onto the matter variables, at O(N_b)
   per Newton step (the paper's §3).
