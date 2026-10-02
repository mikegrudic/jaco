# GIZMO legacy ISM microphysics vs the jaco `starforge` model: sync audit

Date: 2026-10-02. GIZMO: `gizmo_jaco_dev` = `origin/starforge_dev` at `8da21be3`. jaco: branch `model_sync`, based on
`f6d0c76` (`origin/slop_experiments`). The scope is the physics compiled into GIZMO's non-RT cooling tests
(`test/gmc_cooling`, `test/two_temperature`). The Kim+23 nebular term is RT-only but was added anyway.

Units: GIZMO's `CoolingRate()` returns rates per n_H^2. Every comparison below is of volumetric rates
(erg cm^-3 s^-1), i.e. GIZMO's value times n_H^2.

## 1. Flags enabled in the non-RT configs

Expanded with `cpp -dM` over `declarations/precompiler_logic.h`, using the GIZMO_config.h that
`prepare-config.perl` writes for each Config.sh.

- **Cooling/EOS/chemistry flags that are on:** `COOLING`, `COOL_LOW_TEMPERATURES`, `COOL_METAL_LINES_BY_SPECIES`,
  `COOL_MOLECFRAC_NONEQM`, `SIMPLE_STEADYSTATE_CHEMISTRY`, `EOS_SUBSTELLAR_ISM`, `EOS_GENERAL`, `EOS_GAMMA_VARIABLE`,
  `METALS` (NUM_METAL_SPECIES = 11), `GALSF`, `GALSF_FB_FIRE_STELLAREVOLUTION=3`, `RT_ISRF_BACKGROUND`,
  `RT_USE_TREECOL_FOR_NH=6`, `SINGLE_STAR_SINK_DYNAMICS`, `MAGNETIC` (which also enables Compton-like synchrotron losses),
  `CONDUCTION_SPITZER`, `VISCOSITY_BRAGINSKII`, `OUTPUT_MOLECULAR_FRACTION`.
- **Off in these configs:** `RT_CHEM_PHOTOION`, `RT_PHOTOELECTRIC`, `RT_INFRARED`, `RT_NUV`, `RADTRANSFER`,
  `SINGLE_STAR_FB_*`, `GALSF_FB_FIRE_RT_*`, `GALSF_ISMDUSTCHEM_*`, `COSMIC_RAY_FLUID`, `CHIMES`, `COOL_LOWTEMP_THIN_ONLY`,
  `COOLING_OPERATOR_SPLIT`.
- `SINGLE_STAR_STARFORGE_DEFAULTS` alone does not turn on the FB/RT flags. Those come from the hybrid block, which
  these configs do not use.
- `TWO_TEMPERATURE_PLASMA=1` is referenced nowhere in this tree's sources (its test skips when it is not
  implemented). The two configs therefore have identical cooling physics; `gmc_cooling` only adds `OPENMP_GPU_OFFLOAD`.
- Run-time settings: `ComovingIntegrationOn=0`, `InterstellarRadiationFieldStrength=1`, and `Redshift_RT_Background`
  unset (0). With TREECOOL present, `J_UV != 0`, so GIZMO's UVB and FIRE-3 metal-line branches are both active.

## 2. GIZMO cooling/EOS commits since the model was written (`4bdd08bc..HEAD -- cooling/ eos/`)

| commit | change | active in non-RT tests? | jaco action |
|---|---|---|---|
| dde244c1, bad73174 | Kim+23 nebular cooling and its CIE taper | no (`RT_CHEM_PHOTOION && METALS`) | **added**, gated by `f_neb` (9da5dfb) |
| b06e1503 | nebular emission routed into the NUV band | no (`RT_NUV`) | n/a: the model has no radiation bands |
| 12ec77bd | GALSF_ISMDUSTCHEM sync; in cooling.cc only adds `&& METALS` to the nebular guard | no | none |
| e82d7287, 174d05aa, 1eded91d | H2 rotational partition series: iteration cap and icpx miscompile workaround | yes (numerical only) | none: the math is unchanged, and jaco's fit matches GIZMO to <0.3% in E/kT and <1% in c_v over 5 K to 1e4 K |
| 70663252 | `gas_dust_heating_coeff` moved to eos.cc | yes (code motion) | none: identical |
| d4e01a38 | cache Gamma from the freshly computed T | yes (GIZMO-internal EOS) | none |
| 05f0bf48, 85eb89ca, 582f0d92, 693aba19, ee167170 | buffer sizes, guards, vec3, conduction NaN guard, jaco API fix | n/a | none |

## 3. Gap table

The Status column uses: **same** (identical expression and coefficients); **FIXED <hash>** (changed in jaco on
`model_sync`); **differs** or **MISSING** (reported only, see section 6). File references are cooling.cc unless stated.

### 3.1 H/He ionization balance and primordial channels

| term | GIZMO | jaco | status | difference |
|---|---|---|---|---|
| H0 collisional excitation | 1475, 991 | `starforge/line_cooling.py` `LineCoolingSimple("H")` | same | KWH96 |
| He+ collisional excitation | 1476, 992 | `LineCoolingSimple("He+")` | same | |
| collisional ionization rates (H, He, He+) | 1487-1489 | `processes/ionization.py` | same | |
| collisional ionization cooling | 995-997 | `Ionization.heat` | same | 13.6/24.59/54.42 eV |
| recombination rates (VF96 radiative, He+ dielectronic) | 1482-1486 | `processes/recombination.py` | same | |
| H+ recombination cooling | 1001 | `recombination.py` | same | 1.036e-16 T alpha |
| He+ recombination cooling | 1002, 1004 | `recombination.py` | **FIXED 39bfa73** | was 1.55e-26 T^-0.3647 (exponent sign; ~800x low at 1e4 K); the dielectronic term 6.526e-11 alpha_d was missing (4e-23 at 1e5 K) |
| He++ recombination cooling | 1003 | `recombination.py` | **FIXED 39bfa73** | was overwritten with 4x the H+ value, which is 0.57-0.77x GIZMO's over 1e4-3e5 K; now uses its own VF96 rate |
| free-free | 1478 (1.43e-27), 1300 (1.42e-27 above Tmax) | `processes/freefree_emission.py` (1.42e-27) | differs 0.7% | GIZMO itself uses two values; left alone |
| H/He ionization state | 901-930: KWH equilibrium including UVB photoionization x Rahmati shieldfac, no CR ionization of H | steady state of collisional + CR (`cosmic_ray_ionization("H")`) + recombination, no UVB | **differs** (structural) | UVB missing in jaco; jaco has CR ionization of H that GIZMO's H balance lacks |
| electrons from metals | 940-948; `simple_chemistry.cc` 97-170 (heavy ions, alkali, C+ via Gong17 balance, O+ = x_H+ x 3.2e-4 Z_O, molecular ions) | C+ only, as a fixed species with the Tielens-form f_C+ | **MISSING** (structural) | would need new fixed species |
| UVB photoheating | 1218 | none | **MISSING** (structural) | needs Gamma/epsilon (x shieldfac) parameters |
| UVB self-shielding (Rahmati+12) | 2473-2487 | none | **MISSING** | only meaningful together with the UVB |

### 3.2 Nebular cooling (Kim, Gong, Kim & Ostriker 2023, Eq. 47)

| term | GIZMO | jaco | status | difference |
|---|---|---|---|---|
| nebular forbidden-line cooling, CIE taper, collisional de-excitation | 1093-1106 (`RT_CHEM_PHOTOION && METALS`, so not in these builds) | `starforge/nebular_cooling.py` | **ADDED 9da5dfb** | matches GIZMO to 1e-7 for T >= 3e3 K. GIZMO's hard 2e3 K window edge is a logistic 0.01 dex wide here (the fit is still 8% of its peak at 2e3 K, so a hard step would be a real discontinuity). The fit is clamped to [2e3, 5e4] K. n_e has a 1e-20 cm^-3 floor inside the de-excitation power, to keep the derivative finite. jaco's usual C_2 clumping factor applies. |

Metallicity scaling: `Z_d`. GIZMO's term uses `Metallicity[0]/SolarAbundances[0]`, and `gizmo_to_jaco` sets
`Z_d = max(1e-4, Z/Z_sun)`. That is the same ratio, and differs only below 1e-4 Z_sun where the term is negligible.
In jaco, dust processes carry their depletion through `Z_d*f_d`, and this term uses `Z_d` alone. If `Z_d` is ever
made dust-specific, switch this term to `x_O,tot / x_solar("O")`.

### 3.3 Tabulated metal lines (`COOL_METAL_LINES_BY_SPECIES`, FIRE-3 branch)

| term | GIZMO | jaco | status | difference |
|---|---|---|---|---|
| element abundance multiplying each table | 1118, `GetCoolingRateWSpecies` 1866 (total Z_k) | `starforge/metal_line_cooling.py` | **FIXED 1dc26a2** | C used the network's neutral C (x_C,tot - x_C+ - x_CO) and O used x_O,tot - x_CO. Carbon table cooling vanished wherever C+ dominates. The tables are per total element abundance. |
| taper below the tables' 100 K edge | 1124 | | **FIXED 1dc26a2** | jaco held the 100 K value (clamped lookup); now x exp(-(2-log T)^2/0.1) |
| CMB bath (T-T_cmb)/(T+T_cmb), applied to the summed rate only when it is positive | 1125-1127, 1272 | | **FIXED 5bc57ea, 90cfcf3** | the tables carry UVB photoheating (negative entries up to log T ~ 6-7). GIZMO adds net negative LambdaMetal to Heat without the CMB factor. |
| applied at all T, including negative entries | 1112 (`J_UV != 0` for FIRE>2; `logT > 4` only for FIRE<=2) | applied at all T | same | not a discrepancy for these configs (coordinator item 2). If `J_UV == 0` (no TREECOOL, or z beyond the table), GIZMO turns metal lines off entirely; jaco has no switch for that. |
| per-element solar normalization | GIZMO: Z_k/(Z_sun,k x 0.0127/Z_sun) | conversion script: x_k/x_k,Wiersma | **differs** (reported) | jaco/GIZMO at solar: C 1.07, N 0.78, O 0.98, Ne 0.84, Mg 1.14, Si 0.92, S 0.70, Ca 0.94, Fe 1.09. jaco's normalization is the physically consistent one; GIZMO approximates the tables' solar pattern by its own pattern rescaled to Z = 0.0127. |
| table interpolation | bilinear in (log n_H, log T), clamped | C header `jaco_table.h`: same | same in generated code | Python-side `TableInterp2D._evaluate` interpolates linearly in linear n_H and T and extrapolates (`RegularGridInterpolator` on 10**log axes), so Python solves disagree with GIZMO and with the C code |

### 3.4 Low-temperature block (`COOL_LOW_TEMPERATURES`, detailed branch for METAL_LINES && FIRE>2)

| term | GIZMO | jaco | status | difference |
|---|---|---|---|---|
| C+ fine structure (e-, H colliders) | 1168 | `starforge/line_cooling.py` | **FIXED 669aedc** | per-ion coefficients divided by `SolarAbundances.x("C")` = 2.11e-4 (computed as f/(1-f)/A, not f/(A X_H)); GIZMO's Z_C scaling corresponds to 2.95e-4, so C+ cooling was 1.40x high. Still different: GIZMO multiplies by the neutral fraction x_H0, so there is no C+ cooling in ionized gas (effect <1% at T > 1e4 K); GIZMO's H collider is all neutral nuclei including H2, jaco's is atomic H only; e- coefficient 4890e-27 vs 4888e-27. |
| [CI] 609 um, 2.08e-29 exp(-23.6/T) Z_C | 1168 | none | **MISSING** (reported) | GIZMO weights it by the C+ fraction (physically it should be neutral C). It dominates C+ below ~20 K. |
| C+ fraction 1/(1+(n/340 G0)^2/sqrt T) | 1166-1167 | `starforge.py` `x_Cplus_expr` | same form | jaco floors T at 10 K. G0 differs through the interface (section 5). |
| CO cooling | 1169-1175 (HM79 recalibrated, x (1-f_C+), LVG cap 4.42e-28 dv/dr T^4) | `starforge/CO_cooling.py` (Whitworth & Jaffa 2018) | **FIXED b73afbf** (sign); otherwise differs | the rate was passed as a positive heat coefficient, so CO *heated* the gas. Different model otherwise: WJ18's high-density limit equals GIZMO's LVG cap when x_H2 = 0.5; at low density jaco is ~0.6x GIZMO in fully molecular solar gas. jaco's x_CO = x_C,tot (1-f_C+) x_H2 caps at half the carbon. |
| H2 cooling (GA08 thin, HM79 thick, n/n_crit) | 1178-1196 | `starforge/H2_cooling.py` | same fits; colliders differ | GIZMO weights colliders by mass fractions (X_H, Y_He) instead of n_X/n_H: its H2-He term is ~2.9x too strong and the others 0.7x; jaco's weighting is the correct one. jaco caps log T at 4 and T/1e3 at 10. |
| HD abundance | 1191: min(0.00126 x_H2, 4e-5 x_H0) | `starforge.py`: 2.527e-5 x_H2 | **differs** (reported) | GIZMO's HD/H is 3.2x jaco's in molecular gas. jaco's own `n_HD_prescription()` in H2_cooling.py matches GIZMO but is unused. "All D in HD" would be 2.527e-5 x 2 x_H2. |
| block cut-off above 10^4.5 K: exp(-((log T-4.5)/0.2)^2) on molecular, C+, CO and dust | 1157, 1199, 1203 | | **FIXED 5bc57ea** | GIZMO also skips the block above 10^5.3 K, where jaco's factor is already < 1.2e-7 |
| CMB bath (T-T_cmb)/(T+T_cmb) on LambdaMol (H2/HD, C+, CO) | 1207 | | **FIXED 5bc57ea** | T_cmb = 2.73 (1+z) from the existing `z` parameter |
| gas-dust coupling coefficient | eos.cc 248 | `gas_dust_collisions.py` | same | jaco adds C_2 and an a_grain factor (= 1) |
| gas-dust sputtering cut-off above 3e5 K | 1202 | | **FIXED 5bc57ea** | |
| dust temperature | 2168 (`SINGLE_STAR_SINK_DYNAMICS` branch: `rt_eqm_dust_temp` with CMB + ISRF absorption, optical extinction and gas coupling) | parameter `Td` | **MISSING** (structural) | GIZMO's interface passes Td = 10 K (jaco.cc 364) |
| dust-to-metals ratio | eos.cc 217 | `f_d` | same | equals 1 under `SINGLE_STAR_SINK_DYNAMICS`, and GIZMO passes f_d = 1 |
| low-T fallback fit (2.8958629e-26 ...) | 1149-1156 | none | n/a | replaced by the detailed branch in these configs |
| optically-thick cap on net cooling (Rafikov 2007 photosphere) | 1392-1418 | none | **MISSING** (structural) | a cap on the total rate, plus an opacity model |

### 3.5 Heating

| term | GIZMO | jaco | status | difference |
|---|---|---|---|---|
| photoelectric (BT94/Wolfire) | 1281-1293 | `photoelectric_heating.py` | same expression; **FIXED b2372aa** | GIZMO's `T < 1e6` condition was missing. The efficiency fit grows as T^0.7, which heated 1e7 K gas at roughly 20% of its cooling rate (estimate at n_H ~ 1e-3-1e-2). The step is smoothed over ~0.05 dex. G0 differs through the interface. |
| CR heating | 1269, `cosmic_ray_utilities.cc` 1829 (Guo & Oh form, x e_CR/(0.01+n_H)) | 20 eV per CR ionization of atomic H (`cosmic_ray_ionization.py`) | **differs** (reported) | GIZMO gives ~2.5e-28 erg/s per H in any neutral gas; jaco gives ~5e-28 per H atom and **zero in molecular gas** (the CR H2 channels carry no heat) |
| CR ionization rate | `cosmic_ray_utilities.cc` 1775-1823 | `starforge/symbols.py` | **FIXED 1fea3dd** | attenuation was Min(1, 1e21/N exp(-N/1e24)). GIZMO's is Sigma_0/Sigma with Sigma_0 = 2.23e-3 g cm^-2 (1.33e21 in the Sigma/m_p units GIZMO passes N_H in) and an exponential cut at 100 g cm^-2 (6e25, not 1e24). Also added the radionuclide floor 1e-21 Z/Z_sun + 1e-19 Fe/Fe_sun. |
| Compton off the CMB | 2315-2325 | `starforge/compton.py` | **FIXED 72fe6db** | was 5.41e-36 n_e T (1+z)^4; now 2.16e-35 x 0.262 (1+z)^4 n_e (T - T_cmb). The generic `inv_compton_cooling` is unchanged. |
| Compton off the UVB, the fixed MW ISRF (0.31 eV/cc at 30 K, 0.66 at 5800 K), and synchrotron (`MAGNETIC`) | 2321, 2371, 2387 | none | **MISSING** | needs radiation and magnetic energy-density parameters. Negligible except in hot, diffuse gas. |
| hydro (PdV) work | 1430 | `pdv_work` parameter | same | |

### 3.6 H2 chemistry (GIZMO's single-variable implicit update vs jaco's network)

| term | GIZMO | jaco | status | difference |
|---|---|---|---|---|
| overall rate convention | 1969-2148 | `h2_chemistry/*` | differs (GIZMO-side) | GIZMO evolves the mass fraction f = 2n_H2/n_H0 with per-molecule rates times 1/2, so its formation and destruction are both half speed: same equilibrium, 2x relaxation time. jaco is per molecule. |
| dust formation (GJ07/HM79) | 2052 | `grain_formation.py` | same coefficient | collider product n_HI^2 vs GIZMO's n_H0 n_HI (/2, as above) |
| H- channel k1, k2, k5, k17, photodetachment R51 | 2061-2072 | `radiative_association.py`, `associative_detachment.py`, `mutual_neutralization.py`, `H2_chemistry.py` | same | |
| H- collisional detachment k15 (e-), k16 (H) | 2060-2069 | `collisional_detachment.py` | **FIXED 7665512** | the T_eV fits were evaluated at ln(T/8.617e-5) = ln T + 9.36 instead of ln T - 9.36, and the low-T k16 branch used the T_eV coefficient with T in K (1.75e7x too fast) |
| 3-body formation (Forrey 2013) | 2075 | `three_body.py` | same coefficient | |
| collisional dissociation by H, H2, He, e-, H+ (v=0/LTE interpolation) | 2033-2039 | `collisional_dissociation.py` | same fits | GIZMO's n_crit weights use mass fractions. H2-He+ (and D/D+ x 1e-10) are not in jaco (minor). |
| LW photodissociation 3.3e-11 G0 x GD14 self-shielding | 2018-2030, 2076 | `photochemistry.py` | same form | GIZMO uses the unattenuated ISRF (+UVB) here, but 1.7 ISRF exp(-500 Z Sigma) for photoelectric and C+; jaco has a single G_0 |
| CR dissociation | 2079 | `cosmic_ray_dissociation.py` | same | zeta per H2 |

### 3.7 EOS / thermodynamics

| term | GIZMO | jaco | status | difference |
|---|---|---|---|---|
| H2 internal energy (3:1 frozen ortho:para, rotation + vibration, Boley 2007) | `eos/hydrogen_molecule.cc` 14-123 | `eos/H2_partition_function.py` (fit) | same to <0.3% (E/kT), <1% (c_v) | |
| u(T) of the mixture | 551-605 (`EOS_SUBSTELLAR_ISM`) | `eos/eos.py` | equivalent | GIZMO counts metals as one pseudo-species of mass 16+12 f_mol; jaco counts each metal species. Negligible. |
| gamma / mean molecular weight | `set_eos_pressure` | GIZMO-side (called in `jaco_to_gizmo`) | n/a | |

## 4. New parameter symbols

| symbol | C field | meaning | what GIZMO's `gizmo_to_jaco` must pass |
|---|---|---|---|
| `f_neb` | `pr->f_neb` | switch for nebular cooling | `1.0` under `RT_CHEM_PHOTOION` (with `METALS`), `0.0` otherwise. `Params pr = {}` zero-initializes, so leaving it unset already gives GIZMO's non-RT behaviour. |

No solve variable or existing parameter was added, renamed or removed. Generated `Params` now has 25 fields; it is
alphabetical, so `f_neb` sits between `f_d` and `grad_v`. That only matters to code indexing `pr->data[k]` by number.

## 5. GIZMO-side interface issues in `cooling/jaco.cc` (not editable from jaco)

1. **H2 initial value off by 2x** (jaco.cc 347): `sv->x_H_2 = fmol` and `pr->x_H_2_initial = sv->x_H_2`, but
   `x_H_2` is n_H2/n_H and `jaco_to_gizmo` stores `MolecularMassFraction = 2 x_H_2`. The backward-Euler reference
   therefore doubles every step. The PRIMORDIAL branch correctly uses `0.5 * fmol`.
2. **Dust temperature** is hard-coded to 10 K (jaco.cc 364). The legacy path computes `rt_eqm_dust_temp` from CMB +
   ISRF absorption.
3. **Radiation field in non-RT builds:** `G_0 = 1` and `ISRF = 1` (jaco.cc 368-369). Legacy GIZMO uses
   G0 = 1.7 x ISRF x exp(-500 Z Sigma) (`get_FUV_G0`, 2249-2268) for photoelectric heating and the C+ fraction,
   the unattenuated ISRF for H2 photodissociation, and sqrt(InterstellarRadiationFieldStrength) for CRs. Suggested:
   under `RT_ISRF_BACKGROUND`, pass `ISRF = All.InterstellarRadiationFieldStrength` and `G_0 = get_FUV_G0(i, shieldfac, 0)`.
4. **`z` for non-comoving runs** is 0, while GIZMO's CMB temperature uses `Redshift_RT_Background`. GIZMO's own
   Compton term ignores that parameter, so the two disagree inside GIZMO as well.
5. **N_H convention:** passed as Sigma/m_p (nucleons, not H nuclei). The model now relies on that convention in the
   CR attenuation; it was already what the H2 self-shielding assumed.
6. **Codegen entry point:** `python -m jaco.codegen.gizmo.gizmo` (GIZMO Makefile) is not in any committed jaco branch
   (removed in 546e6aa), and neither are the KWH and PRIMORDIAL models that jaco.cc references.

## 6. Reported, not implemented

- **UVB photoionization, photoheating and Rahmati shielding:** structural (new parameters for Gamma and epsilon).
- **Electrons from heavy ions, alkali, molecular ions and O+; Gong17 C+ balance:** structural (new fixed species).
- **CR heating model:** a model choice. jaco ties heat to atomic-H ionization only, so molecular gas gets none. Fix
  options: give the CR H2 channels a heat per event (~10-13 eV, Glassgold+12), or adopt GIZMO's thermal form.
- **Optically-thick cap and dust temperature:** structural (a cap on the total rate with an opacity model; an
  implicit dust energy balance).
- **[CI] 609 um; the C+ x_H0 factor and H2 collider:** ambiguous (GIZMO's weighting is unphysical) and small.
- **HD abundance and the CO model:** prescription choices that differ deliberately.
- **Metal-table per-element normalization:** convention. jaco's is the consistent one, and it belongs to the
  table-conversion script.
- **Free-free 1.42e-27 vs 1.43e-27:** GIZMO is internally inconsistent; 0.7%.
- **Compton UVB/ISRF/synchrotron:** needs radiation and B-field parameters; negligible.
- **GIZMO's half-speed H2 rates and mass-fraction collider weights:** these are GIZMO-side; jaco is correct.
- **H2-He+ dissociation:** minor.

## 7. Pre-existing jaco issues found

- **NaNs in the generated Jacobian/RHS** (present at f6d0c76; **FIXED 41c2e72, 9413220**). Evaluating the
  compiled generated C on a grid (T = 3 K to 1e9 K, n_H = 1e-4 to 1e8, x_H+ from 1e-6 to 1) gave non-finite
  entries in 3406 of 5382 states before the fixes and in 0 after:
  - H2 collisional dissociation `Max(k, 0)**f0 * Max(k_LTE, 0)**(1-f0)`: d/dT carried log(0) wherever a rate
    underflowed (He collider below ~100 K) or the Savin+04 H+ fit went negative (below 76 K). Above ~5e6 K the He
    critical density underflowed to 0, and 1/n_cr did the same. Rates are now assembled in log space, the H+ fit is
    frozen below 100 K, and the 1/n_cr,He exponent is capped at 500.
  - The GD14 self-shielding f_H2 = 2n_H2/(2n_H2 + n_H) was 0/0 at x_H+ = 1. It is now written in terms of the H2
    column, 2 x_H2 N_H / Sigma_0: the neutral fraction cancels algebraically, as it does in GIZMO.
- **Python vs C table interpolation** disagree (see 3.3).
- **`SolarAbundances.get_abundance_per_H`** uses f/(1-f)/A. This affects the generic `processes/line_cooling.py` C+
  copy and the default `y` prescription in `EquationSystem.solve` (0.0925 vs 0.0951). The starforge model now uses
  `starforge.symbols.x_solar`, which follows GIZMO's packing convention.
- **Reactions that do nothing:** `grain_assisted_recombination("C+")`, `cosmic_ray_ionization("C")` and
  `cosmic_ray_photoionization("C")` act only on species that the reductions fix or eliminate, so they have no effect.
- **C_2 clumping on cooling:** jaco applies the C_2 factor to every 2-body process including cooling; GIZMO applies
  clumping only in its H2 network.
- **Untracked data files:** `tests/chianti_H_abundances.npy` was never committed, so `tests/test_CIE.py` fails;
  `spcool_tables.hdf5` is untracked.

## 8. Tests

`PYTHONPATH=src python3 -m pytest tests src -q`, run with `spcool_tables.hdf5` symlinked from `jaco_gizmo`:

- **before (f6d0c76):** 50 passed, 1 failed (`tests/test_CIE.py`, missing data file).
- **after:** 340 passed, 2 skipped (H- fits below GIZMO's exp(-90) floor), 1 failed (the same `test_CIE`).
  With `-x` the run stops at `test_CIE` both before and after.
- **New GIZMO-parity tests** in `src/jaco/models/starforge/tests/`: `test_nebular_cooling`, `test_CO_cooling`,
  `test_Hminus_detachment`, `test_recombination_cooling`, `test_Cplus_cooling`, `test_metal_line_cooling` (skips
  without the hdf5), `test_gizmo_lowtemp_factors`, `test_photoelectric_heating`, `test_cosmic_ray_rate`, `test_compton`.
  `test_nan_regularization` evaluates its expressions with the generated C's floating-point semantics
  (`tests/c_semantics.py`).
- **Full-model C codegen** (`generate_code(..., minimal=False)`) succeeds before and after: 29 s and 24 parameters
  before, 34 s and 25 after. The compiled code from `python -m jaco.codegen.gizmo.gizmo starforge` (on
  `gizmo_integration` plus these commits) has finite RHS and Jacobian at every grid state (section 7). At 7 regular
  states it agrees with the pre-fix build to 7e-12 in the RHS and 1.2e-10 in the Jacobian; the largest difference is
  in an entry 1e-9 of its row's scale.
