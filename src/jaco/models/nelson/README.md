# Nelson & Langer (1999) model

This Jaco model implements the carbon-oxygen chemistry and thermal balance of low-mass molecular clouds and Bok
globules from Nelson & Langer (1999, NL99, ApJ 524, 923). Its thermal terms draw on NL99, Nelson & Langer
(1997, NL97, ApJ 482, 796), and Hollenbach, Takahashi & Tielens (1991, HTT91, ApJ 377, 192). The CO shielding
table comes from van Dishoeck & Black (1988, vDB88, ApJ 334, 771).

```python
from jaco.models.nelson import make_model
model = make_model()          # solves u, T and eight species; use with Model.solve or generated code
```

When using `Model.solve`, set `JAX_ENABLE_X64=1` before starting Python. In JAX's default 32-bit mode, the coupled
solver can return the initial guess without advancing the temperature or abundances.

To reproduce the paper's elemental composition, pass the values from `jaco.models.nelson.PAPER_PARAMETERS` to the solver.

## Files

| File | What it holds |
|---|---|
| `chemistry.py` | The 23 reactions of NL99 Table 2 with their rate coefficients. |
| `thermal.py` | Nine heating and cooling terms, the HTT91 dust temperature and grain-charging factor. |
| `shielding.py` | vDB88 Table 5 CO shielding, Theta(N(CO), N(H2)), as an expression of two column densities. |
| `model.py` | `make_model()`, the species, solved variables, parameters and `PAPER_PARAMETERS`. |
| `tests/test_nelson.py` | 116 tests (with parametrization), listed below. |

## The model

**Species.** NL99 evolves 12 species: He+, H3+, OHx, CHx, CO, C, C+, HCO+, O, M, M+ and e-. H2 and He are the
background (the cloud is taken to be fully molecular), which makes the 14 species of an explicit network. jaco reads a species' element content from its name, which the composites do not have, so:

| Paper | jaco name | Note |
|---|---|---|
| CHx | `CH` | stands for CH + CH2 |
| OHx | `OH` | stands for OH + H2O + O2 |
| M, M+ | `Mg`, `Mg+` | stands for Mg, Fe, Ca, Na |
| H2 | `H_2` | fixed at 1/2 per H nucleus |

With these names jaco's reductions give the paper's conservation conditions automatically: charge neutrality eliminates
e-, and C, O, He and Mg are eliminated through their totals. The H and H+ products the paper does not follow are left
out of the reactions (so those reactions do not conserve charge or H, which does not matter: e- follows from charge
neutrality and H2 is fixed).

The solver returns `u`, `T`, H3+, He+, CH, OH, CO, C+, HCO+, and Mg+. It advances the eight species and the gas energy
with backward Euler and obtains `u` and `T` consistently through the EOS. `PAPER_PARAMETERS` supplies the elemental
totals per H nucleus (y = 0.14, x_C,tot = 5e-5, x_O,tot = 1e-4, x_Mg,tot = 5e-8). These are half the paper's values
per H2 molecule.

### Host inputs

The host supplies `n_Htot`, the starting thermodynamic state, and any overrides for the defaults below.

| Parameter | Default | Meaning |
|---|---|---|
| `G_0` | 1 | FUV field in Habing units |
| `A_V` | 0 | visual extinction: photoreactions and dust temperature |
| `tau_uv` | 0 | UV optical depth of the photoelectric and surface terms (see below) |
| `N_CO`, `N_H2` | 1 | column densities (cm^-2) for CO shielding; 1 = unshielded |
| `dv_dr` | 3.24e-14 | velocity gradient (s^-1) of the C I escape probabilities, NL99's 1 km/s/pc |
| `pdv_work` | 0 | compressional heating (erg cm^-3 s^-1) |

## Chemistry (`chemistry.py`)

Rate coefficients are NL99 Table 2 as printed, in T (K), not T/300: two-body in cm^3 s^-1 and the cosmic-ray and photo
rates in s^-1. Photoreactions are `alpha * G_0 * exp(-gamma * A_V)`, and the CO one also carries the vDB88 shielding.
The cosmic-ray ionization of H2 uses the paper's 1.2e-17 s^-1, while the heating (below) uses NL97's 2.0e-17 s^-1, as
the papers do.

## Thermal terms (`thermal.py`)

| Term | Source | Note |
|---|---|---|
| Molecular line cooling | NL97 eq. 7 | fit to Goldsmith & Langer 1978, in n(H2) |
| C+ cooling | NL97 eq. 8 | n_crit = 3e3 cm^-3 |
| O I cooling | NL97 eq. 9 | |
| C I cooling | NL99 eq. 7, Genzel 1991 | collision rates gamma = A / n_crit (Genzel Table 11), Sobolev escape probabilities for `dv_dr`, closed-form two-pass iteration |
| Surface cooling | NL97 eq. 10 | artificial, keeps the low-density edge below ~100 K |
| Cosmic-ray heating | NL97 eq. 11 | 2.0e-17 s^-1 x 20 eV per H2 (6.4e-28 n) |
| H2 formation heating | NL97 eq. 12 | |
| Gas-dust exchange | NL97 eq. 13 | dust temperature from HTT91 eq. 5, floored at 10 K |
| Photoelectric heating | NL99 eq. 8, HTT91 B7-B9 | 4.86e-26 G_0 n_H f exp(-tau_uv), f = (1-x)^2/x, per H nucleus |

The EOS is jaco's own: an ideal gas from the declared species, which gives gamma = 5/3 at these temperatures as NL99
assumes.

## Model assumptions and approximations

1. **Photoelectric heating per n_H.** The model uses HTT91's coefficient of 4.86e-26 per H nucleus. NL97 gives
   6.0e-26 per H2 molecule. In a fully molecular cloud, n_H is twice n(H2).
2. **`tau_uv` is host-supplied.** NL97 averages exp(-tau) over the six coordinate directions (eq. 17), so a host should
   pass the effective tau_uv = -ln(mean exp(-tau_i)), with tau_i = 1.25e-21 N_H per direction.
3. **C I cooling** follows NL99 eq. 7 with Genzel (1991). The collision rate coefficients with H2 are the constants
   A / n_crit from Genzel's Table 11: for the 609 um line n_crit = 1e4 cm^-3 is given for H2; for the 370 um line the
   table lists only 1.2e3 (for H), which is used as printed (NL99 do not say what they did). Genzel gives no formula
   for the optical depth: the model uses the standard Sobolev expression, tau = A lambda^3 / (8 pi dv/dr) n(C I)
   (f_l g_u/g_l - f_u), local for a velocity gradient `dv_dr` (NL99: 1 km/s/pc), and no 2-0 collisions.

   The escape probabilities depend on the level populations. Starting from the ground-state estimate, each of
   `C_I_PASSES` iterations solves the populations in closed form and then updates the escape probabilities. For
   T = 8-100 K and n(C I)/n(H2) <= 3e-5, two passes give a worst-case relative error of 3% and a median of 5e-7
   against the converged solution. One pass gives a worst-case error of 16% and a median of 5e-4, with a Jacobian
   about seven times smaller.
4. **Grain yield** Y = 0.1 in HTT91 B6 (their B9 corresponds to 0.103).
5. **Photoelectric coefficient and grain charging.** The model uses the HTT91 coefficient of 4.86e-26 and the
   HTT91 fit to the grain-charging cubic in eq. B8. NL99 prints the coefficient as 4.86e26.
6. **CO shielding** uses bilinear interpolation in ln Theta with clamping at the table edges. Both log-column axes
   have a wide first interval followed by unit spacing. Jaco's `table_interp_2d` requires uniform spacing and
   lacks the JAX implementation needed by `Model.solve`. The SymPy piecewise expression preserves the published
   grid and works in the Python solver and in generated code.
7. **Bergin et al. 1995** is the paper NL99 cites for its rate coefficients and shielding, but it contains neither the
   rates nor the shielding table (it uses vDB88's rates); it is not used.

## Tests (`tests/test_nelson.py`)

- Declarations: `jaco.check` is clean; the 12 species plus H2 and He; 23 reactions.
- Every rate coefficient against an independent transcription of Table 2, at three temperatures.
- Each thermal term against its paper equation (comparisons are purely relative, `abs=0`: pytest's default absolute
  tolerance of 1e-12 would make them pass for the 1e-17 to 1e-30 values here); C I cooling in its low-density (linear) and LTE (saturated) limits,
  against Genzel's two-level formula (eq. 24), against an independent converged solution of the rate equations and
  escape probabilities, and against the GL78 C I peak values (read from its plot).
- Dust temperature (HTT91 eq. 5 and the 10 K floor), grain-charging factor and its limits.
- vDB88 shielding: every Table 5 node, printed values, interpolation, clamping, monotonicity, and the CO rate.
- Through `Model.solve`: He, C, O, Mg and charge conservation (solved and eliminated species), CO increasing with A_V
  (the reactions alone at fixed T), and a thermal-balance root between 3 and 300 K (the thermal terms alone, reached
  from 3 K and from 300 K). These tests enable JAX 64-bit precision.
- One coupled timestep with the full model at n_H = 1e3 cm^-3. The test checks each evolved species balance, the
  energy balance, and the EOS at the returned state.

Run with `python -m pytest src/jaco/models/nelson` from the `jaco` directory.
