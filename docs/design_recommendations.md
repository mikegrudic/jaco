# jaco model abstraction: design recommendations

Reviewed 2026-10-02 at jaco `gizmo_integration` ee50bff and GIZMO `gizmo_jaco_dev` 7921120a.
Two independent reviews (Fable 5.1 and a blind Opus 5.5 pass) converged on the same list. Every
defect below was reproduced with a script against the code; file references are relative to
`src/jaco/` unless noted.

## Verdict

Keep the core. A `Process` is the atom, `ChemicalReaction("H + H -> H_2", k, ...)` is the DSL,
`+` composes, and species strings drive the conservation reductions. Most of the starforge model
is one-liners at that level, which is what frictionless looks like.

The friction is in what wraps the core: composition has side effects, model configuration is
spread across module globals and ad-hoc dict attributes and is computed too early, parameters
are inferred from leftover symbols, and the host contract is maintained by hand. Two changes
(a real model object, declared parameters and species) alter what a downstream developer's day
looks like. The rest is debt.

## Do not change

- `+` composition, the string reaction syntax, the symbolic Jacobian, the reduction order.
- The one-file-per-physics layout and the intent of `models/model.md`.
- The table digest check between compiled code and `jaco_tables.hdf5`.
- No YAML/JSON network format, no plugin registry, no deeper class hierarchy, no investment in
  the Julia/Fortran/CUDA printers until a second consumer exists.

## Verified defects (fix regardless of the redesign)

| # | Defect | Evidence | Fix |
|---|--------|----------|-----|
| D1 | `+` mutates its operands | `EquationSystem.__add__` reads `other[k]`; `__getitem__` creates missing keys. A `ThermalProcess` with keys `['heat']` has `['H','H+','e-','heat']` after being added to a recombination. Module-level objects (`compton_cooling`, `H2_chemistry`) accumulate every species they were ever summed with; `chemical_species` and the EOS follow. | Non-creating lookup in `__add__`. Regression test: build A then B equals fresh B. |
| D2 | Rate setters accumulate | `Recombination.rate_coefficient = k` calls `update_network()`, which does `-=` on the existing rows. Reassigning the same rate doubles the term. `Ionization.rate` same. `ChemicalReaction` ignores a post-construction rate change while `.rate` reports the new value. | Immutable processes; rates only in constructors. |
| D3 | Clumping is applied inconsistently | `CollisionalIonization` (via `Ionization`) has no `C_2`; `GasPhaseRecombination` (via `NBodyProcess`) multiplies by `C_2`. With `C_2 = 1 + Mach^2/4` this is a physics asymmetry (26x at Mach 10). | Explicit per-reaction `clumping=` keyed on material reactants only. **Decided 2026-10-02:** every 2-body reaction between material species carries `C_n` by default, ionization included; a separate `starforge_legacy` model, which replicates the terms of GIZMO's legacy cooling solver, overrides ionization to no clumping. |
| D4 | Stale closures survive composition | `starforge.py:110-119` derives the H- steady-state expression from the summed network; `make_model() + new_Hminus_reaction` keeps the stale expression and the new rate vanishes from the reduced system. | Closure *rules* evaluated at reduction time (see R1). |
| D5 | Metals enter the EOS by accident | Fe, Mg, ... are in `chemical_species` only because metal-line rates mention `x_Fe`. Removing metal cooling changes mu and `n_Htot`. | Declared species (R2). |
| D6 | Silent 2-variable system | `Process.generate_code()` defaults to `["u","T"]`; the model's `SOLVE_VARS` is a module global it never sees. | Solve variables live on the model (R1). |
| D7 | Table name collision | `metal_line_cooling` names tables `_z{int(round(z))}`; z=0 and z=0.3 collide, last registered wins. | Name by slice index; registry removal (R4). |
| D8 | Import-order dependent model lookup | `models/__init__.__getattr__("starforge")` shadows the subpackage of the same name. | Delete the `__getattr__`; CLI already uses `import_module`. |
| D9 | Host drops two couplings | `gizmo_jaco_dev/cooling/cooling.cc:180` returns after `call_jaco`, before `CR_cooling_and_losses` (235) and the cooling-radiation return to RT bins (243-290). | Add `COSMIC_RAY_FLUID` and `RADTRANSFER` to the `#error` list at `jaco.cc:26` now. |
| D10 | Misc | H2 photodissociation process is named "Photodissociation of H-"; `Recombination` stores `colliding_species` as a set (self-collision collapses); `sanitize_symbols` is quadratic (`subs` in a loop, use one `xreplace`); undeclared `h5py`/`scipy` deps, declared `sphinx_rtd_theme`/`matplotlib` runtime deps; bibliography warning fires per construction. | One-liners. |

## Recommendations

### R1. One declarative `Model` object

**Problem.** `make_model()` returns a `Process` and pokes `fixed_species`, `equilibrium_overrides`
and `derived_params` onto the `EquationSystem` dict as attributes. `SOLVE_VARS`/`TIME_DEPENDENT`
are module globals. The `Model` dataclass in `models/model.py` is unused. `fixed_species` and
`equilibrium_overrides` carry identical content for H_2+ and HD; the overrides path is dead because
`eliminate_steady_state` is never enabled, while starforge hand-rolls the same linear solve for H-.

**Fix.** `Model` owns: a *keyed* process collection, solve variables, time-dependent set,
`steady_state=["H-", ...]` (worked out at reduction time with the linear elimination already in
`do_conservation_reductions`), fixed and derived species, and parameter declarations (R2).
Operations: `+` (error on duplicate id), `without(id)` (error if absent), `replace(id, new)`.
Merging two models with conflicting declarations raises. Codegen reports every equation it discards.

### R2. Declared parameters and species

**Problem.** `Params` is "whatever free symbols remain" (`_generate_full`). A typo (`G0` for
`G_0`) or a 4-body reaction's `C_4` becomes a new struct field the host zero-initialises. `solve()`
hardcodes defaults for `y`, `Z`, `C_2` in a `prescriptions` dict. Species are identified by element
parsing: `species_counts("NUV")` is N+U+V with 122 electrons, `"He1"`/`"He2"` are helium compounds,
and `"photon_*"` is massless only because lowercase p is not an element.

**Fix.** `Parameter(name, units, default, doc)` and a species declaration with a *kind*
(material, radiation, energy reservoir) checked before element parsing. Kind drives EOS membership,
conservation participation, clumping order and the n-to-x conversion. Codegen refuses any
undeclared free symbol. The generated header emits `JACO_HAS_PARAM_<name>` / `JACO_HAS_VAR_<name>`
so one hand-written `jaco_pack_params` fills each quantity once under `#ifdef`, replacing the
per-model `#ifdef JACO_MODEL_*` ladder (which currently carries PRIMORDIAL and KWH branches for
models that do not exist). Debug builds NaN-poison `Params` before packing and check afterwards.

### R3. Immutable, flat processes

**Problem.** Three construction idioms for one concept (D2, D3). `Ionization` and
`Recombination` resync the network from property setters with name-mangled state.

**Fix.** Everything is a `Reaction` (stoichiometry, rate coefficient, heat per reaction, explicit
clumping, name, bibliography) or a `ThermalTerm`, immutable after construction. Ionization and
recombination become thin factories. Factories take prescription keywords
(`f_selfshield_H2(prescription=)` is currently not reachable from `photodissociation()`;
`three_body` picks rate 5 at import). Expose lists, not pre-summed bundles (`H2_chemistry`,
`LineCoolingSimple(collider=None)`), so one reaction can be removed. Unique names.

Clumping default: all 2-body reactions between material species carry `C_n`. A model may override
per reaction or per reaction class (`starforge_legacy`: no clumping on collisional ionization, matching
GIZMO). The override is a model-level rule (R1), not an edit to the process.

### R4. Remove the global table registry

**Problem.** `_TABLE_REGISTRY` in `interpolation.py` is process-wide state filled as a side
effect; `generate_funcjac_code` regex-searches the emitted C for `&name` to recover which tables
a model uses.

**Fix.** Tables are already sympy atoms (`TableInterp2D`) carrying their name; collect by walking
`expr.atoms(TableInterp2D)` on the reduced expressions and carry the data on the atom or in a
per-model dict. Delete the registry and the regex.

### R5. Outputs channel instead of pruning

**Problem.** `prune_decoupled` deletes `photon_rec`, `photon_assoc,H` and `dust heat`, which
are exactly the per-band sources GIZMO's `Rad_Je` wants. A `photon_rec` product does not appear
anywhere in the generated C.

**Fix.** Emit the decoupled rows as a second, Jacobian-free function evaluated at the converged
state. Under backward Euler, dt times the RHS at the final state is the integral over the step, so
the host applies it directly. Covers band emission, gas-dust energy into the IR band, cosmic-ray
energy loss, optically thin absorption. Most RT/CR coupling is outputs, not new unknowns.

### R6. Solver reads metadata from the generated header

**Problem.** `jaco_solver.cc:94` has `is_time_dependent(k) { return k == IDX_x_H_2; }`; abundance
floors/ceilings assume x at most 1 (a photon abundance per H above 1 is physical); the eliminated-
abundance bound `1 - x_H+ - 2 x_H2 >= 0` is hand-coded. Promoting C+ to a solve variable leaves
`x_C` unbounded with no complaint.

**Fix.** Codegen emits per-variable time dependence, floor, ceiling, scale, and the eliminated-
abundance bounds (with gradients) from the substitution list it already builds. Prerequisite for
any non-starforge solve variable.

### R7. `jaco.check(process)` harness

`models/starforge/tests/c_semantics.py` + `test_nan_regularization.py` is the right idea
(C float semantics: clamped `exp`, `fmax`/`fmin`, underflow) but bespoke to two expressions.
Generalise: evaluate every rate and its partials over a standard (T, n, x) grid, compare the
symbolic Jacobian to finite differences, flag undeclared symbols. Contributors find `0 * inf`
before GIZMO does.

### R8. Split `EquationSystem` (last, optional)

Container, reduction pipeline, JAX solver driver and C/CUDA/Python/Julia/EOS codegen share one
1100-line class; species detection is `str(s)[:2] == "n_"`. Move codegen to `codegen/` and the
solver to `numerics/` as free functions. Mechanical. Do it only if R1-R7 have not already made the
reduction pipeline readable; a second blind review warned this risks premature abstraction.

## Photons and cosmic rays as species

Treating photons and CRs as species does not change the list above; it adds R5, R6 and the
species-kind half of R2, and it bounds the framing:

- **Holds for ionizing bands.** `H + photon_EUV -> H+ + e-` already reduces to a backward-Euler
  unknown with its `_initial`, passes codegen, stays out of the EOS. Only the `C_2` factor is wrong.
- **Fails for cosmic rays.** `CR + H -> CR + H+ + e-` gives a CR row of exactly zero (catalyst);
  the energy loss per event is inexpressible. GIZMO integrates CR losses exactly as `exp(-Lambda dt)`,
  so Newton gains nothing. Keep CR energy a parameter plus a loss output (R5).
- **Fails for dust-absorbed bands** (opacity x energy, with a separate IR radiation temperature)
  and **Lyman-Werner** (6-7 photons per dissociation; `species_and_coeffs` accepts integers only).
- **GIZMO's physics is not simple stoichiometry.** Chemistry ionizes with c, RT removes photons
  with reduced c-tilde, so a photon row carries a runtime c-tilde/c. Heating per photoionization is
  the cross-section-weighted `rt_ion_G_HI`; `rt_chem.cc:36-42` says explicitly it is not
  nu_eff - 13.6 eV, so automatic reaction enthalpy from band-averaged photon energy is wrong.
  Band sigma, G, h nu_eff are computed at runtime from `star_Teff` and must be `Params`.
- **Bands.** The model declares the bands it needs; GIZMO maps them with static asserts. Not
  `make_model(config)` parsing GIZMO's band list.
- **Gate.** Iliev test 1 (`RT_ILIEV_TEST1` hook exists) with c-tilde sigma n_H dt >> 1,
  I-front radius vs analytic, outputs-only vs coupled. If outputs-only passes with absorbed capped
  at available, never make photons Newton unknowns. Only if it fails: a transported species role
  (backward-Euler unknown + host transport source on the `pdv_work` pattern + runtime c-tilde/c),
  ionizing bands only.

## Order of work

1. Golden test: hash of the generated starforge sources, plus regression tests for D1, D2, D4.
2. D1 (`__add__`), D8, D9, D10. Generated C stays byte-identical.
3. R3 immutable processes and factories; apply the D3 decision (clumping default on, `starforge_legacy` override).
4. R1 `Model` with reduction-time closure rules.
5. R2 declared parameters and species kinds; `JACO_HAS_*` macros; codegen refuses undeclared symbols.
6. R4 table registry removal.
7. R5 outputs channel and R6 solver metadata. Needs a GIZMO test round.
8. R7 check harness.
9. R8 only if still warranted.

Steps 1-2 are small and keep output byte-identical. Steps 3-5 change the contributor experience.
Each step lands independently against `gmc_cooling[jaco]` and `isodisk_thermalfb[jaco]`.
