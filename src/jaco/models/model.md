# Implementing a model in `jaco`

A model is a set of physical processes, the assumptions that close them into a solvable network, and the declarations
that tell the host code what the network needs and produces. The code generator turns it into the C sources GIZMO
compiles under `JACO=<model>`; the same object is what `jaco.check` evaluates and what the Python/JAX solver solves,
within its limits (section 9).

The guiding principle is the opposite of the single-file, 8000-line Fortran 77 cooling routine that solves all of
reality in one function. Each physical process is one object with a name and a bibliography; a model is a list of them
plus declarations; everything else (the conservation sums, the steady-state eliminations, the Jacobian, the EOS, the
host interface) is derived.

This document walks through a minimal model, then lists every piece of the contract, then says what the GIZMO host
actually supplies, which is where new models go wrong. File references are relative to `src/jaco/` unless they name
`cooling/jaco.cc`, which is GIZMO's host layer (`gizmo_jaco_dev`, commit e34e4f3b when this was written).

## 1. A minimal model

A model lives in a subpackage `jaco.models.<name>` whose `__init__` exposes `make_model()` returning a `Model`.
The smallest real one is `models/wind_comparison` (an atomic-H cooling curve, a heating floor, PdV work; solves `u`
and `T` only). Here is one with a reaction as well:

```python
# src/jaco/models/toy/__init__.py
import sympy as sp
from jaco.model import Model
from jaco.processes import Reaction, ThermalTerm
from jaco.declarations import Species
from jaco.symbols import T, n_, n_Htot, piecewise_powerlaw
from jaco.math import logistic

EV = 1.602e-12
cooling = ThermalTerm(-piecewise_powerlaw([10, 1e4, 1e5, 1e8, 1e10], [1e-30, 1e-23, 1e-21, 1e-22, 1e-21], T,
                                          extrapolate=True, name="lambda_toy") * n_Htot ** 2,
                      name="toy cooling", bibliography=["made up"])
heating = ThermalTerm(2e-26 * n_("H") * logistic(-(T - 1.5e4) / 1e3), name="toy heating",
                      bibliography=["made up"])  # a heating floor that switches off above ~1.5e4 K
pdv = ThermalTerm(sp.Symbol("pdv_work"), name="PdV work")  # required: the host's compressional heating
ionization = Reaction("H + e- -> H+ + e- + e-", 5.85e-11 * sp.sqrt(T) * sp.exp(-157809.1 / T),
                      heat_per_reaction=-13.6 * EV, clumping=1,
                      name="H collisional ionization", bibliography=["1992ApJS...78..341C"])
recombination = Reaction("H+ + e- -> H", 2.59e-13 * (T / 1e4) ** -0.7, clumping=1,
                         name="H+ recombination", bibliography=["1992ApJS...78..341C"])
cosmic_rays = Reaction("H -> H+ + e-", 1e-16, heat_per_reaction=10 * EV,
                      name="cosmic-ray ionization of H", bibliography=["made up"])  # an electron source besides H+ itself


def make_model():
    return Model([cooling, heating, pdv, ionization, recombination, cosmic_rays],
                 solve_vars=["u", "T", "H+"], time_dependent=["T"],
                 species=[Species("H"), Species("H+"), Species("e-")])
```

Three things in it are not obvious and are explained below. The PdV term is mandatory (codegen refuses a model whose
heat has no `pdv_work`). `clumping=1` is explicit because the default clumping factor `C_2` is a parameter GIZMO does
not fill (section 5); the heating rate is a literal rather than a `Parameter` for the same reason. And the cosmic-ray
line is there because a network whose only electron source is its own ion has a second steady state at `x_e = 0`
(ionization and recombination both vanish with `n_e`), and the solver lands on it from hot starts; every real model
carries cosmic-ray or photo-ionization or `fixed_electrons`. Without that line the toy fails 115 of the driver's 900
cases.

What the author did not write: `e-` is eliminated by charge conservation and `H` by H conservation
(`x_H = 1 - x_H+`); the EOS `u(T, x)` from the material species; the backward-Euler row for `T` with its implied
parameters `Δt` and `u_initial`; the Jacobian; the host struct layouts and the `JACO_HAS_*` macros. `H+` is solved in
steady state because it is not in `time_dependent`.

Check it, generate it, run it through the solver driver:

```sh
export PYTHONPATH=/path/to/jaco/src    # or pip install -e the jaco tree
python -c "import jaco; from jaco.models.toy import make_model; r = jaco.check(make_model()); print(r.summary()); assert r.ok"
python -m jaco.codegen.gizmo.gizmo toy --language c --ext .cc --output-dir /tmp/toy   # what GIZMO's Makefile runs
make -C test/jaco_solver MODEL=toy test      # in the GIZMO tree; models not named starforge* or sfx_* get the generic driver
```

What to expect, measured on gizmo_jaco_dev e34e4f3b: `jaco.check` passes (54 expressions at 7800 points); the codegen
writes seven files; the generic driver reports `900 cases: ... FAILED 0; 0 check failures`. In that summary `FAILED`
counts cases where every solver tier failed and `check failures` counts every failed case, including those. A sound
model passes all 900 (`wind_comparison` does too). One known false positive: the generic driver checks an answer the
subcycling tier produced against the full step's system rather than the last substep's, which the starforge driver
does through `JacoSolveInfo.substep_params`; tier-3 answers can therefore fail its check while being right.

GIZMO builds the model with `JACO=toy` in `Config.sh` (requires `COOLING`); the Makefile runs the codegen through
`JACO_PYTHON` (default `python3`). The codegen stamp `cooling/.jaco_codegen_stamp` depends on the Config file only, so
after editing the model delete the stamp (or `make clean`) before rebuilding.

## 2. Processes

Processes are immutable. Their rows are rates per unit volume and time: a species row in cm^-3 s^-1, the gas `heat`
row and energy-reservoir rows in erg cm^-3 s^-1. Abundances are per H nucleus (`x_X = n_X / n_Htot`); the symbols are
`jaco.symbols.T`, `n_("X")`, `x_("X")`, `n_Htot`, and any declared parameter as `sp.Symbol("name")`. Units are CGS.
Python symbol names are sanitized into C identifiers (`H+` becomes `Hplus`, `Δt` becomes `Delta_t`, `x_C,tot` becomes
`x_C_tot`); write the Python form in expressions and declarations.

- `Reaction(equation, rate_coefficient, heat_per_reaction=0, *, rate=None, heat_rate_coefficient=None, clumping=None,
  row_factors=None, name="", bibliography=())` (`processes/chemical_reaction.py`). The equation string
  (`"3H -> H_2 + H"`, `"H+ + e- -> H"`; species are split on `" + "`, so keep the spaces) drives the stoichiometry; species on both sides keep
  their rows; the parser does not check that the equation balances. Events per unit volume and time are
  `clumping * k * prod(n_reactant)`, or `clumping * rate` for a non-mass-action rate. Each event removes the reactants,
  adds the products and gives the gas `heat_per_reaction` (erg; negative for cooling). Clumping defaults to `C_k` for a
  mass-action rate with `k >= 2` reactants (reactants a model declares as radiation do not count) and to 1 for a
  `rate=`; pass `clumping=` to say otherwise. `row_factors` scales individual rows: a radiation band transported at a
  reduced speed of light, or counted in energy rather than photons. `name` defaults to the equation; a missing
  bibliography warns.
- `ThermalTerm(heat, name="", bibliography=(), *, clumping=1, reservoir=None)` (`processes/thermal_process.py`): a
  heating (positive) or cooling (negative) term of the gas; `reservoir="dust heat"` sends what the gas loses to that
  energy row, so the two sum to zero. It needs a name. Do not use `PdV_heating` from the same module for the host's
  PdV work: it is a gravitational-compression estimate with its own undeclared symbols; the host's term is
  `ThermalTerm(sp.Symbol("pdv_work"), name="PdV work")`.
- `Transfer(rate, rows, name="", bibliography=())`: power moved between rows that are not species abundances (gas heat,
  energy reservoirs, bands counted in energy), each row getting `factor * rate`.
- `a + b` composes processes; the sum keeps `subprocesses`, so nothing is lost when a model lists a composite.

Ready-made families with their references live in `processes/` (collisional ionization, recombination, line cooling,
free-free, ...). Rates with tabulated data use `piecewise_linear`, `piecewise_powerlaw` (1-D, inlined as static data)
or `table_interp_2d/3d` (`symbols.py`); 2-D/3-D tables are written to `jaco_tables.hdf5` at codegen and loaded at run
time, with a digest compiled into the code so a stale table file is refused.

Every process is identified by its name, which must be unique in a model, and should carry a `bibliography` (bibcodes
or a note). The point is that `model.processes` is a complete, citable inventory of the physics.

## 3. The `Model`

`Model(processes, *, solve_vars, time_dependent, steady_state, fixed, derived, intermediates, fixed_electrons, rules,
parameters, species, outputs, variables)` (`model.py`; `Rule` is in `jaco.model`, `Species`, `Parameter`, `Variable`
and `Output` in `jaco.declarations`). Composite processes are split into their atoms; a duplicate id raises. The model
is immutable; `model + process` adds one, `model + other_model` merges two (conflicting declarations raise),
`model.without(*ids)` and `model.replace(id, new)` take out or swap. Closures, rules and reductions are worked out when
the network is assembled, so a process added later is treated like the others.

- `solve_vars`: the solve variables in index order, always `("u", "T", <species or variables>...)`; the solver asserts
  `IDX_u == 0`, `IDX_T == 1`.
- `time_dependent`: the solve variables that get a backward-Euler term over the step; the rest are solved in steady
  state. `T` must be time-dependent; its row is `ρ(x) (u(T, x) - u_initial) / Δt` minus the net heating, with
  `ρ = Σ m_s n_s` over the material species. A time-dependent species `v` gets the implied parameter `x_<v>_initial`.
- `steady_state=[...]`: species eliminated by their own rate equation at assembly (`H-` in starforge). The equation must
  be linear in the species' density; the abundance is clamped at 0.
- `fixed={species: expr}`: a prescribed abundance. `derived={parameter: expr}`: an expression of the solve variables
  substituted for a parameter before differentiation (starforge derives its clumping factor `C_2` from `∇v` this way,
  which is how it gets a nonzero `C_2` in GIZMO; see section 5).
- `intermediates=[(symbol, expr), ...]`: named quantities evaluated once per call before the equations; use them for
  anything shared by many rows (an optical depth, a cap factor) so the generated code and its Jacobian do not expand it
  into every term.
- `fixed_electrons=expr`: free electrons per H nucleus from species outside the network (metal ions in starforge).
- `rules=[Rule(name, apply, exempt=(), only=None)]`: rewrites of the processes at assembly (a different density in the
  rates; clumping only in cold gas). Prefer a rule over editing processes: it applies to processes added later too.
- `variables=[Variable(name, row, kind="temperature", floor, ceiling, scale)]`: solve variables that are neither `u`,
  `T` nor a species abundance, each fixed by the steady state of one row. The dust temperature is one:
  `Variable("Td", "dust heat", ...)` with `"dust heat"` an energy species, whose row is then the equation and not an
  output. A balance variable must be in `solve_vars` and not time-dependent (`Model` raises otherwise); the solver
  accepts only a stable root of its row (decreasing through zero) and brackets it where Newton fails.
- `outputs=[Output(name, expr=0, units="", doc="", heat_of=..., sums=...)]`: quantities evaluated at the converged
  state without derivatives (under backward Euler, `Δt` times the value is the integral over the step). `heat_of`
  sums the heat of named processes (process id -> weight; `(id, row)` takes a row instead), which is how the energy
  ledgers and per-band emissions are booked. `name` must be a C identifier (`Outputs.<name>`, `IDX_OUT_<name>`,
  `JACO_HAS_OUTPUT_<name>`). The rows of radiation and energy species that are not solved for are outputs already. An
  output the host does not book by name is computed and dropped.

## 4. Declarations: the contract

Species must be declared (codegen and `jaco.check` refuse a model without them); every row of the network must be
declared, and code generation refuses any undeclared symbol, so a typo or a stray `C_4` fails at codegen instead of
becoming a struct field the host zero-fills.

- `Species(name, kind="material", doc="", floor=None, ceiling=None, scale=1.0)` (`declarations.py`). Kinds:
  - `material`: has an element composition parsed from its name (`H_2`, `He++`, `CO`, `e-`); enters the EOS, the charge
    and element conservation sums, and the clumping order of its reactions. The EOS is built from exactly these.
  - `trace`: material too dilute to count (HD): in the clumping order, out of the EOS and the sums.
  - `radiation`: photons of a band, per H nucleus like a species, no mass, charge or clumping (`photon_EUV`).
  - `energy`: an energy reservoir whose row is an energy density (`dust heat`).
  `floor` (default 1e-20 per H), `ceiling` (1 for material and trace, none for radiation and energy) and `scale` (the
  abundance at which the species starts to matter; absolute tolerances are fractions of it) bound and scale the solver
  where the species is a solve variable.
- `Parameter(name, units="", default=None, doc="")`: an input, e.g. `Parameter("Gamma_0", "erg s^-1", 1e-26, "heating
  per H atom")` with `sp.Symbol("Gamma_0")` in a rate. The default is used by `jaco.check` and the Python solver only;
  the C side has no defaults (section 5). `jaco.check` reports a parameter without a default and fails; give one, or
  pass `jaco.check(model, params={"name": value})`.
- Parameters the model gets without declaring them:
  - core (`CORE_PARAMETERS`): `n_Htot`, `T` (where not solved), `pdv_work`, `y` (He nuclei per H), `Y`, `Z`, `C_2`,
    `C_3`;
  - implied by the declarations: `x_<species>` for every non-energy species (an input where the species is neither
    solved nor eliminated), the element totals `y` and `x_<El>,tot` for every element of a material species, and for
    the time-dependent variables `Δt`, `u_initial` and `x_<v>_initial`. Only the symbols the reduced system uses end
    up in the generated `Params`.
- The generated header defines `JACO_HAS_VAR_<name>` for every solve variable, `JACO_VAR_TIME_DEPENDENT_<name>` for
  the time-dependent ones, `JACO_HAS_PARAM_<name>` for every parameter and `JACO_HAS_OUTPUT_<name>` for every output,
  so host code fills and reads what the model has under `#ifdef`; plus `JACO_MODEL_<NAME>` and
  `JACO_FAMILY_<FAMILY>` (`GIZMO_FAMILY` in the model module, defaulting to the model name).

## 5. What the GIZMO host supplies

`cooling/jaco.cc` packs `Params` from the cell in `jaco_pack_params` (what the EOS needs: density, abundances, `y`,
the metal totals) and `gizmo_to_jaco` (the rates' inputs). It knows a fixed set of names. **A declared parameter with
any other name compiles and is zero in production**; `jaco.check` cannot see this because it uses the declared default.
Who supplies what:

| parameter | GIZMO | `jaco.check` | generic driver |
|---|---|---|---|
| `n_Htot` | the cell | grid | grid |
| `Δt`, `u_initial` | the cell's step (times the subcycle fraction under `TRANSPORT_SUBCYCLE_COOLING`) and `InternalEnergy` | 1e10 s; `u(T, x)`, so the backward-Euler term vanishes | grid; `u(T0)` |
| `pdv_work` | the cell's, 0 when the step is operator-split | 0 | grid |
| `y` | the cell | solar | 0.0994 |
| `x_H` | the constant 1 | 1 | 1 |
| `Y`, `Z`, `C_2`, `C_3` | **no (zero)** | default (`C_2 = C_3 = 1`) | 0 |
| `x_<v>_initial` of a time-dependent species | `H+`, `He+`, `He++` from the cell's `Ne` and `H_2` from `MolecularMassFraction`, only when all three ions are solved (`JACO_SOLVES_IONS`); `H+` from `HII` under `JACO_TRACKS_HII`; bands from `Rad_E_gamma`; **0 for any other** | start state | seed |
| `Delta_x`, `f_d`, `f_metal`, `f_neb`, `G_0`, `G_LW`, `grad_v`, `grad_v_tf`, `ISRF`, `kappa_thick`, `N_H`, `N_H_cool`, `Td`, `X`, `z`, `Z_d`, `x_C,tot`, `x_O,tot`, `x_Ca`, `x_Fe`, `x_Mg`, `x_N`, `x_Ne`, `x_S`, `x_Si` | yes, each under `#ifdef JACO_HAS_PARAM_` | default | 0 |
| `A_dust`, `rho` | only for a model that solves `Td` without an IR band | default | 0 |
| `rsol`, `eps_HI` | only for a model with bands | default | 0 |
| `T_rad`, `rho`, `Td_initial` | only for a model with an IR band | default | 0 |
| anything else | **zero** | default | 0 |

So a model with two-body chemistry must either derive `C_2` (as starforge does from `∇v`) or pass `clumping=1`; a
time-dependent species the cell does not store restarts from zero every step; and a new physical input needs a line in
the packer before it does anything: `gizmo_to_jaco` for a rate's input, `jaco_pack_params` for anything the EOS reads
(a fixed abundance, `fixed_electrons`), since the EOS calls pack only there. Build once with `JACO_DEBUG_PARAMS` in
`Config.sh` (not in `Template_Config.sh`): every `Params` field starts as a NaN; the solve path aborts naming the
first field the packer left unset, the EOS path aborts without naming it. It cannot catch the time-dependent-species
case, since that field is written.

Models with radiation bands have more unconditional traffic: `jaco_pack_bands` writes `sigma_HI`, `hnu_EUV`, `T_CMB`,
`Z_metals`, `gamma_eos`, `G_LW_bg`, `gamma_12_UVB`, `eps_H0_UVB` without `#ifdef`, so a band model that lacks any of
them does not compile; and the host reads the outputs `photoionization_rate` and, with an IR band, `T_rad_new`.

What the host requires of the variables, enforced by the `#error` guards at the top of `jaco.cc`:
- `u` and `T`, `T` time-dependent (the energy equation; `u_initial` is the cell's energy at the start of the step).
- `H+`, `He+` and `He++` together as solve variables define `JACO_SOLVES_IONS` (and then `y` must be a parameter):
  the ions are seeded from a CIE table and the cell's `Ne`, and `Ne` is written back. Without all three the JACO path
  never writes `Ne`. Under `RT_CHEM_PHOTOION` the model must solve the ions with `H+` time-dependent
  (`JACO_TRACKS_HII`: the cell's `HII` is the start value and is written back).
- `H_2` as a solve variable writes `MolecularMassFraction` back (the write-back reads `x_H+`, so `H_2` without `H+`
  does not compile) but reads it in as the start value only under `JACO_SOLVES_IONS`: with `H_2` time-dependent and
  fewer than three ions, it restarts from zero every step.
- `Td`: required with `RT_INFRARED`; a model that solves `Td` without an IR band is refused under
  `COOL_LOW_TEMPERATURES` unless the configuration uses the SINGLE_STAR dust balance it reproduces.
- Bands `photon_<EUV,FUV,NUV,ONIR,IR>` must match the RT flags exactly when `RADTRANSFER` is on, and need
  `RT_CHEM_PHOTOION`.
- The model's EOS becomes GIZMO's hydro EOS for this gas (T, P, γ come from `jaco_eos`): a model with H only makes the
  hydro treat the gas as pure hydrogen.
- The solver runs T in [1, 1e10] K with no other floor; GIZMO's `MinGasTemp` is applied to the solved state in the
  host, not inside the model. `jaco.check`'s standard grid spans 3 K to 1e9 K.

## 6. Generated code

`python -m jaco.codegen.gizmo.gizmo <model> --language c --ext .cc --output-dir cooling` writes
`microphysics_func_jac.{h,cc}` (RHS and analytic Jacobian of the backward-Euler system; the `SolveVars`/`Params`/
`Outputs` unions, the index enums, the `JACO_HAS_*` macros and the solver metadata: variable kinds, floors, ceilings,
scales, the time-dependent table and the conservation budgets), `microphysics_outputs.cc`, `jaco_eos.cc`
(`jaco_T_to_u`, `jaco_u_to_T`, `jaco_eos_pressure`, the electron abundance), `jaco_interp.h`, `jaco_tables.h`,
`jaco_util.cc` and, for models with 2-D/3-D tables, `jaco_tables.hdf5`, which must sit in the run directory.
`--solve-vars` and `--time-dependent` override the model's own. The solver metadata is checked by `static_assert`s in
`cooling/jaco_solver.cc` (`T` time-dependent with `u_initial`, species from index 2, balance variables in steady
state); most violations are caught earlier, by `Model` itself or by an undeclared `PARAM_u_initial`.

## 7. Checks and tests

In order of cost:

1. `jaco.check(model)` (`model_check.py`): every rate, its partials and the outputs over a `(T, n_Htot, composition)`
   grid with the floating-point semantics of the generated C (clamped `exp`, IEEE underflow, `0 * inf`), the symbolic
   Jacobian against finite differences, undeclared symbols, parameters without defaults and non-finite values.
   `jaco.check(model, "quick")` for a coarse grid. Run it before generating anything.
2. Unit tests next to the model (`models/starforge/tests/`): each rate against its source, and conservation tests
   (`test_dust_heat.py` checks symbolically that two rows cancel; `test_H2_chemical_heat.py` sums rows at fixed states;
   `test_Hminus_detachment.py` at random ones); `c_semantics.py` evaluates expressions as the C does.
3. `tests/test_golden_codegen.py`: hashes of the generated sources for every subpackage of `jaco.models` with a
   `make_model()`; a new model, or an intended change to an old one, regenerates them with
   `UPDATE_GOLDEN=1 python -m pytest tests/test_golden_codegen.py` (slow; all models concurrently).
4. The solver driver in the GIZMO tree, `make -C test/jaco_solver MODEL=<model> test` (it always regenerates; it takes
   `PYTHON=`, not `JACO_PYTHON`). Models named `starforge*` or `sfx_*` get the full driver: sweeps over density,
   temperature, step, PdV work and seed, warm starts, an equilibrium reference check, and recorded cells (states that
   broke the solver in runs, replayed from a `JACO FATAL` dump with `replay`). Every other model gets the generic
   driver: sweeps over density, temperature, step and PdV work only, with every parameter but `n_Htot`, `Δt`, `x_H`,
   `y`, `pdv_work` and `u_initial` at zero. So it exercises the solver, not your parameters.
5. GIZMO tests under `JACO=<model>`: the closed-box ledger (`test/rt_closed_box`, models with bands), the
   thermochemistry benchmarks (`gmc_cooling`, `gmc_cooling_rt`), `second_collapse` for the chemical, radiative and
   thermal bookkeeping, and the full-physics runs. A solve that fails every tier is fatal in production; a new model
   earns its place by surviving these.

## 8. Using a model without GIZMO: the Python/JAX solver

`model.network` is an `EquationSystem` (`equation_system.py`) and `network.solve(knowns, guesses, time_dependent=[],
dt=None, tol=1e-3, careful_steps=20)` solves it with a vectorized Newton iteration in JAX (`numerics/solvers.py`),
one solve per element of the arrays in `knowns` and `guesses`. With `T` in `knowns` it solves the chemistry at that
temperature; without `T` it solves for thermochemical equilibrium (`heat + pdv_work = 0`) as well. Declared parameters
take their defaults when not given. The README's quickstart is this path. What it can and cannot do, measured on
gizmo_integration ae5fbb0 with jax 0.10.1:

- Chemistry at fixed `T`: works, for the README's process list and for a declared model (the toy's `H+` at 1e4 and
  2e4 K: 0.0041 and 0.95).
- Thermochemical equilibrium: works for a model whose rates are closed-form expressions; a `piecewise_powerlaw` or
  `piecewise_linear` table, as in the toy's cooling curve, fails inside the JAX lambdify (`NameError:
  PiecewiseLinearInterp`), as do the 2-D/3-D tables.
- A model with `intermediates` is refused (`NotImplementedError: solve() does not support intermediates; use the
  generated code`). The starforge models have them.
- A backward-Euler step (`time_dependent=["T"]`, `dt=<astropy quantity>`, `u_initial` in `knowns`) fails
  (`NameError: name 'u' is not defined`): the time-dependent path is broken at this commit.
- Results are float32 unless you set `jax.config.update("jax_enable_x64", True)` before importing jaco (the line is
  commented out at the top of `equation_system.py`).

So the Python solver is for equilibrium questions about closed-form networks, as in the README; for a production model,
or for anything time-dependent, the solver is the generated C through the standalone driver (section 7.4) or GIZMO.
`network.solver_functions(solve_vars, time_dependent, return_jac=True)` and `network.generate_code(...)` (README,
"Generating code") give the RHS and Jacobian for another host.

## 9. Style

- One physics per file. A process type gets its own module (H2 formation on grains); families with one functional form
  share one (grain-assisted recombination, WD01, all elements).
- Declare everything: species, parameters, outputs. Never rely on a symbol being inferred.
- Assumptions on the model, never by editing processes: `steady_state`, `fixed`, `derived`, `rules`.
- No model-specific branches in shared code; a variant is `make_model(option=...)` or `model.without`/`replace`.
- Numpy-style docstrings on functions; the process `name` and `bibliography` carry the provenance; comments say what a
  formula is and where it comes from, not what the code does.
- Keep rates finite everywhere `jaco.check` looks: clamp exponents, avoid `0/0` closures, and remember the host compiles
  with `-ffast-math` (no NaN semantics; a non-finite rate is a model bug, not something to absorb).
