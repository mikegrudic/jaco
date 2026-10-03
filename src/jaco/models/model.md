# Implementing a model in `jaco`

A model consists of a set of physical processes along with a set of assumptions that are made for various symbols.

## Step 1: Implement the processes

The guiding principle for implementing chemistry and microphysics in `jaco` is to have an experience more-or-less opposite of your advisor's advisor's single-file, 8000-line, incomprehensible Fotran 77 code of doom that solves all of reality in one function. You know what I'm talking about. A model implementation should be: 

1. Modular
Ideally, you should keep your implementation as modular as possible. Typically, a certain type of process belongs in its own `.py` file, but naturally there will sometimes be groups of processes that are similar enough in their form that they should go together. For example, formation of $H_2$ on the surface of dust grains is a pretty special reaction with a particular implementation (e.g. Hollenbach & McKee 1979), and thus should probably have its own file. However e.g. grain-assisted recombination of ions can and should be implemented in one place because the available fitting functions all have the same form and differ only in their parameters (e.g. Weingartner & Draine 2001).

2. Extensible
Implement things in as general a way as is reasonable, to save time in case you or somebody else wants to extend the model.

3. Documented
It should be easy for an interested colleague to know exactly what is in the model and where its data or models came from. When instantiating the process, please supply a bibcode `bibliography` argument, and a descriptive name of the process in the `name` argument. The goal is to lose none of this information as we assemble the network, so that it is easy to recall what exactly is in the network, and where it came from. Please numpy-style docstrings on any functions.

Processes are immutable: write a reaction as `Reaction("H + H -> H_2", k, heat_per_reaction=..., name=..., bibliography=[...])`
(clumping `C_k` is applied by default to `k >= 2` reactants; pass `clumping=` to say otherwise) and a heating or cooling
term as `ThermalTerm(heat, ...)`. Give every process a unique name: it is the process's id in the model.

## Step 2: Collect the processes in a `Model`

`make_model()` returns `Model(processes, solve_vars=..., time_dependent=...)`. A composite process is split into its
atoms. `model + process` adds one (a duplicate id raises), `model.without(id)` and `model.replace(id, new)` take one out
or swap it, and `model + other_model` merges two (conflicting declarations raise).

## Step 3: Implement the assumptions

Declare them on the model, never by editing processes: `steady_state=["H-"]` closes a species from its own rate
equation when the network is assembled; `fixed={species: expression}` prescribes an abundance; `derived={parameter:
expression}` substitutes an expression of the solve variables for a parameter before differentiation; `rules=[Rule(...)]`
rewrite the processes (e.g. a different density in the rates) at assembly, so processes added later get them too.

## Step 4: Declare the contract

`species=[Species("H"), ..., Species("photon_EUV", "radiation"), Species("dust heat", "energy")]` lists every species
with its kind: material species (composition parsed from the name) make up the EOS and the conservation sums, trace
species stay out of both, radiation is not a collider for the default clumping. `parameters=[Parameter("G_0", "Habing",
1.0, "FUV field"), ...]` declares every input besides the core ones (`n_Htot`, `pdv_work`, `y`, ...) and those the
species imply (abundances, element totals, start-of-step values). Code generation refuses any undeclared symbol, so a
typo or a stray `C_4` fails at codegen instead of becoming a struct field the host zero-fills. The generated header
defines `JACO_HAS_VAR_<name>` and `JACO_HAS_PARAM_<name>` for each field, so host code fills what the model has under
`#ifdef`.

### Style guidelines
