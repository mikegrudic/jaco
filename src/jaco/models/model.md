# Implementing a model in `jaco`

A model consists of a set of physical processes along with a set of assumptions that are made for various symbols.

## Step 1: Implement the processes

The guiding principle for implementing chemistry and microphysics in `jaco` is to have an experience more-or-less opposite of your advisor's advisor's single-file, 8000-line, incomprehensible Fotran 77 code of doom that solves all of reality in one function. You know what I'm talking about. A model implementation should be: 

1. Modular
Ideally, you should keep your implementation as modular as possible. Typically, a certain type of process belongs in its own `.py` file, but naturally there will sometimes be groups of processes that are similar enough in their form that they should go together. For example, formation of $H_2$ on the surface of dust grains is a pretty special reaction with a particular implementation (e.g. Hollenbach & McKee 1979), and thus should probably have its own file. However e.g. grain-assisted recombination of ions can and should be implemented in one place because the available fitting functions all have the same form and differ only in their parameters (e.g. Weingartner & Draine 2001)

2. Extensible
Implement things in as general a way as is reasonable, to save time in case you or somebody else wants to extend the model.

3. Documented
It should be easy for an interested colleague to know exactly what is in the model and where its data or models came from. When instantiating the process, please supply a bibcode `bibliography` argument, and a descriptive name of the process in the `name` argument. The goal is to lose none of this information as we assemble the network, so that it is easy to recall what exactly is in the network, and where it came from. Please numpy-style docstrings on any functions.

## Step 2: Sum the process to obtain the full system

## Step 3: Implement the assumptions

### Style guidelines
