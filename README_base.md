# reaxion

[![Python package](https://github.com/mikegrudic/reaxion/actions/workflows/test.yml/badge.svg)](https://github.com/mikegrudic/reaxion/actions/workflows/test.yml)
[![Readthedocs Status][docs-badge]][docs-link]
[![codecov](https://codecov.io/github/mikegrudic/reaxion/graph/badge.svg?token=OWJQMWGABZ)](https://codecov.io/github/mikegrudic/reaxion)

[docs-link]:           https://reaxion.readthedocs.io
[docs-badge]:          https://readthedocs.org/projects/reaxion/badge

`reaxion` is a flexible, object-oriented implementation for systems of ISM microphysics and chemistry equations, with numerical solvers implemented in JAX, and interfaces for embedding the equations and their Jacobians into other codes.

## Do we really need yet another ISM code?

`reaxion` might be interesting because it combines two powerful concepts:
1. **Object-oriented implementation of microphysics and chemistry via the `Process` class**, which implements methods for representing physical processes, composing them into a network in a fully-symbolic `sympy` representation. OOP is nice here because if you want to add a new process to `reaxion`, you typically only have to do it in one file. Rate expressions never have to be repeated in-code. Most processes one would want to implement follow very common patterns (e.g. 2-body processes), so class inheritance is also used to minimize new lines of code. 
Once you've constructed your system, `reaxion` can give you the symbolic equations to manipulate and analyze as you please. If you want to solve the equations numerically, `Process` has methods for substituting known values into numerical solvers. It can also automatically generate compilable implementations of the RHS of the system to embed in your choice of simulation code and plug into your choice of solver.
2. **Fast, differentiable implementation of nonlinear algebraic and differential-algebraic equation solvers with JAX**, implemented in its functional programming paradigm (e.g. `reaxion.numerics.newton_rootsolve`). These can achieve excellent numerical throughput running natively on GPUs - in fact, crunching iterates in-place is essentially the best-case application of numerics on GPUs. Differentiability enables sensitivity analysis with respect to all parameters in a single pass, instead of constructing a grid of `N` parameter variations for `N` parameters. This makes it easier in principle to directly answer questions like "How sensitive is this temperature to the abundance of C or the ionization energy of H?", etc.

## Roadmap

`reaxion` is in an early prototyping phase right now. Here are some things I would eventually like to add:
* Flexible implementation of a reduced network suitable for RHD simulations in GIZMO and potentially other codes.
* Dust and radiation physics: add the dust energy equation and evolution of photon number densities to the network.
* Interfaces to convert from other existing chemistry network formats to the `Process` representation.
* Solver robustness upgrades: thermochemical networks can be quite challenging numerically, due to how steeply terms switch on with increasing `T`. In can be hard to get a solution without good initial guesses.
* If possible, glue interface allowing an existing compiled hydro code to call the JAX solvers on-the-fly.

## Installation

Clone the repo and run `pip install .` from the directory.

