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

pls halp.

## Installation

Clone the repo and run `pip install .` from the directory.

# Quickstart: Collisional Ionization Equilibrium

Example of using `reaxion` to solve for collisional ionization equilibrium (CIE) for a hydrogen-helium mixture and plot the ionization states as a function of temperature.


```python
%matplotlib inline
%config InlineBackend.figure_format='retina'
import numpy as np
from matplotlib import pyplot as plt
import sympy as sp
```

## Simple processes
A simple process is defined by a single reaction, with a specified rate.

Let's inspect the structure of a single process, the gas-phase recombination of H+: `H+ + e- -> H + hν` 


```python
from reaxion.processes import CollisionalIonization, GasPhaseRecombination

process = GasPhaseRecombination("H+")
print(f"Name: {process.name}")
print(f"Heating rate coefficient: {process.heat_rate_coefficient}")
print(f"Heating rate per cm^-3: {process.heat}"),
print(f"Rate coefficient: {process.rate_coefficient}")
print(f"Recombination rate per cm^-3: {process.rate}")
print(f"RHS of e- number density equation: {process.network['e-']}")
```

    Name: Gas-phase recombination of H+
    Heating rate coefficient: -1.46719838641439e-26*sqrt(T)/((0.00119216696847702*sqrt(T) + 1.0)**1.748*(0.563615123664978*sqrt(T) + 1.0)**0.252)
    Heating rate per cm^-3: -1.46719838641439e-26*sqrt(T)*n_H+*n_e-/((0.00119216696847702*sqrt(T) + 1.0)**1.748*(0.563615123664978*sqrt(T) + 1.0)**0.252)
    Rate coefficient: 1.41621465870114e-10/(sqrt(T)*(0.00119216696847702*sqrt(T) + 1.0)**1.748*(0.563615123664978*sqrt(T) + 1.0)**0.252)
    Recombination rate per cm^-3: 1.41621465870114e-10*n_H+*n_e-/(sqrt(T)*(0.00119216696847702*sqrt(T) + 1.0)**1.748*(0.563615123664978*sqrt(T) + 1.0)**0.252)
    RHS of e- number density equation: Eq(Derivative(n_e-(t), t), -1.41621465870114e-10*n_H+*n_e-/(sqrt(T)*(0.00119216696847702*sqrt(T) + 1.0)**1.748*(0.563615123664978*sqrt(T) + 1.0)**0.252))


Note that all symbolic representations assume CGS units as is standard in ISM physics.

## Composing processes
Now let's define our full network as a sum of simple processes


```python
processes = [CollisionalIonization(s) for s in ("H", "He", "He+")] + [GasPhaseRecombination(i) for i in ("H+", "He+", "He++")]
system = sum(processes)

system.subprocesses
```




    [Collisional Ionization of H,
     Collisional Ionization of He,
     Collisional Ionization of He+,
     Gas-phase recombination of H+,
     Gas-phase recombination of He+,
     Gas-phase recombination of He++]



Summed processes keep track of all subprocesses, e.g. the total net heating rate is:


```python
system.heat
```




$\displaystyle - \frac{1.55 \cdot 10^{-26} n_{He+} n_{e-}}{T^{0.3647}} - \frac{1.2746917300104 \cdot 10^{-21} \sqrt{T} n_{H} n_{e-} e^{- \frac{157809.1}{T}}}{\frac{\sqrt{10} \sqrt{T}}{1000} + 1} - \frac{1.46719838641439 \cdot 10^{-26} \sqrt{T} n_{H+} n_{e-}}{\left(0.00119216696847702 \sqrt{T} + 1.0\right)^{1.748} \left(0.563615123664978 \sqrt{T} + 1.0\right)^{0.252}} - \frac{9.37661057635428 \cdot 10^{-22} \sqrt{T} n_{He} n_{e-} e^{- \frac{285335.4}{T}}}{\frac{\sqrt{10} \sqrt{T}}{1000} + 1} - \frac{4.9524176975855 \cdot 10^{-22} \sqrt{T} n_{He+} n_{e-} e^{- \frac{631515}{T}}}{\frac{\sqrt{10} \sqrt{T}}{1000} + 1} - \frac{5.86879354565754 \cdot 10^{-26} \sqrt{T} n_{He++} n_{e-}}{\left(0.00119216696847702 \sqrt{T} + 1.0\right)^{1.748} \left(0.563615123664978 \sqrt{T} + 1.0\right)^{0.252}}$



Summing processes also sums all chemical and gas/dust cooling/heating rates. 

## Solving ionization equilibrium

We would like to solve for ionization equilibrium given a temperature $T$, overall H number density $n_{\rm H,tot}$.  We define a dictionary of those input quantities and also one for the initial guesses of the number densities of the species in the reduced network.


```python
Tgrid = np.logspace(3,6,10**6)
ngrid = np.ones_like(Tgrid) * 100

knowns = {"T": Tgrid, "n_Htot": ngrid}

guesses = {
    "H": 0.5*np.ones_like(Tgrid),
    "He": 1e-5*np.ones_like(Tgrid),
    "He+": 1e-5*np.ones_like(Tgrid)
}
```

Note that by default, the solver only directly solves for $n_{\rm H}$, $n_{\rm He}$ and $n_{\rm He+}$ because $n_{\rm H+}$, $n_{\rm He++}$, and $n_{\rm e-}$ are eliminated by conservation equations. So we only need initial guesses for those 3 quantities. By default the solver takes abundances $x_i = n_i / n_{\rm H,tot}$ as inputs and outputs.

The `solve` method calls the JAX solver and computes the solution:


```python
sol = system.solve(knowns, guesses,tol=1e-3)
print(sol)
```

    {'H': Array([9.9999994e-01, 9.9999994e-01, 9.9999994e-01, ..., 6.0612069e-07,
           6.0611501e-07, 6.0610915e-07], dtype=float32), 'He': Array([9.2546351e-02, 9.2546351e-02, 9.2546351e-02, ..., 2.7493625e-09,
           2.7493032e-09, 2.7492439e-09], dtype=float32), 'He+': Array([3.1222429e-13, 3.1222396e-13, 3.1222412e-13, ..., 7.6922215e-06,
           7.6921306e-06, 7.6920396e-06], dtype=float32), 'He++': Array([0.        , 0.        , 0.        , ..., 0.09253865, 0.09253865,
           0.09253865], dtype=float32), 'H+': Array([5.9604645e-08, 5.9604645e-08, 5.9604645e-08, ..., 9.9999940e-01,
           9.9999940e-01, 9.9999940e-01], dtype=float32), 'e-': Array([5.9604957e-08, 5.9604957e-08, 5.9604957e-08, ..., 1.1850843e+00,
           1.1850843e+00, 1.1850843e+00], dtype=float32)}



```python
for i, xi in sorted(sol.items()):
    plt.loglog(Tgrid, xi, label=i)
plt.legend(labelspacing=0)
plt.ylabel("$x_i$")
plt.xlabel("T (K)")
plt.ylim(1e-4,3)
```




    (0.0001, 3)




    
![png](CIE_files/CIE_15_1.png)
    

