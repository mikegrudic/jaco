Quickstart: Collisional Ionization Equilibrium
==============================================

Example of using ``pism`` to solve for collisional ionization
equilibrium (CIE) for a hydrogen-helium mixture and plot the ionization
states as a function of temperature.

.. code:: ipython3

    %matplotlib inline
    %config InlineBackend.figure_format='retina'
    import numpy as np
    from matplotlib import pyplot as plt
    import sympy as sp

Simple processes
----------------

A simple process is defined by a single reaction, with a specified rate.

Let’s inspect the structure of a single process, the gas-phase
recombination of H+: ``H+ + e- -> H + hν``

.. code:: ipython3

    from pism.processes import CollisionalIonization, GasPhaseRecombination
    
    process = GasPhaseRecombination("H+")
    print(f"Name: {process.name}")
    print(f"Heating rate coefficient: {process.heat_rate_coefficient}")
    print(f"Heating rate per cm^-3: {process.heat}"),
    print(f"Rate coefficient: {process.rate_coefficient}")
    print(f"Recombination rate per cm^-3: {process.rate}")
    print(f"RHS of e- number density equation: {process.network['e-']}")


.. parsed-literal::

    Name: Gas-phase recombination of H+
    Heating rate coefficient: -1.46719838641439e-26*sqrt(T)/((0.00119216696847702*sqrt(T) + 1.0)**1.748*(0.563615123664978*sqrt(T) + 1.0)**0.252)
    Heating rate per cm^-3: -1.46719838641439e-26*sqrt(T)*n_H+*n_e-/((0.00119216696847702*sqrt(T) + 1.0)**1.748*(0.563615123664978*sqrt(T) + 1.0)**0.252)
    Rate coefficient: 1.41621465870114e-10/(sqrt(T)*(0.00119216696847702*sqrt(T) + 1.0)**1.748*(0.563615123664978*sqrt(T) + 1.0)**0.252)
    Recombination rate per cm^-3: 1.41621465870114e-10*n_H+*n_e-/(sqrt(T)*(0.00119216696847702*sqrt(T) + 1.0)**1.748*(0.563615123664978*sqrt(T) + 1.0)**0.252)
    RHS of e- number density equation: Eq(Derivative(n_e-(t), t), -1.41621465870114e-10*n_H+*n_e-/(sqrt(T)*(0.00119216696847702*sqrt(T) + 1.0)**1.748*(0.563615123664978*sqrt(T) + 1.0)**0.252))


Note that all symbolic representations assume CGS units as is standard
in ISM physics.

Composing processes
-------------------

Now let’s define our full network as a sum of simple processes

.. code:: ipython3

    processes = [CollisionalIonization(s) for s in ("H", "He", "He+")] + [GasPhaseRecombination(i) for i in ("H+", "He+", "He++")]
    system = sum(processes)
    
    system.subprocesses




.. parsed-literal::

    [Collisional Ionization of H,
     Collisional Ionization of He,
     Collisional Ionization of He+,
     Gas-phase recombination of H+,
     Gas-phase recombination of He+,
     Gas-phase recombination of He++]



Summed processes keep track of all subprocesses, e.g. the total net
heating rate is:

.. code:: ipython3

    system.heat




.. math::

    \displaystyle - \frac{1.55 \cdot 10^{-26} n_{He+} n_{e-}}{T^{0.3647}} - \frac{1.2746917300104 \cdot 10^{-21} \sqrt{T} n_{H} n_{e-} e^{- \frac{157809.1}{T}}}{\frac{\sqrt{10} \sqrt{T}}{1000} + 1} - \frac{1.46719838641439 \cdot 10^{-26} \sqrt{T} n_{H+} n_{e-}}{\left(0.00119216696847702 \sqrt{T} + 1.0\right)^{1.748} \left(0.563615123664978 \sqrt{T} + 1.0\right)^{0.252}} - \frac{9.37661057635428 \cdot 10^{-22} \sqrt{T} n_{He} n_{e-} e^{- \frac{285335.4}{T}}}{\frac{\sqrt{10} \sqrt{T}}{1000} + 1} - \frac{4.9524176975855 \cdot 10^{-22} \sqrt{T} n_{He+} n_{e-} e^{- \frac{631515}{T}}}{\frac{\sqrt{10} \sqrt{T}}{1000} + 1} - \frac{5.86879354565754 \cdot 10^{-26} \sqrt{T} n_{He++} n_{e-}}{\left(0.00119216696847702 \sqrt{T} + 1.0\right)^{1.748} \left(0.563615123664978 \sqrt{T} + 1.0\right)^{0.252}}



Summing processes also sums all chemical and gas/dust cooling/heating
rates.

Solving ionization equilibrium
------------------------------

We would like to solve for ionization equilibrium given a temperature
:math:`T`, overall H number density :math:`n_{\rm H,tot}`. We define a
dictionary of those input quantities and also one for the initial
guesses of the number densities of the species in the reduced network.

.. code:: ipython3

    Tgrid = np.logspace(3,6,10**6)
    ngrid = np.ones_like(Tgrid) * 100
    
    knowns = {"T": Tgrid, "n_Htot": ngrid}
    
    guesses = {
        "H": 0.5*np.ones_like(Tgrid),
        "He": 1e-5*np.ones_like(Tgrid),
        "He+": 1e-5*np.ones_like(Tgrid)
    }

Note that by default, the solver only directly solves for
:math:`n_{\rm H}`, :math:`n_{\rm He}` and :math:`n_{\rm He+}` because
:math:`n_{\rm H+}`, :math:`n_{\rm He++}`, and :math:`n_{\rm e-}` are
eliminated by conservation equations. So we only need initial guesses
for those 3 quantities. By default the solver takes abundances
:math:`x_i = n_i / n_{\rm H,tot}` as inputs and outputs.

The ``solve`` method calls the JAX solver and computes the solution:

.. code:: ipython3

    sol = system.solve(knowns, guesses,tol=1e-3)
    print(sol)


.. parsed-literal::

    {'He': Array([9.2546351e-02, 9.2546351e-02, 9.2546351e-02, ..., 2.7493625e-09,
           2.7493037e-09, 2.7492442e-09], dtype=float32), 'He+': Array([3.1916630e-13, 3.1095203e-13, 3.1464317e-13, ..., 7.6922206e-06,
           7.6921306e-06, 7.6920396e-06], dtype=float32), 'H': Array([9.9999994e-01, 9.9999994e-01, 9.9999994e-01, ..., 6.0612075e-07,
           6.0611501e-07, 6.0610921e-07], dtype=float32), 'He++': Array([0.        , 0.        , 0.        , ..., 0.09253865, 0.09253865,
           0.09253865], dtype=float32), 'H+': Array([5.9604645e-08, 5.9604645e-08, 5.9604645e-08, ..., 9.9999940e-01,
           9.9999940e-01, 9.9999940e-01], dtype=float32), 'e-': Array([5.9604965e-08, 5.9604957e-08, 5.9604957e-08, ..., 1.1850843e+00,
           1.1850843e+00, 1.1850843e+00], dtype=float32)}


.. code:: ipython3

    for i, xi in sorted(sol.items()):
        plt.loglog(Tgrid, xi, label=i)
    plt.legend(labelspacing=0)
    plt.ylabel("$x_i$")
    plt.xlabel("T (K)")
    plt.ylim(1e-4,3)




.. parsed-literal::

    (0.0001, 3)




.. image:: CIE_files/CIE_15_1.png
   :width: 761px
   :height: 729px

