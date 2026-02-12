Turning a jaco symbolic expression into a numerical function
============================================================

``jaco`` abstracts the rates of microphysical processes as symbolic
expressions, but to compute something we must interface this
representation with numerical functions. In this example, we will take
the rate of :math:`H_2` formation as implemented in the
``jaco.models.starforge`` model and create a JAX function to compute it
from array data containing the quantities in the expression.

First we import stuff:

.. code:: ipython3

    from jaco.models.starforge.h2_chemistry.grain_formation import grain_formation
    import sympy as sp
    import numpy as np

The rate per cm^-3 of the process is given by the ``rate`` attribute:

.. code:: ipython3

    grain_formation.rate




.. math::

    \displaystyle \frac{3.0 \cdot 10^{-18} \sqrt{T} Z_{d} f_{d} n_{H}}{\left(1.0 + 10000.0 e^{- \frac{600.0}{Td}}\right) \left(8.0 \cdot 10^{-6} T^{2} + 0.002 T + 0.04 \sqrt{T + Td} + 1.0\right)}



First we must inspect the symbols present in the symbolic object using
``free_symbols``:

.. code:: ipython3

    grain_formation.rate.free_symbols




.. parsed-literal::

    {T, Td, Z_d, f_d, n_H}



These are the quantities we must pass to a function to obtain numerical
values.

We can then construct the numerical function using ``sympy.lambdify``,
specifying JAX as the backend. See the sympy docs for the complete list
of available backends, including e.g. numpy.

.. code:: ipython3

    rate = grain_formation.rate
    symbols = list(sp.ordered(rate.free_symbols))
    rate_lambdified = sp.lambdify(symbols, rate, "jax")
    ?rate_lambdified


.. parsed-literal::

    [0;31mSignature:[0m [0mrate_lambdified[0m[0;34m([0m[0mT[0m[0;34m,[0m [0mTd[0m[0;34m,[0m [0mZ_d[0m[0;34m,[0m [0mf_d[0m[0;34m,[0m [0mn_H[0m[0;34m)[0m[0;34m[0m[0;34m[0m[0m
    [0;31mDocstring:[0m
    Created with lambdify. Signature:
    
    func(T, Td, Z_d, f_d, n_H)
    
    Expression:
    
    3.0e-18*sqrt(T)*Z_d*f_d*n_H/((1.0 + 10000.0*exp(-600.0/Td))*(8.0e-6*T**2 +...
    
    Source code:
    
    def _lambdifygenerated(T, Td, Z_d, f_d, n_H):
        return 3.0e-18*sqrt(T)*Z_d*f_d*n_H/((1.0 + 10000.0*exp(-600.0/Td))*(8.0e-6*T**2 + 0.002*T + 0.04*sqrt(T + Td) + 1.0))
    
    
    Imported modules:
    
    from jax.numpy import sqrt
    from jax.numpy import exp
    [0;31mFile:[0m      ~/code/starforge_tools/microphysics/<lambdifygenerated-3>
    [0;31mType:[0m      function

We see that sympy has created a function whose arguments are the symbols
in the expression. This can be used directly, but we will wrap it in a
function with a ``**kwargs`` interface so that we don’t have to worry
about manually entering the arguments in the correct order:

.. code:: ipython3

    def H2_formation_rate(**kwargs):
        """Returns the H2 formation rate in cm^-3 s^-1
        """    
        return rate_lambdified(kwargs) 





.. parsed-literal::

    Array([2.8339944e-18, 2.8339944e-18, 2.8339944e-18, 2.8339944e-18,
           2.8339944e-18, 2.8339944e-18, 2.8339944e-18, 2.8339944e-18,
           2.8339944e-18, 2.8339944e-18], dtype=float32)



Finally, let’s create some example arrays and plug in the values:

.. code:: ipython3

    num = 10
    T=np.repeat(1.,num)
    Td=np.repeat(1.,num)
    Z_d = np.ones(num)
    f_d = np.ones(num)
    n_H = np.ones(num)
    
    H2_formation_rate(T=T,Td=Td,Z_d=Z_d,f_d=f_d,n_H=n_H)




.. parsed-literal::

    Array([2.8339944e-18, 2.8339944e-18, 2.8339944e-18, 2.8339944e-18,
           2.8339944e-18, 2.8339944e-18, 2.8339944e-18, 2.8339944e-18,
           2.8339944e-18, 2.8339944e-18], dtype=float32)


