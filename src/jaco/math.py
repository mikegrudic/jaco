"""Special mathematical functions"""

import sympy as sp


def logistic(x):
    """Symbolic implementation of the logistic function.

    f(x) = 1/(1+exp(-x))

    Uses tanh which is better behaved numerically in extreme limits
    """
    return 0.5 * (1 + sp.tanh(x / 2))
