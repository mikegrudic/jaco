import sympy as sp
from jaco.math import logistic


def H2_energy_over_kbT(T):
    """Fit for the average energy of a H_2 molecule in units of k_B T"""
    p1 = ((0.0510229560518835 * sp.log(T) - 1.59102802948492) * sp.log(T) + 17.4609732213753) * sp.log(
        T
    ) - 64.5635781542307
    p2 = ((0.102728060235312 * sp.log(T) - 2.42078687298443) * sp.log(T) + 19.4951634944834) * sp.log(
        T
    ) - 51.5605088374882
    return 1.5 + logistic(p1) + logistic(p2)
