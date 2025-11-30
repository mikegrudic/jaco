"""Implementation of models for dust temperature, for when RT info is not available."""

# note: ISRF model can be adopted from Zucconi 2001, as an up-to-date version of the original Black multi-blackbody parameterization
import sympy as sp
from .symbols import A_V, ISRF, T_CMB


def dust_temperature_prescription(model="Hocuk 2017"):
    match model:
        case "Hocuk 2017":  # 2017A&A...604A..58H
            return (11 + 5.7 * sp.tanh(0.61 - sp.log(A_V, 10))) * ISRF ** (1.0 / 5.9)
        case "Hollenbach 1991":  # 1991ApJ...377..192H
            nu_0 = 3e15
            T0 = 12.17 * ISRF ** (1 / 5)
            tau_100 = 2.7e3 * ISRF * T0**-5
            return (
                8.9e-11 * nu_0 * ISRF * sp.exp(-1.8 * A_V)
                + T_CMB**5
                + 3.4e-2 * (0.42 - sp.log(3.5e-2 * tau_100 * T0)) * tau_100 * T0**6
            ) ** 0.2
        case _:
            raise NotImplementedError(f"Dust temeprature model {model} not implemented.")
