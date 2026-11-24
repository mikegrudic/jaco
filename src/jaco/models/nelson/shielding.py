"""CO self-shielding and shielding by H_2 of the CO photodissociation, from Table 5 of van Dishoeck & Black (1988,
ApJ 334, 771), the source NL99 takes its S(N(CO), N(H_2)) from (through Bergin, Langer & Goldsmith 1995, who use it).

Theta[N(CO), N(H_2)] gives the reduction of the 12CO photodissociation rate due to line absorption in the unmodified
Draine (1978) field. The model interpolates bilinearly in ln Theta on the log10 column grid and clamps each column
independently to its tabulated range. Theta is 1 when both columns are at or below 1 cm^-2. H_2 can still shield CO
when only the CO column is at its lower edge. The reaction rate applies dust attenuation separately as exp(-3 A_V).
"""

import numpy as np
import sympy as sp

LOG_N_CO = np.array([0.0, 13, 14, 15, 16, 17, 18, 19])
LOG_N_H2 = np.array([0.0, 19, 20, 21, 22, 23])

# Table 5, 12CO block: rows log N(H_2) = 0, 19, ..., 23; columns log N(CO) = 0, 13, ..., 19
THETA = np.array([
    [   1.000, 9.681e-1, 7.764e-1, 3.631e-1, 7.013e-2, 1.295e-2, 1.738e-3, 9.985e-5],
    [8.215e-1, 7.916e-1, 6.160e-1, 2.749e-1, 5.351e-2, 1.065e-2, 1.519e-3, 8.818e-5],
    [7.160e-1, 6.900e-1, 5.360e-1, 2.359e-1, 4.416e-2, 8.769e-3, 1.254e-3, 7.558e-5],
    [3.500e-1, 3.415e-1, 2.863e-1, 1.360e-1, 2.500e-2, 4.983e-3, 7.151e-4, 3.796e-5],
    [4.973e-2, 4.877e-2, 4.296e-2, 2.110e-2, 4.958e-3, 9.245e-4, 1.745e-4, 8.377e-6],
    [1.310e-4, 1.293e-4, 1.160e-4, 6.346e-5, 1.822e-5, 6.842e-6, 3.622e-6, 3.572e-7],
])


def _log10(x):
    return sp.log(x) / sp.log(10)


def _interp1(x, X, Y):
    """Piecewise-linear interpolant through (X, Y), held at the end values outside, as a plain sp.Piecewise: unlike
    jaco's table interpolants it also evaluates in numpy and JAX (Model.solve)"""
    pieces = [(float(Y[0]), x <= float(X[0]))]
    for k in range(len(X) - 1):
        slope = float((Y[k + 1] - Y[k]) / (X[k + 1] - X[k]))
        pieces.append((float(Y[k]) + slope * (x - float(X[k])), x <= float(X[k + 1])))
    pieces.append((float(Y[-1]), True))
    return sp.Piecewise(*pieces)


def _hat(y, nodes, i):
    """Weight of node i of the piecewise-linear interpolant on `nodes` at y (y inside the range of the nodes)"""
    last = len(nodes) - 1
    right = (nodes[min(i + 1, last)] - y) / (nodes[min(i + 1, last)] - nodes[i]) if i < last else 1
    left = (y - nodes[max(i - 1, 0)]) / (nodes[i] - nodes[max(i - 1, 0)]) if i > 0 else 1
    return sp.Max(0, sp.Min(left, right))


def co_shielding(N_CO, N_H2):
    """Theta[N(CO), N(H_2)] of van Dishoeck & Black (1988) Table 5 for column densities (cm^-2, expressions)"""
    x = _log10(sp.Max(N_CO, 1))
    y = sp.Min(sp.Max(_log10(sp.Max(N_H2, 1)), float(LOG_N_H2[0])), float(LOG_N_H2[-1]))
    rows = [_interp1(x, LOG_N_CO, np.log(THETA[i])) for i in range(len(LOG_N_H2))]
    return sp.exp(sum(_hat(y, [float(v) for v in LOG_N_H2], i) * rows[i] for i in range(len(rows))))
