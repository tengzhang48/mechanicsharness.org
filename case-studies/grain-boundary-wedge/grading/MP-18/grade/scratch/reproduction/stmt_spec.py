#!/usr/bin/env python3
"""Statement-kernel spectrum using the PV-validated hat-basis matrix and a uniform FD D3.

Kernel (as written in the problem):
    S(x,y) = 1/(x-y) - 1/(x+y) - 2 y (x-y)/(x+y)^3 .
The integral operator that multiplies s''' is  M @ s'''  with
    M_ij = (1/pi) Int_0^1 S(x_i,y) phi_j(y) dy     (hat basis phi_j on the uniform grid),
so that  s_tau = M D3 s   (no minus sign; the sign of the spectrum is then physical).
"""
import numpy as np
from numpy.polynomial.legendre import leggauss
from d3_uniform import build_D3_uniform


def F1(x, y):
    """antiderivative in y of the STATEMENT kernel S(x,y)."""
    Lm = np.log(np.abs(y - x) + 1e-300)
    Lp = np.log(np.abs(y + x) + 1e-300)
    # S = 1/(x-y) - 1/(x+y) - 2y(x-y)/(x+y)^3 ; sympy antiderivative
    return (Lm - Lp) - 2.0 * x * (3.0 * y + 2.0 * x) / (y + x) ** 2


def F2(x, y):
    """Int F1 dy (needed for the y-weighting of the hat basis)."""
    Lm = np.log(np.abs(y - x) + 1e-300)
    Lp = np.log(np.abs(y + x) + 1e-300)
    return (-2.0 * x ** 2 + (y + x) * (y * Lm - x * Lm - (y + 7.0 * x) * Lp)) / (y + x)


def _int_lin(x, lo, hi, alpha, beta):
    def G(y):
        return alpha * F1(x, y) + beta * (y * F1(x, y) - F2(x, y))
    return G(hi) - G(lo)


def M_stmt(xcols, nodes):
    N = len(xcols)
    M = np.zeros((N, N))
    for i in range(N):
        z = xcols[i]
        for j in range(N):
            a, b, c = nodes[j], nodes[j + 1], nodes[j + 2]
            val = _int_lin(z, a, b, -a / (b - a), 1.0 / (b - a))
            val += _int_lin(z, b, c, c / (c - b), -1.0 / (c - b))
            M[i, j] = -val / np.pi
    return M


def spectrum(N):
    x = np.linspace(0, 1, N + 2)[1:-1]
    nodes = np.concatenate(([0.0], x, [1.0]))
    M = M_stmt(x, nodes)
    D3 = build_D3_uniform(N)
    A = M @ D3
    ev = np.linalg.eigvals(A)
    return ev[np.argsort(-ev.real)]


if __name__ == "__main__":
    for N in [40, 80, 160, 320, 640]:
        ev = spectrum(N)
        print(f"N={N:4d}  " + "  ".join(f"{v.real:+.5f}{v.imag:+.5f}i" for v in ev[:5]))