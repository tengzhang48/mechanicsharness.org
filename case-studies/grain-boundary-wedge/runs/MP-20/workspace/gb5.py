"""
Weak (symmetric) discretisation of the grain-boundary relaxation operator.

    g_tau = A g,     A g(s) = (1/pi) d^2/ds^2 int_0^1 S(s,x) g_x(x) dx

Testing against v with v(0) = 0 and integrating by parts twice (the natural
boundary conditions at s = 0, 1 kill the boundary terms) gives the symmetric
form

    <v, A g> = -(1/pi) int_0^1 int_0^1 v_s(s) S(s,x) g_x(x) dx ds .

With v_j the nodal Lagrange basis on the grid (first node dropped so that
v_j(0) = 0 on the retained subspace is *not* automatically true -- instead we
work on the full nodal space and accept that the s = 0 row is the natural BC
row) the discretisation is

    A_w = -(1/pi) D^T W D ,   D_{ij} = phi_j'(s_i),
    W_{ij} = S(s_i, s_j) w_j

which is exactly symmetric and negative semi-definite.  W is computed with the
p.v. (singularity-subtracted) quadrature validated in gb3.elasticity_matrix.
"""

import numpy as np
from gb2 import nodal_derivatives


def _S(a, x):
    """S(s,x) = 1/(s-x) - 1/(s+x) - 2x(s-x)/(s+x)^3, a 1D array version."""
    return 1.0 / (a - x) - 1.0 / (a + x) - 2.0 * x * (a - x) / (a + x) ** 3


def weight_matrix(nodes, xq, wq):
    """W[i,j] = (2/pi) int S(s_i,x) phi_j(x) dx  (NOT phi_j') -- the p.v. kernel
    against the *value* basis, using singularity subtraction at x = s_i."""
    from gb2 import lagrange_basis
    phi, _, _ = lagrange_basis(nodes)
    Pq = phi(xq)                      # (Nq, n)
    n = len(nodes)
    W = np.zeros((n, n))
    for i, a in enumerate(nodes):
        at = np.abs(xq - a) > 1e-13
        xa, wa = xq[at], wq[at]
        pa = Pq[i, :]                 # phi_j(a)
        # p.v. part of 1/(a-x)
        pv = wa @ ((Pq[at, :] - pa[None, :]) / (a - xa)[:, None])
        cross = pa * np.log(a / (1.0 - a))
        reg = (-Pq[at, :] / (a + xa)[:, None]
               - 2.0 * xa[:, None] * (a - xa)[:, None] * Pq[at, :]
               / (a + xa)[:, None] ** 3)
        W[i, :] = (2.0 / np.pi) * (pv + cross + wa @ reg)
    return W


def weak_operator(nodes, xq, wq):
    D, _ = nodal_derivatives(nodes)
    W = weight_matrix(nodes, xq, wq)
    return -D.T @ W @ D