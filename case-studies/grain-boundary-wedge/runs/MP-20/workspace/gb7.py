"""
Galerkin (symmetric, well-conditioned) discretisation of the dimensionless
grain-boundary relaxation operator.

Dimensionless system:
    g(s)  = 1 - int_0^1 S(s,x) w_x(x) dx            (elasticity)
    w_tau = g_ss(s)                                  (diffusion)
with g(0) = 0, g_s(1) = 0.

Writing A1 for the symmetric operator (A1 w)(s) = int_0^1 S(s,x) w_x(x) dx,
    g = 1 - A1 w  =>  g_tau = -A1 w_tau = -A1 (D^2 g) = -(A1 D^2) g .

Galerkin with basis {chi_j} (chi_j the nodal Lagrange functions on an interior
grid, so chi_j(0) = 0 automatically for j >= 1 and gauge fixed by dropping the
constant):
    M_{ij} = (chi_i, chi_j),    K_{ij} = (chi_i, A1 D^2 chi_j) = -(chi_i', (A1 chi_j)')
so that  M c' = -K c  and the generalised eigenproblem is  K c = lambda M c .

Both M and K are built from FIRST derivatives of the (validated) elasticity
operator, hence well conditioned.
"""

import numpy as np
from numpy.polynomial.legendre import leggauss
from gb2 import lagrange_basis, nodal_derivatives


def _A1_mat(nodes, xq, wq, basis='value'):
    """Matrix A1 with (A1 u)(node_i) = int_0^1 S(node_i,x) sum_j phi_j'(x) u_j dx,
    singularity subtracted.  Same quadrature as gb3.elasticity_matrix but
    without the 2/pi (the dimensionless kernel already has it absorbed)."""
    _, dphi, _ = lagrange_basis(nodes)
    DPq = dphi(xq)
    DPn, _ = nodal_derivatives(nodes)
    n = len(nodes)
    A = np.zeros((n, n))
    for i, a in enumerate(nodes):
        at = np.abs(xq - a) > 1e-13
        xa, wa = xq[at], wq[at]
        dpa = DPn[i, :]
        pv = wa @ ((DPq[at, :] - dpa[None, :]) / (a - xa)[:, None])
        cross = dpa * np.log(a / (1.0 - a))
        reg = (-DPq[at, :] / (a + xa)[:, None]
               - 2.0 * xa[:, None] * (a - xa)[:, None] * DPq[at, :]
               / (a + xa)[:, None] ** 3)
        A[i, :] = pv + cross + wa @ reg
    return A


def galerkin_operator(nodes, nq=800):
    """Return M (mass) and K (= -D1^T A1 D1) for K c = lambda M c.
    K uses (chi_i', (A1 chi_j)') = D1^T A1 D1 with the *first-derivative* form,
    which is automatically symmetric if A1 is (it is, S symmetric)."""
    xg, wg = leggauss(nq)
    xq, wq = 0.5 * (xg + 1.0), 0.5 * wg
    D1, _ = nodal_derivatives(nodes)
    A1 = _A1_mat(nodes, xq, wq)
    A1 = 0.5 * (A1 + A1.T)          # symmetrise the quadrature asymmetry
    K = -D1.T @ A1 @ D1
    K = 0.5 * (K + K.T)
    return K