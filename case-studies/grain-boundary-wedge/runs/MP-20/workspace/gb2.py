"""
Second, independent and better-conditioned discretisation of the transport
operator, built by *singularity subtraction*.

Definition used throughout (matching the problem statement literally):

    A f (s) = (1/pi) int_0^1 S_ss(s,x) f'(x) dx ,
    S_ss(a,b) = d^2 S/da^2 = 2/(a-b)^3 - 2/(a+b)^3 - 12b/(a+b)^4 + 48b^2/(a+b)^5

The 2/(a-b)^3 term is a Hadamard finite part.  For a smooth weight phi with
Taylor expansion about x = a it is handled exactly by subtracting the local
part:

    fp int_0^1 phi(x)/(a-x)^3 dx
      = int_0^1 [phi(x) - phi(a) - phi'(a)(x-a)]/(a-x)^3 dx       (regular)
        + phi(a) J3(a) + phi'(a) J2(a)

with the closed forms (verified against high-precision mpmath)

    J3(a) = fp int_0^1 1/(a-x)^3 dx = 1/(2(1-a)^2) - 1/(2a^2)
    J2(a) =    int_0^1 1/(a-x)^2 dx = 1/(1-a) - 1/a
"""

import numpy as np
from numpy.polynomial.legendre import leggauss

from gb import cheb, S_kernel, S_ss_diag


def J3(a):
    """fp int_0^1 1/(a-x)^3 dx."""
    return 1.0 / (2.0 * (1.0 - a) ** 2) - 1.0 / (2.0 * a ** 2)


def J2(a):
    """(ordinary) int_0^1 1/(a-x)^2 dx."""
    return 1.0 / (1.0 - a) - 1.0 / a


def pv_integral_amx3(phi, dphi, a, xq, wq):
    """fp int_0^1 phi(x)/(a-x)^3 dx for a smooth callable phi (vectorised on xq)
    with derivative callable dphi."""
    p = phi(xq)
    dp = dphi(xq)
    pa, dpa = phi(np.array([a]))[0], dphi(np.array([a]))[0]
    at = np.abs(xq - a) > 1e-14
    g = np.zeros_like(xq)
    g[at] = (p[at] - pa - dpa * (xq[at] - a)) / (a - xq[at]) ** 3
    return np.sum(wq * g) + pa * J3(a) + dpa * J2(a)


def Sss_integral(phi, dphi, a, xq, wq):
    """int_0^1 S_ss(a,x) phi(x) dx, using the decomposition
       S_ss = 2/(a-x)^3 - 2/(a+x)^3 - 12x/(a+x)^4 + 48x^2/(a+x)^5."""
    p = phi(xq)
    dp = dphi(xq)
    reg = -2.0 * p / (a + xq) ** 3 - 12.0 * xq * p / (a + xq) ** 4 \
        + 48.0 * xq ** 2 * p / (a + xq) ** 5
    return 2.0 * pv_integral_amx3(phi, dphi, a, xq, wq) + np.sum(wq * reg)


def lagrange_basis(nodes):
    """Return (phi, dphi, d2phi) callables, each mapping a point array x to an
    array of shape (len(x), len(nodes)) holding the Lagrange basis values and
    its first two derivatives on the given nodes."""
    nodes = np.asarray(nodes, float)
    n = len(nodes)

    def phi(x):
        x = np.atleast_1d(np.asarray(x, float))
        V = np.ones((len(x), n))
        for j in range(n):
            for k in range(n):
                if k != j:
                    V[:, j] *= (x - nodes[k]) / (nodes[j] - nodes[k])
        return V

    def dphi(x):
        x = np.atleast_1d(np.asarray(x, float))
        V = np.zeros((len(x), n))
        for j in range(n):
            for k in range(n):
                if k == j:
                    continue
                t = np.ones(len(x)) / (nodes[j] - nodes[k])
                for l in range(n):
                    if l != j and l != k:
                        t *= (x - nodes[l]) / (nodes[j] - nodes[l])
                V[:, j] += t
        return V

    def d2phi(x):
        x = np.atleast_1d(np.asarray(x, float))
        V = np.zeros((len(x), n))
        for j in range(n):
            idx = [q for q in range(n) if q != j]
            base = np.ones(len(x))
            for k in idx:
                base *= (x - nodes[k]) / (nodes[j] - nodes[k])
            diff = x[:, None] - nodes[None, idx]
            with np.errstate(divide="ignore", invalid="ignore"):
                inv = 1.0 / diff
                s1 = inv.sum(axis=1)
                s2 = (inv ** 2).sum(axis=1)
                V[:, j] = base * (s1 ** 2 - s2)
        return V

    return phi, dphi, d2phi


def nodal_derivatives(nodes):
    """Return D1 and D2, the matrices of phi_j'(x_i) and phi_j''(x_i) for the
    Lagrange basis on `nodes`, computed as products of the verified first
    derivative matrix: D2 = D1 @ D1 (both act on nodal values of polynomials
    in the same space, so the composition is exact)."""
    n = len(nodes)
    D1 = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            if i == j:
                D1[i, j] = np.sum(1.0 / (nodes[j] - np.delete(nodes, j)))
            else:
                idx = [m for m in range(n) if m != j and m != i]
                base = np.prod([(nodes[i] - nodes[m]) / (nodes[j] - nodes[m])
                                for m in idx])
                D1[i, j] = base / (nodes[j] - nodes[i])
    return D1, D1 @ D1


def transport_matrix(nodes, xq, wq):
    """Nodal matrix A[i,j] = (1/pi) int_0^1 S_ss(s_i,x) phi_j'(x) dx, where
    phi_j is the Lagrange basis function on `nodes`.  The Hadamard finite part
    of 2/(a-x)^3 acting on phi_j' is handled by subtracting the Taylor part of
    phi_j' at x = a:

      fp int phi_j'(x)/(a-x)^3 dx
        = int [phi_j'(x) - phi_j'(a) - phi_j''(a)(x-a)]/(a-x)^3 dx
          + phi_j'(a) J3(a) + phi_j''(a) J2(a)
    """
    phi, dphi, d2phi = lagrange_basis(nodes)
    n = len(nodes)
    A = np.zeros((n, n))
    DPq = dphi(xq)          # (Nq, n)   phi_j'(x)
    DPn, D2Pn = nodal_derivatives(nodes)
    for i, a in enumerate(nodes):
        at = np.abs(xq - a) > 1e-13
        xa = xq[at]
        wa = wq[at]
        dpa = DPn[i, :]                     # phi_j'(a)
        d2pa = D2Pn[i, :]                   # phi_j''(a)
        num = DPq[at, :] - dpa[None, :] - d2pa[None, :] * (xa - a)[:, None]
        pv = wa @ (num / (a - xa)[:, None] ** 3)
        # regular tails:  -2/(a+x)^3 - 12x/(a+x)^4 + 48x^2/(a+x)^5  times phi_j'
        reg = (-2.0 * DPq[at, :] / (a + xa)[:, None] ** 3
               - 12.0 * xa[:, None] * DPq[at, :] / (a + xa)[:, None] ** 4
               + 48.0 * xa[:, None] ** 2 * DPq[at, :] / (a + xa)[:, None] ** 5)
        tails = wa @ reg
        A[i, :] = (2.0 * pv + tails + 2.0 * dpa * J3(a) + 2.0 * d2pa * J2(a)) / np.pi
    return A