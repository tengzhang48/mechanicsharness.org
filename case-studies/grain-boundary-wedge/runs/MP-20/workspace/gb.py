"""
Library for: stress relaxation by diffusion into a grain boundary of a thin film.

Non-dimensional variables:  s = z/h in (0,1),  tau = t/t0,
  t0 = 4 pi kT h^3 / (E* delta_gb D_gb Omega),  g = sigma_gb / sigma_0.

  Elasticity :  g(s,tau) = 1 - (2/pi) int_0^1 S(s,x) u_s(x,tau) dx
                S(a,b) = 1/(a-b) - 1/(a+b) - 2b(a-b)/(a+b)^3     (p.v. at a=b)
  Diffusion  :  u_tau = -(1/2) g_ss                      [g_ss = d^2 g/ds^2]
  BC/IC      :  g(0,tau) = 0,  g_s(1,tau) = g_ss(1,tau) = 0,  g(s,0) = 1

Elimination of u gives the linear initial-value problem in g alone:

  (I)  g_tau = A g ,   A g = (1/pi) int_0^1 S_ss(s,x) g_s?? -> no: see below.

Careful derivation
------------------
u_tau = -(1/2) g_ss is exact (linear elasticity: u is a functional of g).
Take d/dtau of the elasticity relation:
   g_tau(s,tau) = -(2/pi) int_0^1 S(s,x) (u_tau)_s(x,tau) dx
                = -(2/pi) int_0^1 S(s,x) ( -(1/2) g_ss(x,tau) )_s dx
                =  (1/pi) int_0^1 S(s,x) g_sss(x,tau) dx
                =  (1/pi) int_0^1 S_ss(s,x) g_s(x,tau) dx
                    (two integrations by parts in s, boundary terms vanish
                     because S and S_s vanish as x -> infinity and the
                     integrations use only s-derivatives at fixed x),
   i.e.  g_tau = A g  with  A = (1/pi) D^2 int_0^1 S(s,x) . dx.

A is the composition of two symmetric positive-semidefinite operators, hence
A >= 0, and A 1 = 0 exactly (the uniform traction is a steady state).

Eigenproblem:  A f = lambda f,  f'(0) = 0 (because g_tau(0) = 0),
A f automatically has f_s(1) = f_ss(1) = 0 to second order (the kernel and one
derivative vanish at s = 1 for symmetric loads).  The kernel vector f = 1 is
the zero eigenvalue.
"""

import numpy as np
from numpy.polynomial.legendre import leggauss


# ---------------------------------------------------------------- kernels ----
def S_kernel(a, b):
    """S(a,b) = 1/(a-b) - 1/(a+b) - 2b(a-b)/(a+b)^3 (p.v. at a=b, b>0)."""
    with np.errstate(divide="ignore", invalid="ignore"):
        return 1.0 / (a - b) - 1.0 / (a + b) - 2.0 * b * (a - b) / (a + b) ** 3


def S_ss_kernel(a, b):
    """d^2 S/da^2  (sympy-verified)  =
         2/(a-b)^3 - 2/(a+b)^3 - 12b/(a+b)^4 + 48b^2/(a+b)^5
    The first term is understood as the Hadamard finite part at a = b."""
    with np.errstate(divide="ignore", invalid="ignore"):
        return (2.0 / (a - b) ** 3
                - 2.0 / (a + b) ** 3
                - 12.0 * b / (a + b) ** 4
                + 48.0 * b ** 2 / (a + b) ** 5)


def S_ss_diag(b):
    """a -> b limit of the second-derivative kernel:
    2/(a-b)^3 term drops (subtracted), the regular part gives
        -2/(2b)^3 - 12b/(2b)^4 + 48b^2/(2b)^5
      = -1/(4b^3) - 12b/(16b^4) + 48b^2/(32b^5)
      = -1/(4b^3) - 3/(4b^4) + 3/(2b^3)
      =  5/(4b^3) - 3/(4b^4)."""
    return 5.0 / (4.0 * b ** 3) - 3.0 / (4.0 * b ** 4)


# --------------------------------------------------------------- chebyshev ---
def cheb(N):
    """Chebyshev differentiation matrix and CGL nodes on [-1,1]."""
    n = np.arange(N + 1)
    x = np.cos(np.pi * n / N)
    c = np.hstack(([2.0], np.ones(N - 1), [2.0])) * (-1.0) ** n
    X = np.tile(x, (N + 1, 1)).T
    dX = X - X.T
    D = np.outer(c, 1.0 / c) / (dX + np.eye(N + 1))
    D = D - np.diag(D.sum(axis=1))
    return D, x


def lagrange_interp(x_src, x_dst):
    """P[a,j] = L_j(x_dst[a]) for Lagrange basis on the nodes x_src."""
    x_src = np.asarray(x_src, float)
    x_dst = np.asarray(x_dst, float)
    n, m = len(x_src), len(x_dst)
    P = np.ones((m, n))
    for j in range(n):
        for k in range(n):
            if k != j:
                P[:, j] *= (x_dst - x_src[k]) / (x_src[j] - x_src[k])
    return P


# ------------------------------------------------------ discretised operators -
def build(Ng, Nc, s_quad=None):
    """Discretise A on Chebyshev-Gauss-Lobatto nodes of the *interpolant*, with
    the elasticity integral evaluated by Gauss-Legendre quadrature.

        A f |_a = (1/pi) int_0^1 S_ss(s_a, x) f'(x) dx
                ~ (1/pi) sum_b w_b S_ss(s_a, x_b) f'(x_b)
        f'(x_b) ~ sum_j dP[b, j] f_j

    Returns dict with nodes, collocation points, and the matrix
        A_mat[a, j] = (1/pi) sum_b w_b S_ss(s_a,x_b) dP[b,j].
    """
    Dcheb, xn = cheb(Ng)
    order = np.argsort(xn)
    xn = xn[order]
    s_nodes = 0.5 * (xn + 1.0)                  # strictly interior in (0,1)

    xg, wg = leggauss(Nc)
    s_col = 0.5 * (xg + 1.0)                    # strictly interior in (0,1)
    w_col = 0.5 * wg
    if s_quad is not None:
        s_quad = np.asarray(s_quad, float)
        w_col = np.ones(len(s_quad)) / len(s_quad)   # crude fallback
    quad_pts = s_col if s_quad is None else s_quad

    P = lagrange_interp(s_nodes, s_col)         # collocation -> node values
    Pq = lagrange_interp(s_nodes, quad_pts)     # quadrature points -> node values
    # nodes and quadrature points are both already in s, so the Lagrange
    # derivative matrix is directly d/ds.
    dPq = _lagrange_deriv(s_nodes, quad_pts)

    K = S_ss_kernel(s_col[:, None], quad_pts[None, :])
    if s_quad is None:
        np.fill_diagonal(K, S_ss_diag(s_col))
    A_mat = (1.0 / np.pi) * (K * w_col[None, :]) @ dPq

    # elasticity operator:  g = 1 - (2/pi) int S(s,x) u_s dx
    #   ->  u_s(x_b) ~ sum_j dPq[b,j] u_j
    Ke = S_kernel(s_col[:, None], quad_pts[None, :])
    E_gu = -(2.0 / np.pi) * (Ke * w_col[None, :]) @ dPq   # g from u
    return dict(s_nodes=s_nodes, s_col=s_col, w_col=w_col, quad_pts=quad_pts,
                P=P, Pq=Pq, dPq=dPq, A_mat=A_mat, E_gu=E_gu, Ng=Ng, Nc=Nc)


def _lagrange_deriv(x_src, x_dst):
    """dP[a,j] = L_j'(x_dst[a]).  Stable via the ratio-of-products form."""
    x_src = np.asarray(x_src, float)
    x_dst = np.asarray(x_dst, float)
    n, m = len(x_src), len(x_dst)
    dP = np.zeros((m, n))
    for j in range(n):
        for k in range(n):
            if k == j:
                continue
            term = np.ones(m) / (x_src[j] - x_src[k])
            for l in range(n):
                if l != j and l != k:
                    term *= (x_dst - x_src[l]) / (x_src[j] - x_src[l])
            dP[:, j] += term
    return dP


def lagrange_deriv2(x_src, x_dst):
    """d2P[a,j] = L_j''(x_dst[a]).

    L_j(x) = prod_{m!=j} (x-x_m)/(x_j-x_m).  With log-derivative building blocks
    phi_j(x) = sum_{m!=j} 1/(x-x_m) and phi2_j(x) = -sum_{m!=j} 1/(x-x_m)^2
    we have L_j'' = L_j (phi_j^2 + phi2_j).
    """
    x_src = np.asarray(x_src, float)
    x_dst = np.asarray(x_dst, float)
    n, m = len(x_src), len(x_dst)
    d2P = np.zeros((m, n))
    for j in range(n):
        idx = [q for q in range(n) if q != j]
        # L_j at x_dst
        Lj = np.ones(m)
        for k in idx:
            Lj *= (x_dst - x_src[k]) / (x_src[j] - x_src[k])
        # phi_j, phi2_j
        diff = x_dst[:, None] - x_src[None, idx]          # (m, n-1)
        inv = 1.0 / diff
        phi = inv.sum(axis=1)
        phi2 = -(inv ** 2).sum(axis=1)
        d2P[:, j] = Lj * (phi ** 2 + phi2)
    return d2P


# ----------------------------------------------- operator on the node space ----
def operator_nodes(Ng, Nc=None):
    """Return the transport operator A as a square matrix acting on the values
    of the unknown at the Chebyshev-Gauss-Lobatto nodes s_nodes.

    The elasticity-diffusion operator maps node-values to collocation values
    (A_mat), so the node-space operator is  Ahat = P^{-1} A_mat  with
    P = interpolation matrix from nodes to collocation points.
    """
    if Nc is None:
        Nc = Ng + 1
    op = build(Ng, Nc)
    A_mat = op["A_mat"]                       # (Nc, n) node-values -> collocation
    P = op["P"]                               # (Nc, n) node-values -> collocation
    n = Ng + 1
    Ahat = np.linalg.solve(P, A_mat)          # (n, n) node-space operator
    # second-derivative matrix at the nodes (for the eigenvalue equation)
    Dcheb, xn = cheb(Ng)
    order = np.argsort(xn)
    Dcheb = Dcheb[np.ix_(order, order)]
    D2 = (2.0 * Dcheb) @ (2.0 * Dcheb)        # d2/ds2 at the nodes
    op["Ahat"] = Ahat
    op["D2"] = D2
    return op