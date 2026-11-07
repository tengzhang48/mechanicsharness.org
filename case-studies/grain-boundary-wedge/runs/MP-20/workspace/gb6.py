"""
Finite-volume discretisation of the dimensionless grain-boundary problem

    g(s)  = 1 - int_0^1 S(s,x) w_x(x) dx        (elasticity, dimensionless)
    w_tau = g_ss(s)                             (diffusion)

on a UNIFORM grid s_i = (i+1/2) h, i = 0..N-1 (cell centres) with h = 1/N.
Uniform grids keep every matrix well conditioned; the p.v. elasticity integral
is evaluated by a composite rule whose singular cell is handled analytically.

Unknowns are the cell-average openings w_i.  The elasticity integral is a
matrix-vector product M @ w with

    M[i,j] = int_{cell_j} S(s_i, x) d/dx[ basis_j ] dx

Using a *linear* (finite-volume) reconstruction of w inside each cell,

    w(x) = w_j + (w_{j+1}-w_j)(x - x_j)/h   on cell j (cell edges x_j = j h),

the derivative w_x is piecewise constant = (w_{j+1}-w_j)/h inside cell j, and

    M[i,j] = (1/h) int_{x_j}^{x_{j+1}} S(s_i, x) dx .

The cell integral of S is done analytically: with x_j <= s_i <= x_{j+1} the
p.v. of 1/(s_i - x) contributes log|...| and the regular part is integrated by
high-order Gauss quadrature on the *two half cells* around s_i.
"""

import numpy as np
from numpy.polynomial.legendre import leggauss


def S_scalar(a, x):
    return 1.0 / (a - x) - 1.0 / (a + x) - 2.0 * x * (a - x) / (a + x) ** 3


def elasticity_fv(N):
    """Return the N x N matrix M with g_i = 1 - (M w)_i, cell centres s_i,
    and the cell-centre abscissae s."""
    h = 1.0 / N
    s = (np.arange(N) + 0.5) * h
    xedges = np.arange(N + 1) * h
    xg, wg = leggauss(24)
    M = np.zeros((N, N))
    for i, a in enumerate(s):
        # locate the cell containing a
        for j in range(N):
            lo, hi = xedges[j], xedges[j + 1]
            if lo <= a <= hi:
                # split [lo, hi] at a and integrate each half with Gauss,
                # handling the 1/(a-x) p.v. by subtracting 1/(a-x) on the
                # half-cells: p.v. int_{lo}^{hi} dx/(a-x) = log((hi-a)/(a-lo))
                for (c0, c1) in ((lo, a), (a, hi)):
                    if c1 - c0 <= 0:
                        continue
                    pts = 0.5 * (c1 - c0) * xg + 0.5 * (c1 + c0)
                    wts = 0.5 * (c1 - c0) * wg
                    reg = (-1.0 / (a + pts)
                           - 2.0 * pts * (a - pts) / (a + pts) ** 3)
                    M[i, j] += wts @ reg
                M[i, j] += np.log((hi - a) / (a - lo))
            else:
                pts = 0.5 * (hi - lo) * xg + 0.5 * (hi + lo)
                wts = 0.5 * (hi - lo) * wg
                M[i, j] += wts @ S_scalar(a, pts)
        M[i, :] /= h
    return M, s


def transport_fv(N):
    """g_tau = -A g with A = M @ Gw where Gw maps g-cell-values to w_x via the
    elasticity inverse.  Instead we time-step w directly:
       w_tau = g_ss = -D2 M w   (D2 the standard 3-point second difference,
                                 cell-centred with Neumann ends)
    Returns D2 (N x N) and M."""
    h = 1.0 / N
    D2 = np.zeros((N, N))
    # cell-centred Laplacian: use ghost cells s_{-1}=h/2 (sym) ... implement
    # Neumann at s=1 and Dirichlet g=0 at s=0 through the ghost at s = -h/2.
    # Simplest consistent choice: standard 2nd difference with one-sided ends.
    for i in range(N):
        D2[i, i] = -2.0
        if i > 0:
            D2[i, i - 1] = 1.0
        if i < N - 1:
            D2[i, i + 1] = 1.0
    # ends: reflect (Neumann) at s=1, and at s=0 the traction g=0 fixes the
    # ghost value g_{-1} = -g_0 so that g_ss(0) uses -2 g_0.
    D2[0, 0] = -2.0
    D2[N - 1, N - 1] = -2.0
    D2[N - 1, N - 2] = 2.0
    D2 /= h * h
    return D2, None