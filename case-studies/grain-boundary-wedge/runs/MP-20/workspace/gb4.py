"""
Robust finite-difference discretisation of the grain-boundary relaxation
problem, using a graded grid clustered towards the free surface s = 0 where the
kernel is singular.

Two-field form (no doubly-differentiated kernel, no second derivative of S):

    g(s) = 1 - (2/pi) int_0^1 S(s,x) u_s(x) dx
    u_tau = -(1/2) g_ss

Discretisation
--------------
Grid  s_i = h * rho(i)^p / rho(N)^p,  i = 0..N  (graded, clustered at s=0).
Unknowns: g_i at the grid points and the "flux" variable
    q_i = u_s(s_i)
The elasticity integral is evaluated with a composite quadrature that uses the
grid itself plus the analytic p.v. treatment at the diagonal cell:

    g_i = 1 - (2/pi) [ sum_{cells k != diag} int_{cell_k} S(s_i,x) u_s(x) dx
                       + int_{diag cell} S(s_i,x) u_s(x) dx ]

Within the diagonal cell, u_s is interpolated linearly and S(s_i,x) integrated
analytically (the p.v. of 1/(s_i-x) over an interval symmetric about s_i is
exactly zero, and the regular remainder is smooth).
"""

import numpy as np


def graded_grid(N, p=2.0, n_end=0.0):
    """s_i in [0,1], clustered at s=0 (or at both ends for p>1)."""
    i = np.arange(N + 1)
    t = (i / N) ** p
    return t