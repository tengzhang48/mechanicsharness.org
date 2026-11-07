"""Check whether the claim of pointwise monotone wedge growth is numerical noise."""
import numpy as np
from scipy.linalg import eig, solve
from eig_spectral import build_matrices, Pderivs

for n in (24, 56, 72):
    a, b, _ = build_matrices(n, 768)
    c0, *_ = np.linalg.lstsq(a, np.ones(a.shape[0]), rcond=None)
    w, v = eig(-solve(a, b))
    amp = solve(v, c0)
    z = np.array([0.889, 0.944])
    r = 1-z
    phi2 = np.empty((2, n+1))
    phi2[:, 0] = -3.75*np.sqrt(r)
    for j in range(1, n+1):
        p, p1, p2, *_ = Pderivs(j-1, z)
        phi2[:, j] = 3.75*np.sqrt(r)*p - 5*r**1.5*p1 + r**2.5*p2
    for tau in (0.0025, 0.005, 0.1):
        c = (v @ (np.exp(w*tau)*amp)).real
        cint = (v @ (np.expm1(w*tau)/w*amp)).real
        print("N", n, "tau", tau, "Z", z, "dimensionless opening rate", -phi2 @ c,
              "dimensionless opening", -phi2 @ cint,
              "late-mode amplitude estimate", c[0]*np.exp(11.742668071309005*tau))
