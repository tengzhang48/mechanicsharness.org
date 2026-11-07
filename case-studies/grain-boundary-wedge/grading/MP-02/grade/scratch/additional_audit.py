"""Resolve remaining accuracy and monotonicity questions in the written answer."""
import sys
from pathlib import Path
import numpy as np
from scipy.linalg import eig, solve
from kernel import LOperator, S0_analytic
from eig_spectral import build_matrices, Pderivs, phi_vals, make_Gfun

print("CONSTANT-KERNEL QUADRATURE ACCURACY")
z_first = 0.5 * (1 - np.cos(np.pi / 73))
for m in (256, 384, 768, 1536):
    op = LOperator(m)
    z = np.array([0.001, z_first])
    got = op.apply(z, lambda y: 2*y, lambda y: 2*np.ones_like(y))
    print("M", m, "Z", z, "error", got - S0_analytic(z))
    gf, gpf = make_Gfun(0)
    print("phi0 third-derivative transform", op.apply(z, gf, gpf))

print("WEDGE GROWTH RATE FROM THE SOLVER'S EVOLUTION")
a, b, _ = build_matrices(40, 384)
c0, *_ = np.linalg.lstsq(a, np.ones(a.shape[0]), rcond=None)
w, v = eig(-solve(a, b))
amp = solve(v, c0)
z = np.linspace(0.0, 1.0, 1001)
r = 1-z
phi2 = np.empty((len(z), 41))
phi2[:, 0] = -3.75*np.sqrt(r)
phi = np.array([phi_vals(j, z) for j in range(41)]).T
for j in range(1, 41):
    p, p1, p2, *_ = Pderivs(j-1, z)
    phi2[:, j] = 3.75*np.sqrt(r)*p - 5*r**1.5*p1 + r**2.5*p2
for tau in (0.0025, 0.005, 0.01, 0.02, 0.05, 0.1):
    c = (v @ (np.exp(w*tau)*amp)).real
    growth = -phi2 @ c
    sig = phi @ c
    i = np.argmin(growth)
    print("tau", tau, "min dimensionless opening rate", growth[i], "at Z", z[i],
          "tip stress", sig[-1], "max stress", sig.max())

print("POLYNOMIAL FITS VS RECIPROCAL N")
data = np.loadtxt(Path(__file__).parent / "submission/eig_poly.stdout.txt", skiprows=1)
for first in (0, 2, 3, 4):
    d = data[first:]
    for degree in (1, 2):
        for col in (1, 2):
            p = np.polyfit(1/d[:, 0], d[:, col], degree)
            print("first N", d[0, 0], "degree", degree, "lambda", col, "limit", p[-1])
