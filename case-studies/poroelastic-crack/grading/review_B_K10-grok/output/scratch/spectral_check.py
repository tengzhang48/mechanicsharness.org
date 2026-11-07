#!/usr/bin/env python3
"""Independent finite-crack spectral solve (Chebyshev-Fourier), after MP-38.

This is not a copy: kernels are assembled from the half-plane symbol
derived in derive_kernel discussion, using even Chebyshev bases.

G=L=1. Laplace variable s, q=sqrt(s) in units c=1.
Step-normalized (multiply by s): equivalent to using q as in MP-38.
"""
import json
import numpy as np
from scipy import special, linalg

NU, NU_U = 0.2, 0.4
ETA = (NU_U - NU) / (1 - NU_U)
DCOEF = 2 * (NU_U - NU) / (1 - NU)  # = 2*ETA/(1+ETA)* something
# d_MP38 = 2(nu_u-nu)/(1-nu) = 0.5
# eta_mod = (nu_u-nu)/(1-nu_u)=1/3
# relation: the mechanical prefactor Gk/(2(1-nu)) vs Gk/(2(1-nu_u))
# undrained extra = eta * Gk/(2(1-nu)) * (m-k)/(m+k) * w
# (m-k)/(m+k) = (m-k)^2 / (m^2-k^2) wait m-k = q^2/(m+k), (m-k)/(m+k)= q^2/(m+k)^2
# k (m-k)/(m+k) = k q^2 / (m+k)^2
# Standard: k^3 /(m(m+k)) form comes from rewriting D(k).

print("ETA", ETA, "d_mp38", DCOEF, "1+ETA", 1 + ETA)


def spectral(n, cutoff, spacing, qs):
    k = np.arange(round(cutoff / spacing) + 1) * spacing
    weights = np.ones(len(k))
    weights[1:-1:2] = 4
    weights[2:-1:2] = 2
    weights *= spacing / (3 * np.pi)
    js = np.empty((2 * n, len(k)))
    js[0] = special.j0(k)
    js[1] = special.j1(k)
    low = k < 2 * n + 20
    high = ~low
    for order in range(2, 2 * n):
        js[order, high] = (2 * (order - 1) / k[high]) * js[order - 1, high] - js[order - 2, high]
        js[order, low] = special.jv(order, k[low])
    orders = np.arange(1, 2 * n, 2)
    phi = np.pi * orders[:, None] * js[1::2] / np.where(k == 0, 1.0, k)
    phi[:, 0] = 0.0
    phi[0, 0] = np.pi / 2
    psi = np.pi * js[::2]
    signs = (-1.0) ** np.arange(n)
    tail_a = np.outer(orders * signs, orders * signs)
    tail_c = np.outer(signs, signs)
    diagonal = np.diag(np.pi * orders / 2.0)
    d = DCOEF
    out = []
    for q in qs:
        m = np.sqrt(k * k + q * q)
        aa = -d * k**3 / (m * (m + k))
        bb = d * q * q * k / (m * (m + k))
        cc = d * q * q / m
        a = diagonal + (phi * (weights * aa)) @ phi.T
        a -= d * tail_a / (np.sqrt(cutoff**2 + q * q) + cutoff)
        b = (phi * (weights * bb)) @ psi.T
        c = (psi * (weights * cc)) @ psi.T
        c += d * q * np.arcsinh(q / cutoff) * tail_c
        schur = a if d == 0 else a - b @ linalg.solve(c, b.T, assume_a="pos")
        rhs = np.zeros(n)
        rhs[0] = np.pi * (1 - NU_U) / 2
        u = linalg.solve(schur, rhs, assume_a="pos")
        volume = np.pi * u[0]
        coeff = 2 * q * (volume - np.pi * (1 - NU_U)) / np.sqrt(np.pi)
        # tip intensity using drained COD-K on the Chebyshev opening:
        # w = 2 * sum u_j * sqrt(1-x^2) U_{2j}(x); near x=1,
        # sqrt(1-x^2) U_{2j} ~ (2j+1) sqrt(2(1-x)) * (-1)^j
        # K_drained / sqrt(pi) = sum (-1)^j (2j+1) u_j / (1-nu)
        tip = np.sum(signs * orders * u) / (1 - NU)
        out.append(
            {
                "q": q,
                "volume": float(volume),
                "A_est": float(coeff),
                "K_est": float(tip),
                "delta": 1.0 / q,
            }
        )
    return out


if __name__ == "__main__":
    rows = []
    for n, cutoff, spacing, qs in [
        (16, 8000, 0.5, [20.0, 40.0, 80.0, 160.0]),
        (24, 12000, 0.4, [40.0, 80.0, 160.0, 320.0]),
        (32, 20000, 0.25, [80.0, 160.0, 320.0]),
    ]:
        vals = spectral(n, cutoff, spacing, qs)
        rec = {"n": n, "cutoff": cutoff, "spacing": spacing, "vals": vals}
        rows.append(rec)
        print(json.dumps(rec))
    with open("spectral_check_out.json", "w") as f:
        json.dump(rows, f, indent=2)
