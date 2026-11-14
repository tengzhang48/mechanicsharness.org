#!/usr/bin/env python3
"""Solve the mixed-mode, coplanar crack-front perturbation problem.

The dimensional quantities cancel in the linearized Rice/Gao--Rice response.
Only Poisson's ratio, the two straight-front SIF ratios, and q enter.
"""

import json
from pathlib import Path

import numpy as np
from scipy.optimize import brentq


NU = 0.3
R2 = 0.6
R3 = 0.8
Q_TOUGH = 0.5
N_MODES = 256
OUT = Path(__file__).with_name("results.json")

CI = 1.0 - NU**2
CIII = 1.0 + NU
C2 = 2.0 / (2.0 - NU)
C3 = 2.0 * (1.0 - NU) / (2.0 - NU)
B3 = (2.0 - NU) / 2.0

# G0 = D*K10**2/E.  The energy response to a(s) is
# L[f] = -P*|D_s|f - B*f'.
D = CI * (1.0 + R2**2) + CIII * R3**2
W = CI * (1.0 + C2 * (R2**2 + R3**2))
P = 2.0 * np.pi * W / D
B = R2 * R3 * (3.0 * CI - (1.0 + NU) * (2.0 - NU)) / D

# g(s) = 2*sum(q**m*cos(2*pi*m*s)).  If
# f = fc*cos + fs*sin, equilibrium gives the two equations
# P*fc + 2*pi*B*fs = -2*q**m/m,
# P*fs - 2*pi*B*fc = 0.
M = np.arange(1, N_MODES + 1, dtype=float)
G = 2.0 * Q_TOUGH ** M
DEN = P**2 + (2.0 * np.pi * B) ** 2
FC = -P * G / (M * DEN)
FS = -2.0 * np.pi * B * G / (M * DEN)


def evaluate(s):
    """Return f, its directional operators, kappas, and energy residuals."""
    scalar = np.ndim(s) == 0
    s = np.atleast_1d(np.asarray(s, dtype=float))
    angle = 2.0 * np.pi * s[:, None] * M[None, :]
    co = np.cos(angle)
    si = np.sin(angle)

    f = co @ FC + si @ FS
    absd_f = co @ (M * FC) + si @ (M * FS)
    fp = 2.0 * np.pi * (-si @ (M * FC) + co @ (M * FS))
    fpp = -(2.0 * np.pi) ** 2 * (co @ (M**2 * FC) + si @ (M**2 * FS))
    absd_fp = 2.0 * np.pi * (-si @ (M**2 * FC) + co @ (M**2 * FS))

    kappa1 = -np.pi * absd_f
    kappa2 = -np.pi * C2 * R2 * absd_f - 1.5 * R3 * fp
    kappa3 = -np.pi * C3 * R3 * absd_f + B3 * R2 * fp
    kappa3_prime = -np.pi * C3 * R3 * absd_fp + B3 * R2 * fpp

    # Two independently assembled expressions for delta(G)/G0.
    energy_from_k = 2.0 * (
        CI * (kappa1 + R2 * kappa2) + CIII * R3 * kappa3
    ) / D
    energy_from_operator = -P * absd_f - B * fp

    result = {
        "f": f,
        "fp": fp,
        "kappa1": kappa1,
        "kappa2": kappa2,
        "kappa3": kappa3,
        "kappa3_prime": kappa3_prime,
        "energy_from_k": energy_from_k,
        "energy_from_operator": energy_from_operator,
    }
    if scalar:
        result = {key: value[0] for key, value in result.items()}
    return result


def toughness_g(s):
    s = np.asarray(s)
    return (1.0 - Q_TOUGH**2) / (
        1.0 - 2.0 * Q_TOUGH * np.cos(2.0 * np.pi * s) + Q_TOUGH**2
    ) - 1.0


def stationary_points(component, grid_size=4097):
    """Locate all sign-changing stationary points on [-1/2, 1/2]."""
    grid = np.linspace(-0.5, 0.5, grid_size)
    values = evaluate(grid)[component]
    roots = []
    for i in range(grid_size - 1):
        if values[i] == 0.0:
            roots.append(grid[i])
        elif values[i] * values[i + 1] < 0.0:
            roots.append(
                brentq(
                    lambda x: evaluate(x)[component],
                    grid[i],
                    grid[i + 1],
                    xtol=1.0e-14,
                    rtol=1.0e-14,
                )
            )
    # Include the seam as a candidate even if no sample landed exactly on it.
    roots.append(-0.5)
    unique = []
    for root in sorted(roots):
        if not unique or abs(root - unique[-1]) > 1.0e-9:
            unique.append(root)
    return np.array(unique)


def part_a():
    result = {
        "C1": -np.pi,
        "S1": 0.0,
        "C2": -np.pi * C2 * R2,
        "S2": 2.0 * np.pi * 1.5 * R3,
        "C3": -np.pi * C3 * R3,
        "S3": -2.0 * np.pi * B3 * R2,
        "CG": -P,
        "SG": 2.0 * np.pi * B,
    }
    # Check the energy coefficients directly from all six SIF coefficients.
    cg_direct = 2.0 * (
        CI * (result["C1"] + R2 * result["C2"])
        + CIII * R3 * result["C3"]
    ) / D
    sg_direct = 2.0 * (
        CI * R2 * result["S2"] + CIII * R3 * result["S3"]
    ) / D
    np.testing.assert_allclose(result["CG"], cg_direct, rtol=0.0, atol=2e-15)
    np.testing.assert_allclose(result["SG"], sg_direct, rtol=0.0, atol=2e-15)
    return result


def hilbert_quadrature_check():
    """Check R[cos(2*pi*s)] = -pi*cos(2*pi*s) by periodic PV quadrature."""
    s0 = 0.0
    # On a period, summing the inverse-square kernel over its images gives
    # pi**2*csc(pi*(s'-s))**2.  Midpoints avoid the removable singular point;
    # their symmetry about s0 also makes the near-singular contributions cancel.
    n = 4096
    sp = (np.arange(n) + 0.5) / n - 0.5
    integrand = (
        0.5
        * np.pi
        * (np.cos(2.0 * np.pi * sp) - np.cos(2.0 * np.pi * s0))
        / np.sin(np.pi * (sp - s0)) ** 2
    )
    numerical = float(np.mean(integrand))
    expected = -np.pi * np.cos(2.0 * np.pi * s0)
    np.testing.assert_allclose(numerical, expected, rtol=0.0, atol=2e-12)


def solve():
    a = part_a()
    samples = np.arange(16, dtype=float) / 16.0
    v = evaluate(samples)

    f_turning = stationary_points("fp")
    f_values = evaluate(f_turning)["f"]
    f_min = float(np.min(f_values))
    f_max = float(np.max(f_values))

    k3_turning = stationary_points("kappa3_prime")
    k3_values = evaluate(k3_turning)["kappa3"]
    imax = int(np.argmax(k3_values))

    result = {
        "partA": {key: float(value) for key, value in a.items()},
        "partB": {
            "K10_scaled": float(np.sqrt(CI / D)),
            "front_peak_to_peak": f_max - f_min,
            "f": v["f"].tolist(),
            "kappa1": v["kappa1"].tolist(),
            "kappa2": v["kappa2"].tolist(),
            "kappa3": v["kappa3"].tolist(),
            "kappa3_max": float(k3_values[imax]),
            "s_kappa3_max": float(k3_turning[imax]),
        },
    }

    # Independent and numerical checks before any output is accepted.
    grid = np.linspace(-0.5, 0.5, 4096, endpoint=False)
    dense = evaluate(grid)
    np.testing.assert_allclose(
        dense["energy_from_k"], toughness_g(grid), rtol=0.0, atol=3e-14
    )
    np.testing.assert_allclose(
        dense["energy_from_k"], dense["energy_from_operator"], rtol=0.0, atol=3e-14
    )
    assert abs(float(np.mean(dense["f"]))) < 2e-15
    np.testing.assert_allclose(
        evaluate(-0.5)["f"], evaluate(0.5)["f"], rtol=0.0, atol=2e-15
    )
    assert f_max >= float(np.max(dense["f"])) - 2e-8
    assert f_min <= float(np.min(dense["f"])) + 2e-8
    assert k3_values[imax] >= float(np.max(dense["kappa3"])) - 2e-8
    hilbert_quadrature_check()

    OUT.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    print(f"Wrote {OUT.name}")
    print(f"front min/max = {f_min:.12f}, {f_max:.12f}")
    print(
        f"kappa3 maximum = {k3_values[imax]:.12f} "
        f"at s = {k3_turning[imax]:.12f}"
    )
    print("All equilibrium, energy, Hilbert, and extremum checks passed.")


if __name__ == "__main__":
    solve()
