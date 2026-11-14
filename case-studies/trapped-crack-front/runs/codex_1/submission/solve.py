#!/usr/bin/env python3
"""Solve PUBLIC_TASK.md at first order and write results.json.

Only NumPy and SciPy are required. Run from this folder with:
    python3 solve.py
All numerical checks run before the result file is written.
"""

import json
import math
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.optimize import minimize_scalar


NU = 0.3
R2 = 0.6
R3 = 0.8
Q = 0.5

# Coefficients in the in-plane crack-front weight-function relation.
P2 = (2.0 - 3.0 * NU) / (2.0 - NU)
P3 = (2.0 + NU) / (2.0 - NU)
B2 = 2.0 / (2.0 - NU)
B3 = 2.0 * (1.0 - NU) / (2.0 - NU)

# G0 = (1-nu^2) K10^2 D / E; the energy perturbation is -L |d/dz| a.
D = 1.0 + R2**2 + R3**2 / (1.0 - NU)
H = 1.0 + P2 * R2**2 + P3 * R3**2 / (1.0 - NU)
L = H / D


def toughness(s):
    """The zero-mean fractional toughness perturbation g(s)."""
    theta = 2.0 * np.pi * np.asarray(s)
    denominator = 1.0 - 2.0 * Q * np.cos(theta) + Q**2
    return (1.0 - Q**2) / denominator - 1.0


def profiles(s):
    """Return f, kappa1, kappa2, and kappa3 at scalar/array positions s."""
    theta = 2.0 * np.pi * np.asarray(s)
    denominator = 1.0 - 2.0 * Q * np.cos(theta) + Q**2
    f = np.log(denominator) / (2.0 * np.pi * L)
    f_prime = 2.0 * Q * np.sin(theta) / (L * denominator)
    line_operator_f = toughness(s) / (2.0 * L)
    kappa1 = line_operator_f
    kappa2 = P2 * R2 * line_operator_f - B2 * R3 * f_prime
    kappa3 = P3 * R3 * line_operator_f + B3 * R2 * f_prime
    return f, kappa1, kappa2, kappa3


def mode3_extrema():
    """Find both stationary points analytically, then select the maximum.

    kappa3(theta) = (a + b sin(theta))/(h - d cos(theta)) - c.
    Its stationary equation is b h cos(theta) - a d sin(theta) = b d.
    """
    c = P3 * R3 / (2.0 * L)
    a = c * (1.0 - Q**2)
    b = 2.0 * Q * B3 * R2 / L
    h = 1.0 + Q**2
    d = 2.0 * Q
    phase = math.atan2(a * d, b * h)
    angle = math.acos(b * d / math.hypot(b * h, a * d))
    candidates = []
    for theta in (-phase + angle, -phase - angle):
        s = (theta / (2.0 * math.pi) + 0.5) % 1.0 - 0.5
        candidates.append((float(profiles(s)[3]), s))
    return max(candidates), candidates


def make_results():
    part_a = {
        "C1": -math.pi,
        "S1": 0.0,
        "C2": -math.pi * P2 * R2,
        "S2": 2.0 * math.pi * B2 * R3,
        "C3": -math.pi * P3 * R3,
        "S3": -2.0 * math.pi * B3 * R2,
        "CG": -2.0 * math.pi * L,
        "SG": 0.0,
    }
    samples = profiles(np.arange(16, dtype=float) / 16.0)
    (kappa3_max, s_kappa3_max), _ = mode3_extrema()
    part_b = {
        "K10_scaled": 1.0 / math.sqrt(D),
        "front_peak_to_peak": math.log((1.0 + Q) / (1.0 - Q))
        / (math.pi * L),
        "f": samples[0].tolist(),
        "kappa1": samples[1].tolist(),
        "kappa2": samples[2].tolist(),
        "kappa3": samples[3].tolist(),
        "kappa3_max": kappa3_max,
        "s_kappa3_max": s_kappa3_max,
    }
    return {"partA": part_a, "partB": part_b}


def line_operator_quadrature(s):
    """Evaluate J[f] by its periodic, symmetrized real-space integral.

    J[f](s) = (pi/2) integral_0^(1/2)
              (f(s+u) + f(s-u) - 2 f(s)) / sin(pi*u)^2 du.

    For the logarithmic profile the numerator is evaluated using log1p,
    avoiding subtraction of almost equal numbers near u=0.
    """
    h = 1.0 + Q**2
    d = 2.0 * Q
    cosine = math.cos(2.0 * math.pi * s)
    denominator = h - d * cosine

    def integrand(u):
        sin2 = math.sin(math.pi * u) ** 2
        cos2 = math.cos(math.pi * u) ** 2
        if sin2 == 0.0:
            return d * (h * cosine - d) / (L * denominator**2)
        relative_change = (
            4.0 * d * sin2 * (h * cosine - d * cos2) / denominator**2
        )
        return math.log1p(relative_change) / (4.0 * L * sin2)

    return quad(integrand, 0.0, 0.5, epsabs=2e-12, epsrel=2e-12)[0]


def run_checks(results):
    """Check closed forms, the singular integral, equilibrium, and extrema."""
    checks = {}

    def check(name, error, tolerance):
        error = float(np.max(np.abs(error)))
        if not math.isfinite(error) or error > tolerance:
            raise AssertionError(f"{name}: error {error}, tolerance {tolerance}")
        checks[name] = error

    part_a = results["partA"]
    energy_cos = 2.0 * (
        part_a["C1"]
        + R2 * part_a["C2"]
        + R3 * part_a["C3"] / (1.0 - NU)
    ) / D
    energy_sin = 2.0 * (
        part_a["S1"]
        + R2 * part_a["S2"]
        + R3 * part_a["S3"] / (1.0 - NU)
    ) / D
    check("Part A energy cosine", energy_cos - part_a["CG"], 1e-13)
    check("Part A energy sine cancellation", energy_sin, 1e-13)
    check("Loading normalization", D * results["partB"]["K10_scaled"] ** 2 - 1, 1e-14)

    s = np.arange(4096, dtype=float) / 4096.0
    closed = np.asarray(profiles(s))
    check("Mean front displacement", np.mean(closed[0]), 1e-14)
    check("Mean SIF perturbations", np.mean(closed[1:], axis=1), 1e-14)
    check("Periodicity", np.asarray(profiles(s + 1.0)) - closed, 2e-14)
    check("Front reflection symmetry", np.asarray(profiles(-s))[0] - closed[0], 1e-14)

    # Sum the toughness Fourier series to 96 harmonics, solve separately for
    # each front Fourier coefficient, and apply the SIF operator mode by mode.
    n = np.arange(1, 97, dtype=float)
    angles = 2.0 * np.pi * np.outer(s, n)
    front_cosine_coeff = 2.0 * Q**n / (n * part_a["CG"])
    fourier_f = np.cos(angles) @ front_cosine_coeff
    fourier_j = np.cos(angles) @ (-np.pi * n * front_cosine_coeff)
    fourier_f_prime = np.sin(angles) @ (-2.0 * np.pi * n * front_cosine_coeff)
    fourier_profiles = np.asarray(
        [
            fourier_f,
            fourier_j,
            P2 * R2 * fourier_j - B2 * R3 * fourier_f_prime,
            P3 * R3 * fourier_j + B3 * R2 * fourier_f_prime,
        ]
    )
    check("96-harmonic solution versus closed forms", fourier_profiles - closed, 2e-14)

    quadrature_s = np.arange(16, dtype=float) / 16.0
    quadrature_values = np.array([line_operator_quadrature(x) for x in quadrature_s])
    check(
        "Real-space singular integral versus Fourier operator",
        quadrature_values - profiles(quadrature_s)[1],
        2e-11,
    )

    energy_profile = 2.0 * (
        closed[1] + R2 * closed[2] + R3 * closed[3] / (1.0 - NU)
    ) / D
    check("Pointwise first-order equilibrium", energy_profile - toughness(s), 2e-14)
    check(
        "Front peak-to-peak",
        np.ptp(closed[0]) - results["partB"]["front_peak_to_peak"],
        1e-14,
    )

    # Locate the global maximum on a fine full-period grid, then optimize in
    # its neighborhood; compare with the two analytical stationary points.
    grid_s = np.arange(65536, dtype=float) / 65536.0 - 0.5
    grid_k3 = profiles(grid_s)[3]
    grid_best = grid_s[np.argmax(grid_k3)]
    radius = 2.0 / len(grid_s)
    optimum = minimize_scalar(
        lambda x: -float(profiles(x)[3]),
        bounds=(grid_best - radius, grid_best + radius),
        method="bounded",
        options={"xatol": 1e-14},
    )
    if not optimum.success:
        raise AssertionError(optimum.message)
    maximum, candidates = mode3_extrema()
    check("Mode III maximum value", -optimum.fun - maximum[0], 1e-12)
    check("Mode III maximum location", optimum.x - maximum[1], 1e-8)
    if not (-0.5 <= maximum[1] < 0.5 and len(candidates) == 2):
        raise AssertionError("Unexpected stationary points")
    if np.max(grid_k3) > maximum[0] + 1e-12:
        raise AssertionError("A grid value exceeds the analytical global maximum")

    # Confirm that the serialized deliverable has exactly the requested keys,
    # finite numbers, and sixteen entries in each sampled profile.
    decoded = json.loads(json.dumps(results, allow_nan=False))
    if set(decoded) != {"partA", "partB"}:
        raise AssertionError("Unexpected top-level keys")
    if set(decoded["partA"]) != {"C1", "S1", "C2", "S2", "C3", "S3", "CG", "SG"}:
        raise AssertionError("Unexpected Part A keys")
    expected_b = {
        "K10_scaled", "front_peak_to_peak", "f", "kappa1", "kappa2", "kappa3",
        "kappa3_max", "s_kappa3_max",
    }
    if set(decoded["partB"]) != expected_b:
        raise AssertionError("Unexpected Part B keys")
    for name in ("f", "kappa1", "kappa2", "kappa3"):
        if len(decoded["partB"][name]) != 16:
            raise AssertionError(f"Incorrect number of samples for {name}")

    for name, error in checks.items():
        print(f"PASS {name}: maximum absolute error {error:.3e}")
    print("PASS result schema and finite JSON numbers")


def main():
    results = make_results()
    run_checks(results)
    output = Path(__file__).resolve().parent / "results.json"
    output.write_text(json.dumps(results, indent=2, allow_nan=False) + "\n")
    print(f"Wrote {output.name}")
    print(f"kappa3_max = {results['partB']['kappa3_max']:.15g}")
    print(f"s_kappa3_max = {results['partB']['s_kappa3_max']:.15g}")


if __name__ == "__main__":
    main()
