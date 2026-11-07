"""Organizer reproduction of Hui and Long (2026); never distribute to agents.

This implements published equations, not an independent reference derivation.
Coordinates: X = pi*x/h, Y = pi*y/h. Stress T = (tau2+i*tau1)/tau_c;
strain E = G_infinity*(gamma2+i*gamma1)/tau_c.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.integrate import quad, quad_vec, solve_ivp
from scipy.optimize import brentq
from scipy.special import digamma, roots_legendre
from scipy.sparse import eye


def log_kernel(xi, epsilon):
    """Logarithmic kernel in paper (51)-(52), evaluated without cancellation."""
    one_minus_alpha2 = -np.expm1(-xi) / -np.expm1(-epsilon)
    alpha = np.sqrt(np.maximum(0.0, 1.0 - one_minus_alpha2))
    return 2.0 * np.log1p(alpha) - np.log(one_minus_alpha2)


def separation(epsilon, speed, modulus_ratio=0.1):
    """Normalized critical separation from paper (51)."""
    if epsilon <= 0:
        return 0.0

    def integrand(t):
        if t == 0:
            return 0.0
        xi = epsilon * t * t
        creep = 1.0 - (1.0 - modulus_ratio) * np.exp(-xi / speed)
        return 2.0 * epsilon * t * creep * log_kernel(xi, epsilon)

    return 2.0 / np.pi**2 * quad(
        integrand, 0.0, 1.0, epsabs=2e-12, epsrel=2e-11, limit=250
    )[0]


def envelope_kernel(speed):
    """Paper (59), summed using a digamma identity derived by the organizer."""
    return 2.0 * speed / np.pi**2 * (
        digamma(0.5 + 1.0 / speed) - digamma(0.5)
    )


def envelope_kernel_quadrature(speed):
    """Direct evaluation of paper (58); numerical alternative to the series."""
    def integrand(t):
        if t == 0:
            return 0.0
        xi = t * t
        q = np.exp(-xi / 2.0)
        log_ratio = np.log1p(q) - np.log(-np.expm1(-xi / 2.0))
        return 2.0 * t * np.exp(-xi / speed) * log_ratio

    return 2.0 / np.pi**2 * quad(
        integrand, 0.0, np.inf, epsabs=2e-12, epsrel=2e-11
    )[0]


def critical_speed(modulus_ratio=0.1, delta_bar=0.2):
    if not (0 < modulus_ratio <= 1 and 0 < delta_bar < 1):
        raise ValueError("Need 0 < modulus_ratio <= 1 and 0 < delta_bar < 1")
    if delta_bar <= modulus_ratio:
        return None
    f = lambda v: 1 - (1 - modulus_ratio) * envelope_kernel(v) - delta_bar
    upper = 1.0
    while f(upper) > 0:
        upper *= 2
    return brentq(f, 1e-10, upper, xtol=2e-11, rtol=2e-12)


def solve_case(speed, modulus_ratio=0.1, delta_bar=0.2):
    if speed <= 0:
        raise ValueError("Use a strictly positive speed; the zero-speed limit is separate")
    if not (0 < modulus_ratio <= 1 and 0 < delta_bar < 1):
        raise ValueError("Need 0 < modulus_ratio <= 1 and 0 < delta_bar < 1")
    maximum_separation = 1.0 - (1.0 - modulus_ratio) * envelope_kernel(speed)
    if maximum_separation <= delta_bar:
        return {"V": speed, "admissible": False}
    upper = 1.0
    while separation(upper, speed, modulus_ratio) < delta_bar:
        upper *= 2
        if upper > 1024:
            raise ArithmeticError("Root too close to the infinite-zone endpoint")
    epsilon = brentq(
        lambda e: separation(e, speed, modulus_ratio) - delta_bar,
        1e-12, upper, xtol=3e-12, rtol=3e-12,
    )
    a = 2.0 / np.pi * np.arccos(np.exp(-epsilon / 2.0))
    return {
        "V": speed, "admissible": True, "epsilon": float(epsilon),
        "a": float(a), "R": float(a * a / delta_bar),
    }


def stress(X, Y, epsilon):
    """Complex stress in paper (12), upper half strip only.

    The omitted far-left stress is below about 1e-13 in normalized units.
    Exact singular points on the crack line are excluded from field sampling.
    """
    X, Y = np.broadcast_arrays(np.asarray(X, float), np.asarray(Y, float))
    z = X.astype(complex) + 1j * Y
    out = np.zeros(z.shape, complex)
    active = X >= -epsilon - 60.0
    if epsilon > 600:
        raise ArithmeticError("Stress evaluation needs rescaling this close to the endpoint")
    zz = z[active]
    argument = np.sqrt((1.0 - np.exp(-zz)) / np.expm1(epsilon))
    out[active] = 1.0 - 2.0 / np.pi * np.arctan(argument)
    return out


def strain_convolution(X, Y, case, modulus_ratio=0.1):
    """Paper (38) integrated by parts to use stress, not its derivative.

    E = lambda*T + (1-lambda)*integral_0^infinity exp(-u) T(X+V*u,Y) du.
    """
    X, Y = np.broadcast_arrays(np.asarray(X, float), np.asarray(Y, float))
    speed, epsilon = case["V"], case["epsilon"]
    if modulus_ratio == 1:
        return stress(X, Y, epsilon)
    breaks = sorted(set(
        v for x in X.ravel() for v in ((-epsilon - x) / speed, -x / speed)
        if v > 0
    ))
    # Beyond u=40 the exponential weight contributes < 5e-18 times bounded
    # sampled stress. Splitting prevents missed sharp changes at large V.
    breaks = [0.0] + [v for v in breaks if v < 40] + [40.0]
    mean = np.zeros(X.shape, complex)
    for lo, hi in zip(breaks[:-1], breaks[1:]):
        value, _ = quad_vec(
            lambda u: np.exp(-u) * stress(X + speed * u, Y, epsilon),
            lo, hi, epsabs=3e-12, epsrel=3e-11, limit=500,
        )
        mean += value
    return modulus_ratio * stress(X, Y, epsilon) + (1 - modulus_ratio) * mean


def history_solution(Y, case, modulus_ratio=0.1, tail=18.0, rtol=2e-9):
    """Alternative strain evaluation by the local SLS history equation (82).

    Filter A obeys dA/dX=(A-T)/V, integrated from relaxed far ahead to behind.
    This shares the published stress field with the convolution route.
    """
    Y = np.atleast_1d(np.asarray(Y, float))
    epsilon, speed = case["epsilon"], case["V"]
    right, left = 32.0, -epsilon - max(32.0, tail * speed)
    initial = stress(right, Y, epsilon)
    sol = solve_ivp(
        lambda x, a: (a - stress(x, Y, epsilon)) / speed,
        (right, left), initial, method="BDF", dense_output=True,
        jac=eye(len(Y), format="csc") / speed,
        rtol=rtol, atol=rtol * 0.01,
    )
    if not sol.success:
        raise ArithmeticError(sol.message)
    return sol, left, right


def bulk_dissipation(case, modulus_ratio=0.1, delta_bar=0.2, ny=32,
                     tail=18.0, rtol=2e-9):
    """Spatial integral of paper (84)-(86), without using R-1 as input.

    Transformed Gaussian quadrature Y=pi*t^2 resolves the crack-line layer.
    SLS history is solved along each ordinate, then power density is integrated
    over X with adaptive quadrature. The lower half is supplied by symmetry.
    """
    if modulus_ratio == 1:
        return {"value": 0.0, "ny": ny, "quadrature_error": 0.0}
    nodes, weights = roots_legendre(ny)
    t, wt = (nodes + 1.0) / 2.0, weights / 2.0
    Y, wy = np.pi * t**2, 2.0 * np.pi * t * wt
    sol, left, right = history_solution(Y, case, modulus_ratio, tail, rtol)
    epsilon, speed = case["epsilon"], case["V"]

    def density_integral(x):
        T = stress(x, Y, epsilon)
        E = modulus_ratio * T + (1 - modulus_ratio) * sol.sol(x)
        power = np.abs(T - E)**2 / (1 - modulus_ratio)
        return float(np.dot(wy, power))

    value, error = 0.0, 0.0
    for lo, hi in zip([left, -epsilon, 0.0], [-epsilon, 0.0, right]):
        val, err = quad(
            density_integral, lo, hi, epsabs=2e-10 * max(speed, 0.001),
            epsrel=2e-7, limit=700,
        )
        value += val
        error += err
    factor = 2.0 / (np.pi * delta_bar * speed)
    return {
        "value": float(factor * value), "ny": ny,
        "X_min": float(left), "X_max": float(right),
        "tail_times": tail, "ode_rtol": rtol, "ode_evaluations": sol.nfev,
        "quadrature_error": float(factor * error),
    }


def field_rows(case, sample_X, sample_y_over_h, modulus_ratio=0.1):
    X, yh = np.meshgrid(sample_X, sample_y_over_h, indexing="ij")
    E = strain_convolution(X.ravel(), np.pi * yh.ravel(), case, modulus_ratio)
    T = stress(X.ravel(), np.pi * yh.ravel(), case["epsilon"])
    return [
        {"X": float(x), "y_over_h": float(y), "T1": float(t.imag),
         "T2": float(t.real), "E1": float(e.imag), "E2": float(e.real)}
        for x, y, t, e in zip(X.ravel(), yh.ravel(), T, E)
    ]


def build_reference(inputs, energy_ny=48):
    lam, d = inputs["lambda"], inputs["delta_bar"]
    cases = []
    for speed in inputs["speeds"]:
        case = solve_case(speed, lam, d)
        if case["admissible"]:
            case["fields"] = field_rows(
                case, inputs["sample_X"], inputs["sample_y_over_h"], lam
            )
            energy = bulk_dissipation(case, lam, d, ny=energy_ny)
            case["bulk_dissipation_ratio"] = energy["value"]
            case["energy_integration"] = energy
        cases.append(case)
    return {
        "basis": "published-solution organizer reproduction; not independent review",
        "lambda": lam, "delta_bar": d,
        "critical_speed": critical_speed(lam, d), "cases": cases,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("input", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--energy-ny", type=int, default=48)
    args = parser.parse_args()
    result = build_reference(json.loads(args.input.read_text()), args.energy_ny)
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    print(json.dumps({"critical_speed": result["critical_speed"], "cases": [
        {k: c[k] for k in ("V", "epsilon", "a", "R", "bulk_dissipation_ratio")}
        for c in result["cases"]
    ]}, indent=2))
