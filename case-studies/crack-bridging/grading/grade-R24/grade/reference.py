"""Finite symmetric crack: analytic reduction checked with crack Green functions.

Run from any directory with Python, NumPy and SciPy. No paper or data files needed.
All internal units are SI. delta is the full face-to-face opening.
"""

import json
from math import acos, exp, expm1, log, pi, sinh, sqrt

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq


E = 130.0e9
NU = 0.28
EP = E / (1.0 - NU**2)
A = 2.4e-3  # original traction-free half-length
SIGMA_B = 160.0e6
DELTA_C = 9.0e-6
GAMMA_0 = 22.0


def analytic_solution(a=A, ep=EP, sigma_b=SIGMA_B, delta_c=DELTA_C,
                      gamma_0=GAMMA_0):
    """Solve for u = log(c/a) using the finite-crack opening and tip condition."""
    k0 = sqrt(ep * gamma_0)
    k = k0 / (sigma_b * sqrt(pi * a))
    target = ep * delta_c / (sigma_b * a)

    def residual(u):
        return 4.0 * k * sqrt(2.0 * sinh(u)) + 8.0 * u / pi - target

    upper = 1.0
    while residual(upper) < 0.0:
        upper *= 2.0
    u = brentq(residual, 0.0, upper, xtol=1.0e-14, rtol=1.0e-14)
    c = a * exp(u)
    sigma = sigma_b * (2.0 * acos(exp(-u)) / pi + k * exp(-u / 2.0))
    # Check that this endpoint precedes a maximum of the load along the branch.
    slope = 2.0 / (pi * sqrt(expm1(2.0 * u))) - 0.5 * k * exp(-u / 2.0)
    return {"sigma_Pa": sigma, "L_m": a * expm1(u), "c_m": c,
            "log_eta": u, "tip_K_Pa_sqrt_m": k0,
            "positive_load_branch_slope": slope}


def green_kernel(x, s):
    """Dimensionless exact crack-opening Green kernel; -1 < x,s < 1."""
    numerator = 1.0 - x * s + sqrt((1.0 - x*x) * max(0.0, 1.0 - s*s))
    return log(numerator / abs(x - s))


def opening_by_green(x, c, sigma, a, ep=EP, sigma_b=SIGMA_B):
    """Integrate the face-pressure kernel, independently of the opening formula.

    Face pressure is sigma on |s| < a and sigma - sigma_b on a < |s| < c.
    Symmetry combines the positive and negative halves of the crack.
    """
    r = x / c
    split = a / c

    def kernel_sum(t):
        return green_kernel(r, t) + green_kernel(r, -t)

    # x=a in the main calculation, so the logarithmic singularity is at the
    # integration boundary. quad handles this as an integrable singularity.
    edges = sorted(set([0.0, split, 1.0, abs(r)]))
    integral = 0.0
    for left, right in zip(edges[:-1], edges[1:]):
        if right <= left:
            continue
        pressure = sigma if (left + right) / 2.0 < split else sigma - sigma_b
        value = quad(kernel_sum, left, right, epsabs=2.0e-11,
                     epsrel=2.0e-11, limit=200)[0]
        integral += pressure * value
    return 4.0 * c * integral / (pi * ep)


def tip_by_weight(c, sigma, a, sigma_b=SIGMA_B):
    """Numerical symmetric mode-I weight-function integration.

    s=c*sin(theta) removes the square-root endpoint singularity.
    """
    theta_a = np.arcsin(a / c)
    inner = quad(lambda theta: sigma, 0.0, theta_a, epsabs=1.0e-7)[0]
    wing = quad(lambda theta: sigma - sigma_b, theta_a, pi / 2.0,
                epsabs=1.0e-7)[0]
    return 2.0 * sqrt(c / pi) * (inner + wing)


def green_solution():
    """Solve the same physical conditions without using the analytic opening."""
    k0 = sqrt(EP * GAMMA_0)

    def stress_for_tip(c):
        k_per_pa = tip_by_weight(c, 1.0, A, sigma_b=0.0)
        closing_k = -tip_by_weight(c, 0.0, A)
        return (k0 + closing_k) / k_per_pa

    def residual(eta):
        c = A * eta
        sigma = stress_for_tip(c)
        return opening_by_green(A, c, sigma, A) / DELTA_C - 1.0

    high = 2.0
    while residual(high) < 0.0:
        high *= 2.0
    eta = brentq(residual, 1.0 + 1.0e-8, high, xtol=1.0e-12, rtol=1.0e-12)
    c = A * eta
    sigma = stress_for_tip(c)
    return {"sigma_Pa": sigma, "L_m": c - A, "c_m": c,
            "opening_m": opening_by_green(A, c, sigma, A),
            "tip_K_Pa_sqrt_m": tip_by_weight(c, sigma, A)}


def main():
    first = analytic_solution()
    second = green_solution()
    relative_sigma = abs(first["sigma_Pa"] / second["sigma_Pa"] - 1.0)
    relative_length = abs(first["L_m"] / second["L_m"] - 1.0)
    assert relative_sigma < 2.0e-8
    assert relative_length < 2.0e-8
    assert first["positive_load_branch_slope"] > 0.0
    assert abs(second["opening_m"] / DELTA_C - 1.0) < 2.0e-8
    assert abs(second["tip_K_Pa_sqrt_m"] / sqrt(EP * GAMMA_0) - 1.0) < 2.0e-8

    # Independent normalization check: a uniformly pressurized crack has an
    # elliptical full opening, 4*p*sqrt(c^2-x^2)/E'.
    uniform_c = 6.0e-3
    uniform_x = 0.37 * uniform_c
    uniform_p = 1.0e6
    green_uniform = opening_by_green(uniform_x, uniform_c, uniform_p,
                                     uniform_c, sigma_b=0.0)
    exact_uniform = 4.0 * uniform_p * sqrt(uniform_c**2 - uniform_x**2) / EP
    uniform_error = abs(green_uniform / exact_uniform - 1.0)
    assert uniform_error < 2.0e-8

    # Zero intrinsic crack-front energy recovers the closed rectangular-law
    # finite-crack solution, with log(c/a) = pi*E'*delta_c/(8*sigma_b*a).
    zero_tip = analytic_solution(gamma_0=0.0)
    u0 = pi * EP * DELTA_C / (8.0 * SIGMA_B * A)
    special_sigma = SIGMA_B * 2.0 * acos(exp(-u0)) / pi
    special_length = A * expm1(u0)
    special_error = max(abs(zero_tip["sigma_Pa"] / special_sigma - 1.0),
                        abs(zero_tip["L_m"] / special_length - 1.0))
    assert special_error < 1.0e-11

    # Large-crack limit at zero intrinsic energy: Sigma -> sqrt(E'*W/(pi*a)).
    ell = EP * DELTA_C / SIGMA_B
    large_a = 1.0e5 * ell
    large = analytic_solution(a=large_a, gamma_0=0.0)
    large_sigma = sqrt(EP * SIGMA_B * DELTA_C / (pi * large_a))
    large_error = abs(large["sigma_Pa"] / large_sigma - 1.0)
    assert large_error < 1.0e-5

    lefm_sigma = sqrt(EP * (GAMMA_0 + SIGMA_B * DELTA_C) / (pi * A))
    result = {
        "inputs_SI": {"E": E, "nu": NU, "E_prime": EP, "a": A,
                      "sigma_b": SIGMA_B, "delta_c": DELTA_C, "Gamma_0": GAMMA_0},
        "answer": {"Sigma_MPa": first["sigma_Pa"] / 1.0e6,
                   "L_mm_per_end": first["L_m"] * 1.0e3},
        "analytic": first,
        "independent_Green_integrals": second,
        "checks": {"relative_sigma_difference": relative_sigma,
                   "relative_L_difference": relative_length,
                   "uniform_pressure_Green_relative_error": uniform_error,
                   "zero_tip_special_solution_relative_error": special_error,
                   "large_crack_LEFM_relative_error": large_error},
        "diagnostics": {"material_length_mm": ell * 1.0e3,
                        "a_over_material_length": A / ell,
                        "L_over_a": first["L_m"] / A,
                        "unbridged_LEFM_prediction_MPa": lefm_sigma / 1.0e6,
                        "LEFM_overprediction_percent":
                            100.0 * (lefm_sigma / first["sigma_Pa"] - 1.0)}
    }
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
