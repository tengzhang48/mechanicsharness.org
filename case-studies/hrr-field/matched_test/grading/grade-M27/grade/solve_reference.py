"""Antiplane crack: angular IVP versus an independent phase solution.

Run: python V1/solve_reference.py
Requires only numpy and scipy. Units: mm, MPa, N.
"""

import json
import numpy as np
from scipy.integrate import solve_ivp, quad
from scipy.optimize import brentq

N = 3.5
TAU0 = 80.0
GAMMA0 = 0.002
J = 0.12
R1, R2 = 0.02, 0.08
GAMMA_STAR = 0.01


def angular_ivp(n):
    lam = 1.0 / (n + 1.0)
    exponent = 1.0 / n - 1.0
    b = lam * (1.0 - lam)

    def rhs(theta, y):
        f, g = y
        q2 = (lam * f) ** 2 + g**2
        gprime = -f * (b * q2 + exponent * lam**2 * g**2) / (
            q2 + exponent * g**2
        )
        return g, gprime

    sol = solve_ivp(rhs, (0.0, np.pi), (0.0, 1.0), method="DOP853",
                    rtol=2e-12, atol=2e-13, max_step=0.025,
                    dense_output=True)
    if not sol.success:
        raise RuntimeError(sol.message)
    return sol.sol


def phase_solution(n, theta):
    """Independent closed parametrization, normalized by f'(0)=1."""
    lam = 1.0 / (n + 1.0)
    if theta <= 0:
        psi = 0.0
    elif theta >= np.pi:
        psi = np.pi / 2
    else:
        psi = brentq(lambda v: v + np.arctan2(n * np.sin(v), np.cos(v))
                     - theta, 0.0, np.pi / 2, xtol=5e-15)
    d = np.cos(psi)**2 + n**2 * np.sin(psi)**2
    q = d ** (-n / (2 * (n + 1)))
    return np.array([q * np.sin(psi) / lam, q * np.cos(psi)])


def flux_integrand(n, theta, angular):
    lam = 1.0 / (n + 1)
    p = (n + 1.0) / n
    f, g = angular(theta)
    q = np.hypot(lam * f, g)
    return ((n / (n + 1)) * q**p * np.cos(theta)
            - q**(p - 2) * lam * f *
            (lam * f * np.cos(theta) - g * np.sin(theta)))


def compute(n, tau0, gamma0, j, angular, integral):
    lam = 1.0 / (n + 1)
    p = (n + 1.0) / n
    kappa = tau0 * gamma0 ** (-1.0 / n)
    a = (j / (kappa * integral)) ** (1.0 / p)

    def gradient(r, theta):
        f, g = angular(theta)
        return a * r**(lam - 1) * np.array([lam * f, g])

    def shear(r, theta):
        grad = gradient(r, theta)
        return kappa * np.linalg.norm(grad)**(1.0 / n - 1) * grad

    def contour(theta):
        f, g = angular(theta)
        q = np.hypot(lam * f, g)
        return (a * q / GAMMA_STAR) ** (1 / (1 - lam))

    values = {
        "slip_0.02_mm": float(2 * a * R1**lam * angular(np.pi)[0]),
        "slip_0.08_mm": float(2 * a * R2**lam * angular(np.pi)[0]),
        "tau_rz_0.02mm_60deg_MPa": float(shear(R1, np.pi / 3)[0]),
        "tau_thetaz_0.02mm_60deg_MPa": float(shear(R1, np.pi / 3)[1]),
        "gamma_contour_0deg_mm": float(contour(0.0)),
        "gamma_contour_90deg_mm": float(contour(np.pi / 2)),
        "gamma_contour_180deg_mm": float(contour(np.pi)),
    }

    def circle_j(r):
        def integrand(theta):
            grad = gradient(r, theta)
            tau = shear(r, theta)
            w = (n / (n + 1)) * np.dot(tau, grad)
            ux = grad[0] * np.cos(theta) - grad[1] * np.sin(theta)
            return r * (w * np.cos(theta) - tau[0] * ux)
        return 2 * quad(integrand, 0, np.pi, epsabs=2e-12,
                        epsrel=2e-12)[0]

    return values, {"lambda": lam, "a": a, "I": integral,
                    "face_fprime": float(angular(np.pi)[1]),
                    "J_by_circle": {str(r): circle_j(r)
                                    for r in (0.005, 0.02, 0.08)},
                    "gamma_max_at_r1": float(np.linalg.norm(gradient(R1, 0)))}


def main():
    ode = angular_ivp(N)
    i_num = 2 * quad(lambda t: flux_integrand(N, t, ode), 0, np.pi,
                     epsabs=2e-12, epsrel=2e-12)[0]
    phase = lambda t: phase_solution(N, t)
    i_exact = np.pi * (N + 1) / (2 * N)
    va, da = compute(N, TAU0, GAMMA0, J, ode, i_num)
    vb, db = compute(N, TAU0, GAMMA0, J, phase, i_exact)
    errors = {k: abs(va[k] - vb[k]) / abs(vb[k]) for k in va}
    grid = np.linspace(0, np.pi, 2001)
    angular_error = max(float(np.max(np.abs(ode(t) - phase(t)))) for t in grid)

    # Known special solution: n=1 is linear antiplane elasticity.
    # For f'(0)=1, f=2 sin(theta/2), and I=pi.
    elastic = angular_ivp(1.0)
    elastic_f_error = max(abs(elastic(t)[0] - 2 * np.sin(t / 2)) for t in grid)
    elastic_i = 2 * quad(lambda t: flux_integrand(1, t, elastic), 0, np.pi)[0]
    # J=K_III^2/(2G) independently gives the displacement coefficient.
    gmod, k3 = 26000.0, 35.0
    j_el = k3**2 / (2 * gmod)
    a_el = np.sqrt(j_el / (gmod * elastic_i))
    r_el = 0.03
    slip_ode = 2 * a_el * np.sqrt(r_el) * elastic(np.pi)[0]
    slip_known = 2 * k3 / gmod * np.sqrt(2 * r_el / np.pi)

    out = {
        "route_A_angular_ODE": va, "route_B_phase_formula": vb,
        "route_A_diagnostics": da, "route_B_diagnostics": db,
        "relative_errors_A_vs_B": errors,
        "max_angular_absolute_error": angular_error,
        "linear_limit": {
            "max_f_absolute_error": float(elastic_f_error),
            "I_numerical": float(elastic_i), "I_exact": float(np.pi),
            "slip_ODE_mm": float(slip_ode),
            "slip_known_modeIII_mm": float(slip_known),
            "relative_slip_error": float(abs(slip_ode / slip_known - 1))},
    }
    assert max(errors.values()) < 1e-9
    assert angular_error < 1e-9
    assert abs(elastic_i / np.pi - 1) < 1e-10
    assert abs(slip_ode / slip_known - 1) < 1e-10
    assert abs(da["face_fprime"]) < 1e-10
    assert all(abs(v / J - 1) < 1e-10 for v in da["J_by_circle"].values())
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
