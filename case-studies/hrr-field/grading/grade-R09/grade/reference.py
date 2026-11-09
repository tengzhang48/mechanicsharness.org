"""Standalone notch reference: mixed stress/strain BVP, explicit IVP check.

Run: python T2/reference.py
No imports from T1 or scratch. Requires only NumPy and SciPy.
"""
from pathlib import Path
import json

import numpy as np
from scipy.integrate import solve_bvp, solve_ivp
from scipy.optimize import brentq, root


def invert(g, v, n, s):
    """Solve g = D [3 D^2/4 + 3(s-1)^2 v^2]^((n-1)/2)."""
    if n == 1:
        return g
    target = np.abs(g)
    a, b, m = 0.75, 3*(s - 1)**2, (n - 1)/2
    # An upper bound, followed by Newton iteration of a monotone scalar law.
    magnitude = (target/a**m)**(1/n)
    for _ in range(45):
        q2 = np.maximum(a*magnitude**2 + b*v**2, 1e-26)
        factor = q2**m
        residual = magnitude*factor - target
        derivative = factor*(1 + (n - 1)*a*magnitude**2/q2)
        new = np.maximum(magnitude - residual/derivative, 0)
        if np.max(np.abs(new - magnitude)) < 2e-14:
            magnitude = new
            break
        magnitude = new
    return np.sign(g)*magnitude


def mixed_ode(theta, y, n, s):
    f, v, g, gp = y
    h, c, k = s - 1, s*(2 - s), n*(s - 2)
    d = invert(g, v, n, s)
    q2 = np.maximum(0.75*d**2 + 3*h**2*v**2, 1e-26)
    factor = q2**((n - 1)/2)
    vp = d - c*f
    dp = (gp/factor - (n - 1)*d*3*h**2*v*vp/q2)/(1 + (n - 1)*0.75*d**2/q2)
    factor_v_prime = factor*(vp + (n - 1)*v*(0.75*d*dp + 3*h**2*v*vp)/q2)
    gpp = k*(k + 2)*g - 4*h*(k + 1)*factor_v_prime
    return np.array([v, vp, gp, gpp])


def mixed_bvp(n, beta, tol=1e-9):
    theta = np.linspace(0, beta, 301)
    t = theta/beta
    f, v = 1 - 2*t**2 + t**4, (-4*t + 4*t**3)/beta
    w, z = (-4 + 12*t**2)/beta**2, 24*t/beta**3
    s = 2 - 1/(n + 1)
    d, dp = w + s*(2 - s)*f, z + s*(2 - s)*v
    q2 = np.maximum(0.75*d**2 + 3*(s - 1)**2*v**2, 1e-26)
    factor = q2**((n - 1)/2)
    g = factor*d
    gp = factor*(dp + (n - 1)*d*(0.75*d*dp + 3*(s - 1)**2*v*w)/q2)
    def bc(left, right, par):
        return np.array([left[0] - 1, left[1], left[3], right[0], right[1]])
    sol = solve_bvp(lambda x, y, par: mixed_ode(x, y, n, par[0]), bc,
                    theta, np.array([f, v, g, gp]), p=[s], tol=tol, max_nodes=30000)
    if not sol.success:
        raise RuntimeError(sol.message)
    s = sol.p[0]
    f0, v0, g0, gp0 = sol.y[:, 0]
    curvature = float(invert(np.array([g0]), np.array([v0]), n, s)[0] - s*(2 - s)*f0)
    ratio = (s*f0 + curvature)/(s*(s - 1)*f0)
    return sol, dict(s=s, p=2 - s, R=ratio, curvature=curvature,
                     maximum_BVP_residual=float(np.max(sol.rms_residuals)),
                     endpoint_residual=float(max(abs(sol.y[0, -1]), abs(sol.y[1, -1]))))


def explicit_fourth(theta, y, n, s):
    """Independent variables f,f',f'',f''' rather than mixed constitutive inversion."""
    f, v, w, z = y
    h, c, k = s - 1, s*(2 - s), n*(s - 2)
    d, dp = w + c*f, z + c*v
    q2 = max(0.75*d*d + 3*h*h*v*v, 1e-26)
    q2p = 2*(0.75*d*dp + 3*h*h*v*w)
    q2pp_without_f4 = 2*(0.75*(dp*dp + d*c*w) + 3*h*h*(w*w + v*z))
    m = (n - 1)/2
    factor_prime_over_factor = m*q2p/q2
    factor_second_rest = m*q2pp_without_f4/q2 + m*(m - 1)*(q2p/q2)**2
    rest = c*w + 2*factor_prime_over_factor*dp + d*factor_second_rest
    f4 = (k*(k + 2)*d - 4*h*(k + 1)*(w + factor_prime_over_factor*v) - rest)
    f4 /= 1 + (n - 1)*0.75*d*d/q2
    return [v, w, z, f4]


def shooting(n, beta, start):
    def integrate(par):
        s, curvature = par
        sol = solve_ivp(lambda x, y: explicit_fourth(x, y, n, s), (0, beta),
                         [1, 0, curvature, 0], method="DOP853", rtol=2e-12,
                         atol=2e-13, max_step=beta/180, dense_output=True)
        if not sol.success:
            raise RuntimeError(sol.message)
        return sol
    fit = root(lambda par: integrate(par).y[:2, -1], start, tol=1e-10)
    sol = integrate(fit.x)
    if max(abs(sol.y[0, -1]), abs(sol.y[1, -1])) > 1e-8:
        raise RuntimeError(fit.message)
    s, curvature = fit.x
    return sol, dict(s=s, p=2 - s, R=(s + curvature)/(s*(s - 1)), curvature=curvature,
                     endpoint_residual=float(max(abs(sol.y[0, -1]), abs(sol.y[1, -1]))))


def linear_exact(beta):
    # Symmetric traction-free elastic wedge determinant; lambda=s-1.
    lam = brentq(lambda v: np.sin(2*beta*v) + v*np.sin(2*beta), 0.1, 0.99,
                 xtol=1e-14)
    s = 1 + lam
    b_over_a = -np.cos(s*beta)/np.cos((s - 2)*beta)
    aa, bb = 1/(1 + b_over_a), b_over_a/(1 + b_over_a)
    curvature = -aa*s*s - bb*(s - 2)**2
    return dict(p=1 - lam, R=(s + curvature)/(s*(s - 1)))


def main():
    n, beta = 4, 2*np.pi/3
    sol, primary = mixed_bvp(n, beta)
    ivp, secondary = shooting(n, beta, [primary["s"], primary["curvature"]])
    differences = {k: abs(primary[k] - secondary[k]) for k in ("p", "R")}
    theta = np.linspace(0, beta, 5001)
    # The methods also agree pointwise on f, rather than only on the reported scalars.
    shape_error = float(np.max(np.abs(sol.sol(theta)[0] - ivp.sol(theta)[0])))
    linear_sol, linear_computed = mixed_bvp(1, beta)
    exact = linear_exact(beta)
    linear_errors = {k: abs(linear_computed[k] - exact[k]) for k in ("p", "R")}
    crack_sol, crack = mixed_bvp(n, np.pi)
    crack_exact_p = 1/(n + 1)
    historical_sol, historical = mixed_bvp(3, 3*np.pi/4)

    assert max(differences.values()) < 1e-7
    assert shape_error < 1e-7
    assert max(linear_errors.values()) < 1e-7
    assert abs(crack["p"] - crack_exact_p) < 1e-8
    assert abs(historical["p"] - 0.225) < 5e-5
    assert 0 < primary["p"] < 2/(n + 1)  # Finite tip energy.
    assert np.min(sol.sol(theta)[0, :-1]) > 0  # Leading opening branch.
    result = dict(primary=primary, explicit_shooting=secondary,
                  absolute_differences=differences, maximum_shape_difference=shape_error,
                  linear_wedge_check=dict(computed=linear_computed, analytic=exact,
                                           absolute_errors=linear_errors),
                  zero_opening_angle_check=dict(computed=crack, exact_p=crack_exact_p),
                  rounded_90_degree_n3_check=dict(computed=historical, reported_p=0.225),
                  radial_energy_power=2 - (n + 1)*primary["p"])
    serialized = json.dumps(result, indent=2)
    Path(__file__).with_name("reference_results.json").write_text(serialized + "\n")
    print(serialized)


if __name__ == "__main__":
    main()
