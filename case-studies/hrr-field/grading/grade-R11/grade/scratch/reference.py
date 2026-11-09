"""Standalone reference: nonlinear plane-strain crack field (NumPy/SciPy).

Run: python T1/reference.py
The fourth-order angular equation is solved by collocation and by IVP shooting.
J is also evaluated using Cartesian finite differences of displacements, and
the n=1 incompressible elastic solution provides an analytic normalization check.
"""
from pathlib import Path
import json

import numpy as np
from scipy.integrate import simpson, solve_bvp, solve_ivp
from scipy.optimize import root


def ode(theta, y, n, s):
    f, v, w, z = y
    h, c, k = s - 1, s * (2 - s), n * (s - 2)
    d, dp = w + c * f, z + c * v
    a, b = 0.75, 3 * h**2
    q2 = np.maximum(a * d**2 + b * v**2, 1e-26)
    q2prime = 2 * (a * d * dp + b * v * w)
    q2second_rest = 2 * (a * (dp**2 + d * c * w) + b * (w**2 + v * z))
    m = (n - 1) / 2
    logqprime = m * q2prime / q2
    qsecond_over_q_rest = m * q2second_rest / q2 + m * (m - 1) * (q2prime / q2)**2
    coefficient = 1 + (n - 1) * a * d**2 / q2
    rest = c * w + 2 * logqprime * dp + d * qsecond_over_q_rest
    f4 = (k * (k + 2) * d - 4 * h * (k + 1) * (w + logqprime * v) - rest) / coefficient
    return np.array([v, w, z, f4])


def collocation(n, tol=1e-9):
    beta = np.pi
    theta = np.linspace(0, beta, 301)
    t = theta / beta
    guess = np.array([1 - 2*t**2 + t**4, (-4*t + 4*t**3)/beta,
                      (-4 + 12*t**2)/beta**2, 24*t/beta**3])
    def bc(left, right, parameters):
        return np.array([left[0] - 1, left[1], left[3], right[0], right[1]])
    sol = solve_bvp(lambda x, y, par: ode(x, y, n, par[0]), bc,
                    theta, guess, p=[2 - 1/(n + 1)], tol=tol, max_nodes=30000)
    if not sol.success:
        raise RuntimeError(sol.message)
    return sol


def shooting(n, start):
    def integrate(parameters, tight=False):
        s, curvature = parameters
        return solve_ivp(lambda x, y: ode(x, y, n, s), (0, np.pi),
                         [1, 0, curvature, 0], method="DOP853",
                         rtol=2e-12 if tight else 2e-10,
                         atol=2e-13 if tight else 2e-11,
                         max_step=np.pi/180, dense_output=True)
    def residual(parameters):
        sol = integrate(parameters)
        if not sol.success:
            raise RuntimeError(sol.message)
        return sol.y[:2, -1]
    fit = root(residual, start, tol=1e-10)
    if np.max(np.abs(residual(fit.x))) > 1e-8:
        raise RuntimeError(fit.message)
    return fit.x, integrate(fit.x, tight=True)


def angular(n, s, evaluate, theta):
    """Set alpha*eps0*sigma0 = 1 and use f(0)=1; scale cancels below."""
    f, v, w, z = evaluate(theta)
    h, c = s - 1, s * (2 - s)
    d, dp = w + c*f, z + c*v
    rr, tt, rt = s*f + w, s*h*f, -h*v
    q = np.sqrt(0.75*d**2 + 3*rt**2)
    er = 0.75*q**(n - 1)*d
    et = -er
    ert = 1.5*q**(n - 1)*rt
    if n == 1:
        erprime = 0.75*dp
    else:
        qprime = (0.75*d*dp + 3*h**2*v*w)/q
        erprime = 0.75*((n - 1)*q**(n - 2)*qprime*d + q**(n - 1)*dp)
    displacement_power = n*(s - 2) + 1
    g = er/displacement_power
    gp = erprime/displacement_power
    ell = (2*ert - gp)/(displacement_power - 1)
    ellp = et - g
    energy = n/(n + 1)*q**(n + 1)
    return rr, tt, rt, q, g, gp, ell, ellp, energy


def properties(n, s, evaluate, samples=20001):
    theta = np.linspace(0, np.pi, samples)
    rr, tt, rt, q, g, gp, ell, ellp, energy = angular(n, s, evaluate, theta)
    power = n*(s - 2) + 1
    ux_r = power*g*np.cos(theta) - (gp - ell)*np.sin(theta)
    ux_t = power*ell*np.cos(theta) - (ellp + g)*np.sin(theta)
    j = 2*simpson(energy*np.cos(theta) - rr*ux_r - rt*ux_t, x=theta)
    energy_angle = 2*simpson(energy, x=theta)
    radial_energy_power = (n + 1)*(s - 2) + 2
    # H and C have no dependence on the arbitrary Airy amplitude.
    return dict(n=n, s=s, p=2 - s, H=tt[0]/j**(1/(n + 1)),
                C=energy_angle/(radial_energy_power*j),
                curvature=evaluate(np.array([0.0]))[2, 0],
                I_max_effective_stress_one=j/q.max()**(n + 1),
                J_f0_one=j, Q_min=q.min(), Q_max=q.max(),
                face_residual=max(abs(tt[-1]), abs(rt[-1])),
                radial_energy_power=radial_energy_power)


def cartesian_contour_j(n, s, evaluate, radius, delta=1e-5):
    """Second J route: differentiate Cartesian displacement numerically."""
    theta = np.linspace(0, np.pi, 16001)
    x, y = radius*np.cos(theta), radius*np.sin(theta)
    power = n*(s - 2) + 1
    def displacement(xpos, ypos):
        r = np.hypot(xpos, ypos)
        angle = np.arctan2(ypos, xpos)
        *_, g, gp, ell, ellp, energy = angular(n, s, evaluate, angle)
        return np.array([g*np.cos(angle) - ell*np.sin(angle),
                         g*np.sin(angle) + ell*np.cos(angle)])*r**power
    step = radius*delta
    ux = (displacement(x + step, y) - displacement(x - step, y))/(2*step)
    rr, tt, rt, q, g, gp, ell, ellp, energy = angular(n, s, evaluate, theta)
    traction = np.array([rr*np.cos(theta) - rt*np.sin(theta),
                         rr*np.sin(theta) + rt*np.cos(theta)])*radius**(s - 2)
    w = energy*radius**((n + 1)*(s - 2))
    integrand = (w*np.cos(theta) - np.sum(traction*ux, axis=0))*radius
    return float(2*simpson(integrand, x=theta))


def main():
    n = 6.5
    bvp = collocation(n)
    primary = properties(n, bvp.p[0], bvp.sol)
    pars, ivp = shooting(n, [bvp.p[0], bvp.y[2, 0]])
    secondary = properties(n, pars[0], ivp.sol)
    differences = {k: abs(primary[k] - secondary[k]) for k in ("p", "H", "C")}
    j_contours = {str(r): cartesian_contour_j(n, bvp.p[0], bvp.sol, r)
                  for r in (0.6, 1.0, 1.7)}
    j_relative_error = max(abs(j/primary["J_f0_one"] - 1) for j in j_contours.values())

    # Independent known solution: isotropic incompressible linear elasticity.
    elastic = collocation(1)
    elastic_computed = properties(1, elastic.p[0], elastic.sol)
    elastic_exact = dict(p=0.5, H=np.sqrt(2/(3*np.pi)), C=0.25)
    elastic_errors = {k: abs(elastic_computed[k] - value) for k, value in elastic_exact.items()}
    # External-to-our-derivation corroboration: rounded n=3 table, plane strain.
    n3 = collocation(3)
    n3_values = properties(3, n3.p[0], n3.sol)

    assert max(differences.values()) < 1e-7
    assert abs(primary["p"] - 1/(n + 1)) < 1e-8
    assert j_relative_error < 2e-6
    assert max(elastic_errors.values()) < 2e-7
    assert abs(n3_values["I_max_effective_stress_one"] - 5.51) < 0.01
    result = dict(primary=primary, shooting=secondary, absolute_differences=differences,
                  cartesian_contour_J=j_contours, contour_J_max_relative_error=j_relative_error,
                  analytic_elastic_check=dict(computed=elastic_computed, exact=elastic_exact,
                                              absolute_errors=elastic_errors),
                  rounded_n3_table_check=dict(computed_I=n3_values["I_max_effective_stress_one"],
                                               tabulated_I=5.51, computed_En=np.pi*n3_values["C"],
                                               tabulated_En=1.23))
    serialized = json.dumps(result, indent=2)
    Path(__file__).with_name("reference_results.json").write_text(serialized + "\n")
    print(serialized)


if __name__ == "__main__":
    main()
