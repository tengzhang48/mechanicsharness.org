"""Finite-tangent material: closed reduction versus tensor/energy checks.

Run: python V2/solve_reference.py
Requires only numpy and scipy. Units: mm, MPa, N.
"""

import json
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

E = 70000.0
NU = 0.30
SIGMA_Y = 10.0
J = 0.80
CZZ = -0.005
R_GEOM = 0.20
MAX_STRAIN = 0.05


def polar_stress_coefficient(theta, k):
    c, s = np.cos(theta / 2), np.sin(theta / 2)
    return k / np.sqrt(2 * np.pi) * np.array([c * (1 + s*s), c**3, s*c*c])


def cartesian_stress(r, theta, k):
    rr, tt, rt = polar_stress_coefficient(theta, k) / np.sqrt(r)
    q = np.array([[np.cos(theta), -np.sin(theta)],
                  [np.sin(theta), np.cos(theta)]])
    a = np.zeros((3, 3))
    a[:2, :2] = q @ np.array([[rr, rt], [rt, tt]]) @ q.T
    return a


def compliance_tensor(stress, t, retain_intercept=False):
    ident = np.eye(3)
    dev = stress - np.trace(stress) / 3 * ident
    plastic_factor = 1.0
    if retain_intercept:
        se = np.sqrt(1.5 * np.sum(dev * dev))
        plastic_factor = max(1 - SIGMA_Y / se, 0) if se else 0.0
    return ((1 + NU) / E * stress - NU / E * np.trace(stress) * ident
            + 1.5 * (1 / t - 1 / E) * plastic_factor * dev)


def displacement(r, theta, k, t):
    nut = 0.5 - (0.5 - NU) * t / E
    mu = t / (2 * (1 + nut))
    kap = (3 - nut) / (1 + nut)
    factor = k / (2 * mu) * np.sqrt(r / (2 * np.pi)) * (kap - np.cos(theta))
    return factor * np.array([np.cos(theta / 2), np.sin(theta / 2)])


def closed_route():
    l = -CZZ / np.sqrt(J / (2 * np.pi))
    b = (1 - 2 * NU) / E
    sqt = 2 / (l + np.sqrt(l*l + 4*b))
    t = sqt*sqt
    nut = 0.5 - (0.5 - NU) * t / E
    k = np.sqrt(J*t)
    rr, tt, rt = polar_stress_coefficient(np.pi / 3, k)
    rr0, tt0, _ = polar_stress_coefficient(0, k)
    ahead_strains = np.array([(rr0 - nut*tt0) / t,
                              (tt0 - nut*rr0) / t,
                              -nut*(rr0 + tt0) / t])
    rmin = (np.max(np.abs(ahead_strains)) / MAX_STRAIN)**2
    rmax_stress = (k / (np.sqrt(2*np.pi) * 10 * SIGMA_Y))**2
    cy90 = (3 - nut) * np.sqrt(J / (4*np.pi*t))
    out = {"T_MPa": float(t), "nu_effective": float(nut),
           "C_rr_60deg_sqrtmm": float((rr - nut*tt) / t),
           "C_thetatheta_60deg_sqrtmm": float((tt - nut*rr) / t),
           "C_engineering_rtheta_60deg_sqrtmm": float(2*(1 + nut)*rt / t),
           "C_y_90deg_sqrtmm": float(cy90),
           "r_min_mm": float(rmin),
           "r_max_mm": float(min(rmax_stress, R_GEOM))}
    return out


def tensor_route():
    # Calibrate directly through the original 3-D high-stress compliance tensor.
    def residual(t):
        stress = cartesian_stress(1.0, 0.0, np.sqrt(J*t))
        return compliance_tensor(stress, t)[2, 2] - CZZ
    t = brentq(residual, 1.0, E, xtol=1e-10)
    k = np.sqrt(J*t)
    theta = np.pi / 3
    eps = compliance_tensor(cartesian_stress(1.0, theta, k), t)
    q = np.array([[np.cos(theta), -np.sin(theta)],
                  [np.sin(theta), np.cos(theta)]])
    ep = q.T @ eps[:2, :2] @ q

    def stress_criterion(r):
        stress = cartesian_stress(r, 0, k)
        dev = stress - np.trace(stress) / 3 * np.eye(3)
        return np.sqrt(1.5 * np.sum(dev*dev)) - 10*SIGMA_Y

    def strain_criterion(r):
        eps0 = compliance_tensor(cartesian_stress(r, 0, k), t)
        return np.max(np.abs(np.linalg.eigvalsh(eps0))) - MAX_STRAIN

    rmin = brentq(strain_criterion, 1e-6, R_GEOM, xtol=1e-14)
    rmax = min(brentq(stress_criterion, 1e-6, 5.0, xtol=1e-14), R_GEOM)
    # Calculate nu_eff from two independent tensor experiments, not its formula.
    uniaxial = np.diag([1.0, 0.0, 0.0])
    strain_uni = compliance_tensor(uniaxial, t)
    nu_eff = -strain_uni[1, 1] / strain_uni[0, 0]
    cy = displacement(1.0, np.pi / 2, k, t)[1]
    out = {"T_MPa": float(t), "nu_effective": float(nu_eff),
           "C_rr_60deg_sqrtmm": float(ep[0, 0]),
           "C_thetatheta_60deg_sqrtmm": float(ep[1, 1]),
           "C_engineering_rtheta_60deg_sqrtmm": float(2*ep[0, 1]),
           "C_y_90deg_sqrtmm": float(cy),
           "r_min_mm": float(rmin), "r_max_mm": float(rmax)}
    return out


def contour_integral(r, t, k):
    """Independent J evaluation using W and finite-difference u_,x."""
    step = r * 2e-5

    def u_cart(x, y):
        return displacement(np.hypot(x, y), np.arctan2(y, x), k, t)

    def integrand(theta):
        stress = cartesian_stress(r, theta, k)
        eps = compliance_tensor(stress, t)
        w = np.sum(stress * eps) / 2
        nx, ny = np.cos(theta), np.sin(theta)
        x, y = r*nx, r*ny
        dux = (u_cart(x + step, y) - u_cart(x - step, y)) / (2*step)
        traction = stress[:2, :2] @ np.array([nx, ny])
        return r * (w * nx - traction @ dux)
    return quad(integrand, -np.pi, np.pi, epsabs=1e-10, epsrel=1e-10)[0]


def main():
    a, b = closed_route(), tensor_route()
    errors = {key: abs(a[key] - b[key]) / abs(a[key]) for key in a}
    t = a["T_MPa"]
    k = np.sqrt(J*t)
    js = {str(r): float(contour_integral(r, t, k)) for r in (0.015, 0.04, 0.10)}

    # Direct check of epsilon = sym(grad u) at an off-axis point.
    r, theta = 0.04, np.pi / 3
    x, y = r*np.cos(theta), r*np.sin(theta)
    step = r*2e-5
    def u_xy(xp, yp):
        return displacement(np.hypot(xp, yp), np.arctan2(yp, xp), k, t)
    grad = np.column_stack(((u_xy(x+step, y)-u_xy(x-step, y))/(2*step),
                            (u_xy(x, y+step)-u_xy(x, y-step))/(2*step)))
    strain_displacement = (grad + grad.T) / 2
    strain_tensor = compliance_tensor(cartesian_stress(r, theta, k), t)[:2, :2]
    compatibility_error = np.max(np.abs(strain_displacement - strain_tensor))

    # Independent material checks: retain the bilinear intercept in uniaxial tension,
    # and recover the original elastic tensor exactly when T=E.
    uni_stress = 30.0
    uni_strain = compliance_tensor(np.diag([uni_stress, 0, 0]), t, True)[0, 0]
    uni_expected = SIGMA_Y / E + (uni_stress - SIGMA_Y) / t
    test = np.array([[40., 7., 0.], [7., -12., 0.], [0., 0., 3.]])
    elastic_expected = (1+NU)/E*test - NU/E*np.trace(test)*np.eye(3)
    elastic_error = np.max(np.abs(compliance_tensor(test, E) - elastic_expected))

    out = {"route_A_closed_reduction": a, "route_B_full_tensor_and_roots": b,
           "relative_errors_A_vs_B": errors, "J_circle_N_per_mm": js,
           "strain_displacement_max_absolute_error": float(compatibility_error),
           "uniaxial_full_law_strain": float(uni_strain),
           "uniaxial_bilinear_curve_strain": float(uni_expected),
           "T_equals_E_elastic_tensor_error": float(elastic_error)}
    assert max(errors.values()) < 1e-9
    assert all(abs(v/J - 1) < 2e-8 for v in js.values())
    assert compatibility_error < 2e-9
    assert abs(uni_strain - uni_expected) < 1e-14
    assert elastic_error < 1e-14
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
