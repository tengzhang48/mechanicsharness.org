"""Unequal splitting arms with softening connectors and intrinsic front work.

Primary route: a beam-opening boundary-value problem with an unknown length.
Independent route: exact hyperbolic/trigonometric solution plus energy balance.
All internal units are SI; reported M is the total couple on EACH arm.
"""

import json
from math import cos, cosh, pi, sin, sinh, sqrt

import numpy as np
from scipy.integrate import quad, solve_bvp
from scipy.optimize import brentq


WIDTH = 25.0e-3
E1, H1 = 110.0e9, 0.60e-3
E2, H2 = 190.0e9, 0.90e-3
SIGMA_B = 0.20e6
DELTA_C = 0.20e-3
GAMMA_0 = 1.40
D1 = E1 * H1**3 / 12.0  # bending stiffness per unit width, N m
D2 = E2 * H2**3 / 12.0
C = 1.0 / D1 + 1.0 / D2
ELL = (DELTA_C / (C * SIGMA_B))**0.25


def solve_opening(gamma, traction=lambda u: 1.0 - u, tol=1.0e-10):
    """u=delta/delta_c on t=x/L. gamma=Gamma_0/(sigma_b*delta_c)."""
    kappa = sqrt(2.0 * gamma)

    def ode(t, y, parameters):
        lam = parameters[0]
        return np.vstack((y[1], y[2], y[3], -lam**4 * traction(y[0])))

    def bc(left, right, parameters):
        lam = parameters[0]
        return np.array([left[0], left[1], left[2] - kappa * lam**2,
                         right[0] - 1.0, right[3]])

    t = np.linspace(0.0, 1.0, 65)
    lam_guess = 1.7
    a2 = 0.5 * kappa * lam_guess**2
    a3 = 4.0 * (1.0 - a2) / 3.0
    a4 = (a2 - 1.0) / 3.0
    guess = np.vstack((a2*t**2 + a3*t**3 + a4*t**4,
                       2.0*a2*t + 3.0*a3*t**2 + 4.0*a4*t**3,
                       2.0*a2 + 6.0*a3*t + 12.0*a4*t**2,
                       6.0*a3 + 24.0*a4*t))
    sol = solve_bvp(ode, bc, t, guess, p=[lam_guess], tol=tol,
                    bc_tol=tol, max_nodes=10000)
    if not sol.success:
        raise RuntimeError(sol.message)
    lam = float(sol.p[0])
    rear_curvature = float(sol.sol(1.0)[2]) / lam**2  # derivative in x/ell
    m_per_width = DELTA_C * rear_curvature / (C * ELL**2)
    grid = np.linspace(0.0, 1.0, 1001)
    profile = sol.sol(grid)
    assert np.min(profile[0]) >= -1.0e-9
    assert np.max(profile[0]) <= 1.0 + 1.0e-9
    assert np.min(profile[1]) >= -1.0e-9
    assert np.min(profile[2]) >= -1.0e-9
    data = {"lambda": lam, "L_m": lam * ELL,
            "M_each_arm_N_m": WIDTH * m_per_width,
            "delta_mid_m": DELTA_C * float(sol.sol(0.5)[0]),
            "rear_dimensionless_curvature": rear_curvature,
            "boundary_max_abs_residual":
                float(np.max(np.abs(bc(sol.sol(0.0), sol.sol(1.0), sol.p)))),
            "max_RMS_ODE_residual": float(np.max(sol.rms_residuals)),
            "mesh_nodes": int(sol.x.size)}
    return data, sol


def exact_softening_solution(gamma):
    """Solve the first positive exact characteristic root, independently of BVP."""
    kappa = sqrt(2.0 * gamma)

    def characteristic(lam):
        return 1.0 + cosh(lam)*cos(lam) - kappa*sinh(lam)*sin(lam)

    cantilever_root = brentq(lambda z: 1.0 + cosh(z)*cos(z),
                            pi/2.0, 2.0, xtol=1.0e-14)
    lam = brentq(characteristic, 1.0e-8, cantilever_root + 1.0e-10,
                 xtol=1.0e-14, rtol=1.0e-14)
    a3 = ((sinh(lam) + sin(lam)) - kappa*(sinh(lam) - sin(lam))) / (
        cosh(lam) + cos(lam))

    def exact_profile(z):
        return (1.0 - 0.5*(np.cosh(z) + np.cos(z))
                + 0.5*kappa*(np.cosh(z) - np.cos(z))
                + 0.5*a3*(np.sinh(z) - np.sin(z)))

    work = SIGMA_B * DELTA_C * quad(lambda u: 1.0 - u, 0.0, 1.0,
                                    epsabs=1.0e-13, epsrel=1.0e-13)[0]
    intrinsic_work = gamma * SIGMA_B * DELTA_C
    m_energy = sqrt(2.0 * (intrinsic_work + work) / C)
    data = {"lambda": lam, "L_m": ELL * lam,
            "M_each_arm_N_m": WIDTH * m_energy,
            "delta_mid_m": DELTA_C * float(exact_profile(lam/2.0)),
            "dimensionless_tip_third_derivative": a3,
            "characteristic_residual": characteristic(lam),
            "full_connector_work_J_m2": work,
            "zero_intrinsic_energy_root": cantilever_root}
    return data, exact_profile


def main():
    gamma = GAMMA_0 / (SIGMA_B * DELTA_C)
    first, bvp = solve_opening(gamma)
    second, exact_profile = exact_softening_solution(gamma)
    errors = {key: abs(first[key] / second[key] - 1.0)
              for key in ("L_m", "M_each_arm_N_m", "delta_mid_m")}
    assert max(errors.values()) < 1.0e-8
    grid = np.linspace(0.0, 1.0, 1001)
    profile_error = float(np.max(np.abs(bvp.sol(grid)[0] - exact_profile(second["lambda"]*grid))))
    assert profile_error < 1.0e-8

    m = first["M_each_arm_N_m"] / WIDTH
    j = 0.5 * C * m*m
    energy_error = abs(j / (GAMMA_0 + second["full_connector_work_J_m2"]) - 1.0)
    assert energy_error < 1.0e-8

    # Zero intrinsic front work gives the familiar first cantilever root.
    zero_tip, _ = solve_opening(0.0)
    zero_tip_error = abs(zero_tip["lambda"] / second["zero_intrinsic_energy_root"] - 1.0)
    assert zero_tip_error < 1.0e-8

    # Change the connector law to constant pressure and set Gamma_0=0.
    # Then delta(L)=C*sigma_b*L^4/8, so lambda=8^(1/4) exactly.
    rectangle, _ = solve_opening(0.0, traction=lambda u: np.ones_like(u))
    rectangle_error = abs(rectangle["lambda"] / 8.0**0.25 - 1.0)
    assert rectangle_error < 1.0e-8

    # Local first integral throughout the softening zone:
    # 0.5*u''^2 - u'*u''' = gamma + u - u^2/2.
    y = bvp.sol(grid)
    lam = first["lambda"]
    invariant_error = float(np.max(np.abs(
        0.5*(y[2]/lam**2)**2 - (y[1]/lam)*(y[3]/lam**3)
        - gamma - y[0] + 0.5*y[0]**2)))
    assert invariant_error < 1.0e-8

    result = {
        "inputs_SI": {"width": WIDTH, "E1": E1, "h1": H1, "E2": E2, "h2": H2,
                      "sigma_b": SIGMA_B, "delta_c": DELTA_C, "Gamma_0": GAMMA_0},
        "derived": {"D1_per_width_N_m": D1, "D2_per_width_N_m": D2,
                    "C_per_width_inverse_N_m": C, "ell_m": ELL,
                    "gamma": gamma},
        "answer": {"M_each_arm_N_m": first["M_each_arm_N_m"],
                   "L_mm": first["L_m"] * 1.0e3,
                   "delta_at_L_over_2_um": first["delta_mid_m"] * 1.0e6},
        "BVP": first,
        "independent_exact_and_energy_route": second,
        "checks": {"relative_differences": errors,
                   "maximum_dimensionless_profile_difference": profile_error,
                   "energy_balance_relative_error": energy_error,
                   "zero_front_work_cantilever_root_relative_error": zero_tip_error,
                   "constant_traction_length_relative_error": rectangle_error,
                   "maximum_first_integral_abs_error": invariant_error},
        "geometry_diagnostics": {"L_over_thicker_arm": first["L_m"] / H2,
                                 "relative_rotation_at_rear_rad":
                                     DELTA_C / first["L_m"] * float(bvp.sol(1.0)[1])}
    }
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
