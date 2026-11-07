"""Reproduce dimensionless strip checks and requested crack-array results.

No external data are read. The numerical finite-element solution is assembled
from the energy functional, rather than from the analytic strip displacement.
"""

import json
import os
import platform
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import scipy
import sympy as sp
from scipy.integrate import quad, solve_bvp
from scipy.linalg import solve_banded


ROOT = Path(__file__).resolve().parent
NU = 0.3
SAMPLES = [0.5, 1.0, 1.5, 2.0]
CHECKS = []
ROWS = []
LOG = []


def emit(message):
    print(message, flush=True)
    LOG.append(message)
    (ROOT / "CHECKS.txt").write_text("\n".join(LOG) + "\n")


def save_results(status, missing):
    data = {
        "status": status,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "nu": NU,
        "units": {
            "ell": "m",
            "G": "J/m^2",
            "KI": "Pa*sqrt(m)",
            "ell_over_H": "dimensionless",
            "gamma_A": "GA*Ebar/(sigma0^2*ell), dimensionless",
            "gamma_B": "GB*Ebar/(sigma0^2*ell), dimensionless",
            "kI_A": "KI_A/(sigma0*sqrt(ell)), dimensionless",
            "kI_B": "KI_B/(sigma0*sqrt(ell)), dimensionless",
        },
        "formulas": {
            "Ebar": "E/(1-nu^2)",
            "ell": "sqrt(h*Ebar/k)",
            "gamma_A": "tanh(H/(2*ell))",
            "gamma_B": "2*tanh(H/(2*ell))-tanh(H/ell)",
            "gamma_B_stable": "2*T^3/(1+T^2), T=tanh(H/(2*ell))",
            "kI": "sqrt((1-nu^2)*gamma)",
            "KII": "0 by reflection symmetry",
            "sigma_c_A": "sqrt(Gc*Ebar/(ell*gamma_A))",
            "sigma_c_B": "sqrt(Gc*Ebar/(ell*gamma_B))",
        },
        "rows": ROWS,
        "checks": CHECKS,
        "missing_work": missing,
        "verification_scope": (
            "Symbolic equilibrium, numerical BVP, direct energy quadrature, "
            "and variational 1D finite elements share the given elastic "
            "functional and translation-invariant far-field assumption. "
            "No independent 2D crack-tip extraction or experimental validation."
        ),
    }
    (ROOT / "results.json").write_text(json.dumps(data, indent=2) + "\n")


def record(name, outputs, tolerance, passed, detects):
    CHECKS.append({
        "name": name,
        "outputs": outputs,
        "tolerance": tolerance,
        "pass": bool(passed),
        "detects": detects,
    })
    emit(json.dumps(CHECKS[-1], sort_keys=True))
    save_results("checks_in_progress", ["Complete checking run and final report"])


def gammas(w):
    t = np.tanh(w / 2.0)
    return float(t), float(2.0 * t**3 / (1.0 + t**2))


def fe_strip(w, elements):
    # z=y/ell, v=Ebar*u2/(sigma0*ell), energy unit=h*sigma0^2*ell/Ebar.
    # Dimensionless energy: integral(v' + v'^2/2 + v^2/2) dz.
    dx = w / elements
    diag = np.full(elements + 1, 2.0 / dx + 2.0 * dx / 3.0)
    diag[[0, -1]] = 1.0 / dx + dx / 3.0
    off = np.full(elements, -1.0 / dx + dx / 6.0)
    matrix = np.zeros((3, elements + 1))
    matrix[0, 1:] = off
    matrix[1] = diag
    matrix[2, :-1] = off
    rhs = np.zeros(elements + 1)
    rhs[0], rhs[-1] = 1.0, -1.0
    v = solve_banded((1, 1), matrix, rhs)
    left, right = v[:-1], v[1:]
    slope = (right - left) / dx
    energy = np.sum(
        dx * slope + 0.5 * dx * slope**2
        + dx * (left**2 + left * right + right**2) / 6.0
    )
    return float(energy)


def main():
    emit("UTC " + datetime.now(timezone.utc).isoformat())
    emit(f"PID {os.getpid()}; reproduction command uses shell timeout=120s")
    emit(f"Python {platform.python_version()}; NumPy {np.__version__}; "
         f"SciPy {scipy.__version__}; SymPy {sp.__version__}")
    save_results("checks_in_progress", ["Complete checking run and final report"])

    z, a = sp.symbols("z a", real=True)
    v = -sp.sinh(z) / sp.cosh(a)
    residual = sp.simplify(sp.diff(v, z, 2) - v)
    face_plus = sp.simplify(1 + sp.diff(v, z).subs(z, a))
    face_minus = sp.simplify(1 + sp.diff(v, z).subs(z, -a))
    t = sp.symbols("t", real=True)
    positive_identity = sp.simplify(2 * t - 2 * t / (1 + t**2)
                                    - 2 * t**3 / (1 + t**2))
    record("symbolic_equilibrium_faces_and_B_identity", {
        "v_second_minus_v": str(residual),
        "sigma22_plus_face": str(face_plus),
        "sigma22_minus_face": str(face_minus),
        "B_positive_identity_residual": str(positive_identity),
    }, "exact zero", all(q == 0 for q in [residual, face_plus, face_minus,
                                           positive_identity]),
           "Wrong decay length, displacement sign, face loading, or B subtraction")

    w_symbol = sp.symbols("w", positive=True)
    ga_symbol = sp.tanh(w_symbol / 2)
    gb_symbol = 2 * ga_symbol - sp.tanh(w_symbol)
    series_a = sp.series(ga_symbol, w_symbol, 0, 8)
    series_b = sp.series(gb_symbol, w_symbol, 0, 8)
    expected_a = (w_symbol / 2 - w_symbol**3 / 24 + w_symbol**5 / 240
                  - 17 * w_symbol**7 / 40320)
    expected_b = w_symbol**3 / 4 - w_symbol**5 / 8 + 17 * w_symbol**7 / 320
    series_ok = (
        sp.expand(series_a.removeO() - expected_a) == 0
        and sp.expand(series_b.removeO() - expected_b) == 0
    )
    record("symbolic_small_spacing_series", {
        "gamma_A": str(series_a), "gamma_B": str(series_b),
    }, "exact requested coefficients", series_ok,
           "Wrong power or coefficient in dense-array shielding limits")

    quadrature_errors = []
    bvp_errors = []
    bvp_bc_errors = []
    for w in sorted({1.0 / r for r in SAMPLES}
                    | {2.0 / r for r in SAMPLES}):
        a_value = w / 2.0

        def integrand(x):
            vv = -np.sinh(x) / np.cosh(a_value)
            vp = -np.cosh(x) / np.cosh(a_value)
            return vp + 0.5 * vp**2 + 0.5 * vv**2

        energy, quad_error = quad(integrand, -a_value, a_value,
                                  epsabs=1e-12, epsrel=1e-12)
        quadrature_errors.append(abs(energy + np.tanh(a_value)))
        mesh = np.linspace(-a_value, a_value, 41)
        sol = solve_bvp(
            lambda x, state: np.array([state[1], state[0]]),
            lambda left, right: np.array([left[1] + 1, right[1] + 1]),
            mesh, np.zeros((2, len(mesh))), tol=1e-10, max_nodes=20000,
        )
        if not sol.success:
            raise RuntimeError(sol.message)
        grid = np.linspace(-a_value, a_value, 1001)
        computed = sol.sol(grid)
        expected = np.array([-np.sinh(grid) / np.cosh(a_value),
                             -np.cosh(grid) / np.cosh(a_value)])
        bvp_errors.append(float(np.max(np.abs(computed - expected))))
        bvp_bc_errors.append(float(np.max(np.abs(sol.sol(
            [-a_value, a_value])[1] + 1))))

    record("direct_full_energy_quadrature", {
        "widths_H_over_ell": sorted({1.0 / r for r in SAMPLES}
                                    | {2.0 / r for r in SAMPLES}),
        "max_abs_energy_error": max(quadrature_errors),
    }, {"max_abs_energy_error": 1e-11}, max(quadrature_errors) < 1e-11,
           "Missing foundation energy, wrong prestress term, or factor of two")
    record("numerical_strip_boundary_value_problem", {
        "max_abs_field_error": max(bvp_errors),
        "max_abs_face_residual": max(bvp_bc_errors),
    }, {"field_error": 2e-9, "face_residual": 1e-12},
           max(bvp_errors) < 2e-9 and max(bvp_bc_errors) < 1e-12,
           "Wrong strip centering, relaxation sign, or natural face conditions")

    convergence = []
    for elements in [128, 256, 512, 1024, 2048]:
        errors_a, errors_b = [], []
        for r in SAMPLES:
            w = 1.0 / r
            ua, u2a = fe_strip(w, elements), fe_strip(2 * w, elements)
            fe_a, fe_b = -ua, u2a - 2 * ua
            ga, gb = gammas(w)
            errors_a.append(abs(fe_a - ga) / ga)
            errors_b.append(abs(fe_b - gb) / gb)
        convergence.append({
            "elements_per_strip": elements,
            "max_relative_gamma_A_error": max(errors_a),
            "max_relative_gamma_B_error": max(errors_b),
        })
    ratios = [convergence[j - 1]["max_relative_gamma_B_error"]
              / convergence[j]["max_relative_gamma_B_error"]
              for j in range(1, len(convergence))]
    record("variational_finite_element_energy_release", {
        "mesh_convergence": convergence,
        "gamma_B_error_reduction_ratios": ratios,
    }, {"final_relative_A_and_B_error": 2e-6,
         "successive_B_error_ratio": [3.5, 4.5]},
           convergence[-1]["max_relative_gamma_A_error"] < 2e-6
           and convergence[-1]["max_relative_gamma_B_error"] < 2e-6
           and all(3.5 < q < 4.5 for q in ratios),
           "Energy-release signs, B old/new strip accounting, and analytic energy")

    for r in SAMPLES:
        ga, gb = gammas(1.0 / r)
        ROWS.append({"ell_over_H": r, "gamma_A": ga, "gamma_B": gb,
                     "kI_A": float(np.sqrt((1 - NU**2) * ga)),
                     "kI_B": float(np.sqrt((1 - NU**2) * gb))})
    emit("SAMPLING_ROWS " + json.dumps(ROWS, sort_keys=True))

    tiny = 1e-3
    ga_small, gb_small = gammas(tiny)
    ga_large, gb_large = gammas(40.0)
    small_a_error = abs(ga_small / (tiny / 2) - 1)
    small_b_error = abs(gb_small / (tiny**3 / 4) - 1)
    large_error = max(abs(ga_large - 1), abs(gb_large - 1))
    record("numerical_spacing_limits", {
        "small_H_over_ell": tiny,
        "gamma_A_over_H_over_2ell": ga_small / (tiny / 2),
        "gamma_B_over_H_cubed_over_4ell_cubed": gb_small / (tiny**3 / 4),
        "large_H_over_ell": 40.0,
        "large_gamma_A": ga_large, "large_gamma_B": gb_large,
    }, {"small_relative_error": 1e-6, "large_abs_error": 1e-12},
           max(small_a_error, small_b_error) < 1e-6 and large_error < 1e-12,
           "Wrong dense-array scaling or isolated-crack normalization")

    intensity_error = max(abs(row[key]**2 / (1 - NU**2) - row[gamma])
                          for row in ROWS
                          for key, gamma in [("kI_A", "gamma_A"),
                                             ("kI_B", "gamma_B")])
    record("plane_stress_intensity_normalization", {
        "max_abs_gamma_from_kI_error": intensity_error,
    }, {"max_abs_error": 1e-14}, intensity_error < 1e-14,
           "E versus Ebar conversion only; not an independent crack-tip extraction")

    ordered = all(0 < row["gamma_B"] < row["gamma_A"] < 1 for row in ROWS)
    critical_ratios = [float(np.sqrt(row["gamma_A"] / row["gamma_B"]))
                       for row in ROWS]
    record("shielding_and_critical_stress_ordering", {
        "all_0_less_B_less_A_less_1": ordered,
        "sigma_c_B_over_sigma_c_A_at_sampling_rows": critical_ratios,
    }, "strict inequalities", ordered and all(q > 1 for q in critical_ratios),
           "Interchanged geometries or incorrect signs in the split-strip release")

    passed = all(check["pass"] for check in CHECKS)
    save_results("complete" if passed else "check_failure", [] if passed else [
        "Investigate failed checks; do not treat checks as passing"
    ])
    emit("OVERALL " + ("PASS" if passed else "FAIL"))
    if not passed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
