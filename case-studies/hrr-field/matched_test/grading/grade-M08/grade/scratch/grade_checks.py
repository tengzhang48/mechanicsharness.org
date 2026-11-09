"""Grading evidence: displayed-value tolerances and the full leading strain tensor."""
import json
from pathlib import Path

import numpy as np

base = Path(__file__).resolve().parent
reference = json.loads((base / "reference_rerun.json").read_text())["route_A_closed_reduction"]
reported = {
    "T_MPa": 4816.47,
    "C_rr_60deg_sqrtmm": 0.003942,
    "C_thetatheta_60deg_sqrtmm": 0.000633,
    "C_engineering_rtheta_60deg_sqrtmm": 0.005731,
    "C_y_90deg_sqrtmm": 0.009139,
    "r_min_mm": 0.02336,
    "r_max_mm": 0.06133,
}
checks = {}
for key, submitted in reported.items():
    target = reference[key]
    relative_error = abs(submitted - target) / abs(target)
    tolerance = 0.01 if key == "T_MPa" else 0.02
    checks[key] = {
        "submitted": submitted,
        "reference": target,
        "relative_error_percent": 100 * relative_error,
        "tolerance_percent": 100 * tolerance,
        "within_tolerance": bool(relative_error <= tolerance),
    }

E, nu, J, max_strain = 70000.0, 0.30, 0.800, 0.05
T = reference["T_MPa"]
K = np.sqrt(J * T)
A = K / np.sqrt(2 * np.pi)
stress_coeff = np.diag([A, A, 0.0])
dev = stress_coeff - np.trace(stress_coeff) / 3 * np.eye(3)
strain_coeff = (
    (1 + nu) / E * stress_coeff
    - nu / E * np.trace(stress_coeff) * np.eye(3)
    + 1.5 * (1 / T - 1 / E) * dev
)
principal_coeff = np.linalg.eigvalsh(strain_coeff)
a = 1.5 / T + (nu - 0.5) / E
b = (0.5 - nu) / E - 0.5 / T
correct_lower = (np.max(np.abs(principal_coeff)) / max_strain) ** 2
submitted_lower_unrounded = (A * max(abs(a), abs(2 * b)) / max_strain) ** 2
upper = min(J * T / (2 * np.pi * 100.0**2), 0.200)
assert abs(correct_lower - reference["r_min_mm"]) < 1e-14
assert abs(upper - reference["r_max_mm"]) < 1e-14
assert checks["r_min_mm"]["within_tolerance"] is False
assert all(v["within_tolerance"] for k, v in checks.items() if k != "r_min_mm")

out = {
    "displayed_answer_tolerances": checks,
    "leading_principal_strain_coefficients_theta0": np.diag(strain_coeff).tolist(),
    "strain_per_p": {"a": a, "2b": 2 * b, "correct_in_plane_a_plus_2b": a + 2 * b},
    "lower_endpoint_from_original_leading_tensor": float(correct_lower),
    "lower_endpoint_reproducing_submitted_single_omission": float(submitted_lower_unrounded),
    "correct_interval_mm": [float(correct_lower), float(upper)],
    "logarithmic_T_sensitivity_to_abs_Czz": float(
        -2 * reference["nu_effective"] / (1 - reference["nu_effective"])
    ),
}
print(json.dumps(out, indent=2))
