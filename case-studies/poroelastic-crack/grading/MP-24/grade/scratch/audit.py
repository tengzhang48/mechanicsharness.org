"""Local grading calculations; no external sources or submission mutations."""
import json
from pathlib import Path
import numpy as np

root = Path(__file__).resolve().parent.parent
reported = json.loads((root / "submission/results.json").read_text())
rerun = json.loads((root / "scratch/submission/results.json").read_text())
nu, nu_u, alpha, G = 0.2, 0.4, 0.89, 1.0
eta = alpha * (1 - 2 * nu) / (2 * (1 - nu))
face_coefficient = 16 * (1 - nu_u) * (nu_u - nu) / (np.sqrt(np.pi) * (1 - nu))
expected = {
    "K1_0": 0.807593,
    "K1_inf": 1.0,
    "V1_0": np.pi * (1 - nu_u),
    "V1_inf": np.pi * (1 - nu),
    "V1_short_coeff": 2.440220,
    "p_faces_0": 2 * (nu_u - nu) / (alpha * (1 - 2 * nu)),
    "K2_0": 0.0,
    "K2_inf": -eta,
    "V2_inf": -alpha * np.pi * (1 - 2 * nu) / 2,
}

# Exact numerical far-field integral outside a square, using s = 1/r.
# Re Z is evaluated in a rationalized form to avoid cancellation at s = 0.
nodes, weights = np.polynomial.legendre.leggauss(80)
tail_square = 0.0
R = 4.0
for lower, upper in [(0.0, np.pi / 4), (np.pi / 4, np.pi / 2)]:
    theta = lower + (nodes + 1) * (upper - lower) / 2
    theta_weights = weights * (upper - lower) / 2
    s_max = np.maximum(np.cos(theta), np.sin(theta)) / R
    s = (nodes[:, None] + 1) * s_max[None, :] / 2
    s_weights = weights[:, None] * s_max[None, :] / 2
    a = s**2 * np.exp(-2j * theta)[None, :]
    zroot = np.sqrt(1 - a)
    trace = 2 * (a / (zroot * (1 + zroot))).real
    radial_integrals = np.sum(s_weights * trace**2 / s**3, axis=0)
    tail_square += 4 * np.sum(theta_weights * radial_integrals)

# Finite-radius pressure patch: u=a*r inside R, u=a*R^2/r outside.
# Displacement and total radial traction match at R; fields vanish at infinity
# for every finite R. The inner stress persists as R tends to infinity.
lame_lambda = 2 * G * nu / (1 - 2 * nu)
radial_a = alpha / (2 * (lame_lambda + 2 * G))
inner_stress = 2 * (lame_lambda + G) * radial_a - alpha
outer_stress_at_R = -2 * G * radial_a

checks = {
    "numeric_rerun_matches": {
        name: reported[name] == rerun[name] for name in expected
    },
    "complete_results_json_equal": reported == rerun,
    "against_corrected_key": {
        name: {
            "reported": reported[name],
            "key": float(value),
            "pass": bool(abs(reported[name] - value) <= (0.01 if value == 0 else 0.01 * abs(value))),
        }
        for name, value in expected.items()
    },
    "short_time_face_only_key": float(face_coefficient),
    "short_time_tip_key": float(2.440220 - face_coefficient),
    "ANSWER_line_56_product_L_G_one": float(
        alpha * (1 - 2 * nu) / 2 * 2 * 2 * expected["p_faces_0"] * 2 / np.sqrt(np.pi)
    ),
    "ANSWER_line_57_coefficient": float(16 * (nu_u - nu) / np.sqrt(np.pi)),
    "check1_tail_used": float(np.pi / (2 * R**2)),
    "check1_tail_square_exact_quadrature": float(tail_square),
    "check1_tail_square_leading": float((np.pi / 4 + 2 / 3) / R**2),
    "check2_actual_slit_y": float((933 // 2 + 0.5) * 0.003 - 1.4),
    "pressure_patch_limit": {
        "u_over_r": float(radial_a),
        "sigma_rr_inside": float(inner_stress),
        "sigma_theta_theta_inside": float(inner_stress),
        "sigma_rr_outside_at_R": float(outer_stress_at_R),
        "minus_eta": float(-eta),
        "stress_matches_at_R": bool(np.isclose(inner_stress, outer_stress_at_R)),
        "confining_stress_crack_K": float(-eta),
        "confining_stress_crack_V": float(-eta * np.pi * (1 - nu)),
    },
}
(root / "scratch/audit_results.json").write_text(json.dumps(checks, indent=2) + "\n")
print(json.dumps(checks, indent=2))
