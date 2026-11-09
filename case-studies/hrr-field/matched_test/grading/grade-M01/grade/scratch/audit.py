"""Local grading audit; this is grader code, not submitted code."""

import json
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

ROOT = Path(__file__).resolve().parents[1]
E, NU, SY, J = 70000.0, 0.30, 10.0, 0.800


def tensor_law(stress, t, full=False):
    trace = np.trace(stress)
    dev = stress - trace * np.eye(3) / 3.0
    equivalent = np.sqrt(1.5 * np.sum(dev * dev))
    factor = max(1.0 - SY / equivalent, 0.0) if full and equivalent else 1.0
    return (
        (1.0 + NU) * stress / E
        - NU * trace * np.eye(3) / E
        + 1.5 * (1.0 / t - 1.0 / E) * factor * dev
    )


def polar_stress(r, theta, amplitude):
    rr = amplitude * (1.25 * np.cos(theta / 2) - 0.25 * np.cos(1.5 * theta))
    tt = amplitude * (0.75 * np.cos(theta / 2) + 0.25 * np.cos(1.5 * theta))
    rt = amplitude * (0.25 * np.sin(theta / 2) + 0.25 * np.sin(1.5 * theta))
    return np.array([[rr, rt, 0.0], [rt, tt, 0.0], [0.0, 0.0, 0.0]]) / np.sqrt(r)


def cartesian_stress(r, theta, amplitude):
    q = np.array([[np.cos(theta), -np.sin(theta), 0.0],
                  [np.sin(theta), np.cos(theta), 0.0],
                  [0.0, 0.0, 1.0]])
    return q @ polar_stress(r, theta, amplitude) @ q.T


def calibration(t):
    amplitude = np.sqrt(J * t / (2.0 * np.pi))
    return tensor_law(np.diag([amplitude, amplitude, 0.0]), t)[2, 2] + 0.005


t = brentq(calibration, 1.0, E, xtol=1e-10)
nut = 0.5 - (0.5 - NU) * t / E
g = t / (2.0 * (1.0 + nut))
kap = (3.0 - nut) / (1.0 + nut)
amplitude = np.sqrt(J * t / (2.0 * np.pi))
eps60 = tensor_law(polar_stress(1.0, np.pi / 3.0, amplitude), t)
cy90 = amplitude / (2.0 * g) * kap / np.sqrt(2.0)
ahead = np.linalg.eigvalsh(tensor_law(polar_stress(1.0, 0.0, amplitude), t))
calculated = {
    "T_MPa": float(t),
    "C_rr_60deg_sqrtmm": float(eps60[0, 0]),
    "C_thetatheta_60deg_sqrtmm": float(eps60[1, 1]),
    "C_engineering_rtheta_60deg_sqrtmm": float(2.0 * eps60[0, 1]),
    "C_y_90deg_sqrtmm": float(cy90),
    "r_min_mm": float((np.max(np.abs(ahead)) / 0.05) ** 2),
    "r_max_mm": float(min((amplitude / 100.0) ** 2, 0.2)),
}
reported = {
    "T_MPa": 4816.47,
    "C_rr_60deg_sqrtmm": 3.9420e-3,
    "C_thetatheta_60deg_sqrtmm": 6.3318e-4,
    "C_engineering_rtheta_60deg_sqrtmm": 5.7311e-3,
    "C_y_90deg_sqrtmm": 9.1390e-3,
    "r_min_mm": 0.010000,
    "r_max_mm": 0.061325,
}
reference = json.loads((ROOT / "reference_output.json").read_text())["route_A_closed_reduction"]
comparisons = {}
for key, value in reported.items():
    relative_error = abs(value - reference[key]) / abs(reference[key])
    tolerance = 0.01 if key == "T_MPa" else 0.02
    comparisons[key] = {
        "reported": value,
        "reference": reference[key],
        "relative_error_percent": 100.0 * relative_error,
        "tolerance_percent": 100.0 * tolerance,
        "within_tolerance": bool(relative_error <= tolerance),
    }


def direct_j(radius):
    def integrand(theta):
        stress = cartesian_stress(radius, theta, amplitude)
        strain = tensor_law(stress, t)
        energy = 0.5 * np.sum(stress * strain)
        normal = np.array([np.cos(theta), np.sin(theta)])
        v = np.array([np.cos(theta / 2.0), np.sin(theta / 2.0)])
        vp = np.array([-np.sin(theta / 2.0), np.cos(theta / 2.0)])
        f = kap - np.cos(theta)
        # The compatible mode-I Cartesian displacement completes the stated u_r.
        du_dr = amplitude / (4.0 * g * np.sqrt(radius)) * f * v
        du_dtheta = amplitude * np.sqrt(radius) / (2.0 * g) * (
            np.sin(theta) * v + 0.5 * f * vp
        )
        du_dx = np.cos(theta) * du_dr - np.sin(theta) / radius * du_dtheta
        traction = stress[:2, :2] @ normal
        return radius * (energy * normal[0] - traction @ du_dx)

    return quad(integrand, -np.pi, np.pi, epsabs=1e-12, epsrel=1e-12)[0]


full_law_check = {}
for radius in (1e-2, 1e-3, 1e-4, 1e-5):
    stress = polar_stress(radius, 0.0, amplitude)
    full_law_check[str(radius)] = float(np.sqrt(radius) * tensor_law(stress, t, full=True)[2, 2])

bulk = E / (3.0 * (1.0 - 2.0 * NU))
hyd = (1.0 - 2.0 * NU) / (3.0 * E)
output = {
    "comparisons": comparisons,
    "independently_calculated": calculated,
    "auxiliary": {
        "nu_effective": float(nut),
        "G_effective_MPa": float(g),
        "kappa_effective": float(kap),
        "K_A_MPa_sqrtmm": float(amplitude),
        "bulk_modulus_MPa": float(bulk),
        "A_hyd_per_MPa": float(hyd),
        "1_over_3_bulk": float(1.0 / (3.0 * bulk)),
        "1_over_9_bulk": float(1.0 / (9.0 * bulk)),
        "ahead_principal_strain_coefficients": ahead.tolist(),
    },
    "recomputed_J_N_per_mm": {str(r): float(direct_j(r)) for r in (0.05, 0.1, 0.2, 1.0)},
    "recomputed_full_law_Czz": full_law_check,
}
assert all(entry["within_tolerance"] for entry in comparisons.values())
assert all(abs(v - J) < 1e-10 for v in output["recomputed_J_N_per_mm"].values())
assert all(abs(calculated[key] - reference[key]) < 1e-8 for key in calculated)
print(json.dumps(output, indent=2))
