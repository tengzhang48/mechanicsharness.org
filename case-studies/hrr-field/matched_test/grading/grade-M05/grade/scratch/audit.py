import contextlib
import json
import math
import re
import runpy
from pathlib import Path

from scipy.integrate import quad

base = Path(__file__).resolve().parent
reference = json.loads((base / "reference.stdout").read_text())["route_B_phase_formula"]
stdout = (base / "submission.stdout").read_text()
checks = [
    ("slip_0.02_mm", 0.0013181, r"Delta w\(r=0\.0200\) = ([\d.eE+-]+)"),
    ("slip_0.08_mm", 0.0017937, r"Delta w\(r=0\.0800\) = ([\d.eE+-]+)"),
    ("tau_rz_0.02mm_60deg_MPa", 38.59, r"tau_rz  = ([\d.eE+-]+)"),
    ("tau_thetaz_0.02mm_60deg_MPa", 138.01, r"tau_thz = ([\d.eE+-]+)"),
    ("gamma_contour_0deg_mm", 0.046894, r"theta=0 : r = ([\d.eE+-]+)"),
    ("gamma_contour_90deg_mm", 0.025066, r"theta=90 : r = ([\d.eE+-]+)"),
    ("gamma_contour_180deg_mm", 0.013398, r"theta=180 : r = ([\d.eE+-]+)"),
]
out = {"tolerance_checks": {}}
for key, answer, pattern in checks:
    computed = float(re.search(pattern, stdout)[1])
    ref = reference[key]
    answer_error = abs(answer - ref) / abs(ref)
    code_error = abs(computed - ref) / abs(ref)
    assert answer_error <= 0.02 and code_error <= 0.02
    assert math.copysign(1, answer) == math.copysign(1, ref)
    out["tolerance_checks"][key] = {
        "main_answer_value": answer,
        "rerun_value": computed,
        "reference": ref,
        "main_answer_relative_error_percent": 100 * answer_error,
        "rerun_relative_error_percent": 100 * code_error,
        "within_2_percent": True,
    }

with (base / "audit_import.stdout").open("w") as log, contextlib.redirect_stdout(log):
    s = runpy.run_path(str(base / "submission_solve.py"))

out["dimensional_circle_flux_N_per_mm"] = {
    str(r): s["J_circle"](r) for r in (0.005, 0.02, 0.1, 0.3)
}
assert all(abs(j - 0.12) < 1e-12 for j in out["dimensional_circle_flux_N_per_mm"].values())
out["angular_ODE_face_derivative"] = float(s["sol"].y[1, -1])
out["stored_work_absolute_difference"] = abs(s["Wnum"] - s["Wform"])

# Independently evaluate the formula used by the auxiliary debug script:
# its w_,x lacks A. This does not occur in the submitted main solver.
def debug_integrand(theta, r):
    _, f, fp, h = s["field"](theta)
    a, amplitude, n = s["a"], s["A"], s["n"]
    gamma = amplitude * r ** (a - 1) * h
    tau = s["tau0"] * (gamma / s["gam0"]) ** (1 / n)
    work = n / (n + 1) * tau * gamma
    tau_r = tau * a * f / h
    wx_without_amplitude = r ** (a - 1) * (a * f * math.cos(theta) - fp * math.sin(theta))
    return r * (work * math.cos(theta) - tau_r * wx_without_amplitude)

out["auxiliary_debug_flux_missing_A_N_per_mm"] = {
    str(r): 2 * quad(lambda theta: debug_integrand(theta, r), 0, math.pi)[0]
    for r in (1.0, 0.02)
}

(base / "audit.json").write_text(json.dumps(out, indent=2) + "\n")
print(json.dumps(out, indent=2))
