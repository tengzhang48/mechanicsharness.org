"""Record tolerance checks and reproduce the submitted independent checks."""
import contextlib
import io
import json
import runpy
from pathlib import Path

import numpy as np

root = Path(__file__).resolve().parent.parent
reference = json.loads((root / "reference_output.json").read_text())["route_A_closed_reduction"]
rerun = json.loads((root / "scratch/reference_rerun.json").read_text())["route_A_closed_reduction"]
assert all(abs(rerun[k] / reference[k] - 1) < 1e-12 for k in reference)
with contextlib.redirect_stdout(io.StringIO()):
    submitted = runpy.run_path(str(root / "scratch/compute.py"))

# Conservatively use the coarsest values explicitly reported in ANSWER.md's table.
reported = {
    "T_MPa": 4816.0,
    "C_rr_60deg_sqrtmm": 3.942e-3,
    "C_thetatheta_60deg_sqrtmm": 6.332e-4,
    "C_engineering_rtheta_60deg_sqrtmm": 5.731e-3,
    "C_y_90deg_sqrtmm": 9.139e-3,
    "r_min_mm": 0.0100,
    "r_max_mm": 0.0613,
}
computed = {
    "T_MPa": submitted["T"],
    "C_rr_60deg_sqrtmm": np.sqrt(submitted["r"]) * submitted["err"],
    "C_thetatheta_60deg_sqrtmm": np.sqrt(submitted["r"]) * submitted["ett"],
    "C_engineering_rtheta_60deg_sqrtmm": np.sqrt(submitted["r"]) * 2 * submitted["ert"],
    "C_y_90deg_sqrtmm": submitted["Cy"],
    "r_min_mm": submitted["r2"],
    "r_max_mm": submitted["rhi"],
}
checks = {}
for key, value in reported.items():
    tolerance = 0.01 if key == "T_MPa" else 0.02
    error = abs(value - reference[key]) / abs(reference[key])
    checks[key] = {
        "reported": value,
        "reference": reference[key],
        "relative_error_percent": 100 * error,
        "tolerance_percent": 100 * tolerance,
        "within_tolerance": bool(error <= tolerance),
        "code_value": float(computed[key]),
        "code_relative_error": float(abs(computed[key] / reference[key] - 1)),
    }
    assert error <= tolerance
    assert checks[key]["code_relative_error"] < 1e-12

constitutive = np.array([submitted["Crr_f"], submitted["Ctt_f"], submitted["Cgt_f"]])
kinematic = np.sqrt(submitted["r0"]) * np.array([
    submitted["err_fd"], submitted["ett_fd"], 2 * submitted["ert_fd"]
])
fd_error = float(np.max(np.abs((kinematic - constitutive) / constitutive)))
assert fd_error < 1e-6
j_values = {str(rad): float(submitted["Jnum"](
    rad, submitted["K"], submitted["Ee"], submitted["nue"]
)) for rad in (0.001, 0.01, 0.1)}
assert all(abs(value / 0.8 - 1) < 1e-6 for value in j_values.values())
assert submitted["r2"] <= submitted["rhi"]
assert 0 < submitted["T"] < submitted["E"]
audit = {
    "requested_number_checks": checks,
    "fd_strain_max_relative_error": fd_error,
    "J_circle_N_per_mm": j_values,
    "interval_nonempty": True,
    "T_admissible": True,
}
out = json.dumps(audit, indent=2)
(root / "scratch/grading_audit.json").write_text(out + "\n")
print(out)
