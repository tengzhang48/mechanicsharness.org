"""Check all seven table entries against the recorded reference and 2% tolerance."""
import json
from pathlib import Path

root = Path(__file__).resolve().parents[1]
reference = json.loads((root / "reference_output.json").read_text())["route_A_angular_ODE"]
# Transcribed from the seven requested-output entries of submission/ANSWER.md.
submitted = {
    "slip_0.02_mm": 1.31813e-3,
    "slip_0.08_mm": 1.79369e-3,
    "tau_rz_0.02mm_60deg_MPa": 38.5866,
    "tau_thetaz_0.02mm_60deg_MPa": 138.011,
    "gamma_contour_0deg_mm": 0.046894,
    "gamma_contour_90deg_mm": 0.025066,
    "gamma_contour_180deg_mm": 0.013398,
}
checks = {}
for key, value in submitted.items():
    correct = reference[key]
    error = abs(value - correct) / abs(correct)
    checks[key] = {"submitted": value, "reference": correct,
                   "relative_error_percent": 100 * error,
                   "tolerance_percent": 2,
                   "tolerance_interval": [0.98 * correct, 1.02 * correct],
                   "within_tolerance": error <= 0.02}
assert all(entry["within_tolerance"] for entry in checks.values())
result = json.dumps(checks, indent=2) + "\n"
(root / "scratch" / "tolerance_checks.json").write_text(result)
print(result, end="")
