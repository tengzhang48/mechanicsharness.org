"""Compare every reported requested value with the recorded reference."""
import json
from pathlib import Path

root = Path(__file__).resolve().parent.parent
reference = json.loads((root / "reference_output.json").read_text())["route_A_angular_ODE"]
rerun = json.loads((root / "scratch/reference_rerun.json").read_text())["route_A_angular_ODE"]
reported = {
    "slip_0.02_mm": 1.318e-3,
    "slip_0.08_mm": 1.794e-3,
    "tau_rz_0.02mm_60deg_MPa": 38.59,
    "tau_thetaz_0.02mm_60deg_MPa": 138.01,
    "gamma_contour_0deg_mm": 0.0469,
    "gamma_contour_90deg_mm": 0.0251,
    "gamma_contour_180deg_mm": 0.0134,
}
out = {}
for name, value in reported.items():
    expected = reference[name]
    error = abs(value - expected) / abs(expected)
    out[name] = {
        "reported": value,
        "reference": expected,
        "relative_error_percent": 100 * error,
        "within_2_percent": bool(error <= 0.02),
        "reference_rerun_relative_error": abs(rerun[name] - expected) / abs(expected),
    }
    assert error <= 0.02
    assert value * expected > 0
    assert abs(rerun[name] / expected - 1) < 1e-12

print(json.dumps(out, indent=2))
