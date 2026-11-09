"""Compare the submitted answer and field with the supplied grading key."""

import contextlib
import io
import json
import runpy
from pathlib import Path

import numpy as np

root = Path(__file__).resolve().parent.parent
reference = json.loads((root / "reference_output.json").read_text())
rerun = json.loads((root / "scratch/reference_rerun.json").read_text())
key = reference["route_A_angular_ODE"]
answer = {
    "slip_0.02_mm": 1.31813e-3,
    "slip_0.08_mm": 1.79369e-3,
    "tau_rz_0.02mm_60deg_MPa": 38.587,
    "tau_thetaz_0.02mm_60deg_MPa": 138.011,
    "gamma_contour_0deg_mm": 4.68944e-2,
    "gamma_contour_90deg_mm": 2.50661e-2,
    "gamma_contour_180deg_mm": 1.33984e-2,
}
comparisons = {
    name: {
        "submitted": value,
        "reference": key[name],
        "absolute_relative_error_percent": 100 * abs(value / key[name] - 1),
        "within_2_percent": abs(value / key[name] - 1) <= 0.02,
    }
    for name, value in answer.items()
}
with contextlib.redirect_stdout(io.StringIO()):
    submitted = runpy.run_path(str(root / "scratch/solve.py"))
f90, fp90 = submitted["f_fp"](np.pi / 2)
out = {
    "requested_numbers": comparisons,
    "largest_relative_error_percent": max(
        row["absolute_relative_error_percent"] for row in comparisons.values()
    ),
    "reference_rerun_max_change": max(
        abs(key[name] - rerun["route_A_angular_ODE"][name]) for name in key
    ),
    "angular_value_at_90_degrees": {
        "f": f90,
        "fprime": fp90,
        "g": submitted["gfun"](np.pi / 2)[0],
    },
    "face_fprime": submitted["fppi"],
    "J_from_submitted_cartesian_check": {
        str(r): submitted["J_circle"](r) for r in (0.5, 1.0, 2.0)
    },
}
assert all(row["within_2_percent"] for row in comparisons.values())
print(json.dumps(out, indent=2))
