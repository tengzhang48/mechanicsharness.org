"""Grader audit: precise refinement change, branch shape, and requested tolerances."""
from pathlib import Path
import json
import sys
import warnings

import numpy as np
from scipy.optimize import fsolve

sys.path.insert(0, str(Path(__file__).parent / "submission"))
from notch import integrate, integrate_rk4

beta = 2 * np.pi / 3
refinement = []
for steps in (1000, 8000):
    with warnings.catch_warnings(record=True) as caught:
        solution, info, status, message = fsolve(
            lambda x: np.array(integrate_rk4(4.0, x[0], x[1], beta, nsteps=steps)),
            [1.8438247397, -1.26307900], xtol=1e-14, full_output=True,
        )
    refinement.append({
        "steps": steps, "s": float(solution[0]), "curvature": float(solution[1]),
        "endpoint_residual": float(np.max(np.abs(info["fvec"]))),
        "solver_status": int(status), "message": message,
        "warnings": [str(w.message) for w in caught],
    })

s, curvature = refinement[-1]["s"], refinement[-1]["curvature"]
sol = integrate(4.0, s, curvature, beta)
theta = np.linspace(0, beta, 5001)
p = 2 - s
ratio = (s + curvature) / (s * (s - 1))
record = {
    "refinement": refinement,
    "absolute_refinement_change_s": abs(refinement[0]["s"] - s),
    "IVP_success": bool(sol.success),
    "minimum_f_excluding_face": float(np.min(sol.sol(theta)[0, :-1])),
    "radial_energy_power": 2 - 5 * p,
    "rounded_answer_absolute_errors": {
        "p": abs(0.15618 - 0.156175260286),
        "R": abs(0.3733 - 0.373262313312),
    },
    "rerun_requested_values": {"p": p, "R": ratio},
    "rerun_absolute_errors_against_printed_key": {
        "p": abs(p - 0.156175260286),
        "R": abs(ratio - 0.373262313312),
    },
}
serialized = json.dumps(record, indent=2)
Path(__file__).with_name("audit_results.json").write_text(serialized + "\n")
print(serialized)
