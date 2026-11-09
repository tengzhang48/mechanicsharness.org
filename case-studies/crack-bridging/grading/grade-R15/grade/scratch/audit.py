"""Local grading audit; submission stays unchanged."""
import contextlib
import io
import json
import runpy
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

folder = Path(__file__).resolve().parent.parent
with contextlib.redirect_stdout(io.StringIO()):
    submitted = runpy.run_path(str(folder / "scratch" / "solve.py"))
key = json.loads((folder / "reference_results.json").read_text())["answer"]
computed = {
    "M_each_arm_N_m": float(submitted["M"]),
    "L_mm": float(1000 * submitted["L"]),
    "delta_at_L_over_2_um": float(1e6 * submitted["dhaf"]),
}
reported = {
    "M_each_arm_N_m": 0.2126,
    "L_mm": 10.80,
    "delta_at_L_over_2_um": 38.2,
}
more_digits = {
    "M_each_arm_N_m": 0.21263,
    "L_mm": 10.800,
    "delta_at_L_over_2_um": 38.22,
}
comparisons = {}
for name, target in key.items():
    comparisons[name] = {
        "reference": target,
        "computed": computed[name],
        "computed_relative_error_percent": 100 * abs(computed[name] / target - 1),
        "reported": reported[name],
        "reported_relative_error_percent": 100 * abs(reported[name] / target - 1),
        "parenthetical": more_digits[name],
        "parenthetical_relative_error_percent": 100 * abs(more_digits[name] / target - 1),
    }

characteristic_zero_front = lambda z: 1 + np.cosh(z) * np.cos(z)
zero_front_root = brentq(characteristic_zero_front, np.pi / 2, 2)
x_mid = submitted["L"] / 2
ode_rhs = submitted["lam"] ** 4 * submitted["u"](x_mid)
finite_difference_checks = []
for step in (1e-6, 2e-5, 1e-4, 2e-4):
    stencil = submitted["u"](x_mid + step * np.arange(-2, 3))
    fourth = (stencil[0] - 4 * stencil[1] + 6 * stencil[2]
              - 4 * stencil[3] + stencil[4]) / step**4
    finite_difference_checks.append({
        "step_m": step,
        "absolute_residual_m_minus_4": float(fourth - ode_rhs),
        "relative_residual": float(abs(fourth - ode_rhs) / abs(ode_rhs)),
    })

result = {
    "comparisons": comparisons,
    "energy_relative_error": float(abs(submitted["lhs"] / submitted["rhs"] - 1)),
    "statics_absolute_error_N": float(abs(submitted["m_ap"] - submitted["mtip"] - submitted["Ib"])),
    "correct_zero_front_root": float(zero_front_root),
    "characteristic_residual_at_claimed_pi_over_2": float(characteristic_zero_front(np.pi / 2)),
    "finite_difference_checks": finite_difference_checks,
    "starter_cut_full_relative_end_rotation_rad": float(
        submitted["dc"] * abs(submitted["up"](submitted["L"]))
        + submitted["m_ap"] / submitted["dbar"] * 0.020),
}
(folder / "scratch" / "audit_results.json").write_text(json.dumps(result, indent=2) + "\n")
print(json.dumps(result, indent=2))
