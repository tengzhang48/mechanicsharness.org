"""Read-only audit of delivered values, captured checks, and copied reruns."""
import json
import math
from pathlib import Path

import numpy as np

root = Path(__file__).resolve().parent.parent
answer = (root / "submission/ANSWER.md").read_text()
delivered = json.loads((root / "submission/results.json").read_text())
printed = json.loads((root / "scratch/reproduction/solve.stdout.json").read_text())
rerun = [json.loads(line) for line in
         (root / "scratch/reproduction/tip_check.stdout").read_text().splitlines()]
captured = [json.loads(line) for line in answer.splitlines()
            if line.startswith('{"n":')]
matches = []
for actual in rerun:
    expected = next(x for x in captured
                    if x["n"] == actual["n"] and x["ell"] == actual["ell"])
    matches.append({
        "n": actual["n"], "ell": actual["ell"],
        "max_tip_difference": max(abs(actual["K_from_opening"][j] -
                                      expected["K_from_opening"][j])
                                  for j in ("2", "4", "8", "16")),
        "volume_difference": actual["volume"] - expected["volume"],
        "relative_residual": actual["relative_residual"],
    })

nu, nuu, alpha = .2, .4, .89
b = (nuu - nu) / (alpha * (1 - 2 * nu))
eta = alpha * (1 - 2 * nu) / (2 * (1 - nu))
Du = 1 / (1 - nuu)
k = 1.
m = math.sqrt(k*k + 1.)
v, q_y = .3, .4
# Directly eliminate p0 from ANSWER.md:62-65, with q_y = p_y.
p0 = (2*b*Du*k*(m-k)*v - q_y) / m
T_direct = -Du*k*v - 2*eta*k/(m+k)*(p0 - 2*b*Du*k*v)
A = -Du*k*(1 - 4*eta*b*k*k/(m*(m+k)))
B_code = -2*eta*k/(m*(m+k))
C = 2*b*Du*k*(m-k)/m
T_as_named = A*v + B_code*q_y
p_as_named = C*v + q_y/m
T_with_q_minus_py = A*v + B_code*(-q_y)
p_with_q_minus_py = C*v + (-q_y)/m

V50, V25, V125 = 1.9794399541751084, 1.9422132554242526, 1.9213397890115673
evidence = {
    "printed_equals_delivered": printed == delivered,
    "K1_0_relative_error_to_corrected_key": abs(delivered["K1_0"]/.807593 - 1),
    "V1_short_relative_error_to_corrected_key": abs(delivered["V1_short_coeff"]/2.440220 - 1),
    "V1_short_missing_coefficient": 2.440220 - delivered["V1_short_coeff"],
    "claimed_laplace_volume_slope": delivered["V1_short_coeff"]*math.sqrt(math.pi)/2,
    "corrected_key_laplace_volume_slope": 2.440220*math.sqrt(math.pi)/2,
    "fixed_grid_volume_slopes_n2048": [
        (V50-V25)/(.05-.025), (V25-V125)/(.025-.0125)],
    "volume_ratio": delivered["V1_inf"]/delivered["V1_0"],
    "finest_rerun_matches": matches,
    "matrix_flux_sign_audit": {
        "T_direct_from_stated_surface_relations": T_direct,
        "T_matrix_if_q_equals_py_as_written": T_as_named,
        "p0_direct_from_stated_surface_relations": p0,
        "p0_matrix_if_q_equals_py_as_written": p_as_named,
        "T_matrix_if_q_equals_minus_py": T_with_q_minus_py,
        "p0_matrix_if_q_equals_minus_py": p_with_q_minus_py,
    },
}

out = root / "scratch/evidence.json"
out.write_text(json.dumps(evidence, indent=2) + "\n")
print(json.dumps(evidence, indent=2))
