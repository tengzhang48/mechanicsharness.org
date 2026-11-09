"""Local grading audit of the submitted shooting method and auxiliary claims."""
from pathlib import Path
import importlib.util
import json
import sys

import numpy as np
from scipy.optimize import least_squares

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE / "submission"))
import solve as submitted

spec = importlib.util.spec_from_file_location("notch_reference", HERE / "reference/reference.py")
reference = importlib.util.module_from_spec(spec)
spec.loader.exec_module(reference)

n, beta = 4.0, 2 * np.pi / 3
key_p, key_R = 0.156175260286, 0.373262313312
studies = []
for rtol, atol in [(1e-9, 1e-11), (1e-10, 1e-12), (1e-11, 1e-13),
                   (1e-12, 1e-14), (1e-13, 1e-15)]:
    fit = least_squares(
        lambda x: submitted.residual(x[0], x[1], n, beta, rtol, atol),
        [0.1561752603, -1.26307900], xtol=1e-15, ftol=1e-15, gtol=1e-15)
    p, curvature = fit.x
    R = ((2 - p) + curvature) / ((2 - p) * (1 - p))
    studies.append(dict(rtol=rtol, atol=atol, success=bool(fit.success),
                        p=float(p), curvature=float(curvature), R=float(R),
                        face_residual=float(np.max(np.abs(fit.fun)))))

p, curvature = studies[-1]["p"], studies[-1]["curvature"]
sol = submitted.shoot(p, curvature, n, beta, rtol=1e-13, atol=1e-15)
theta = np.linspace(0, beta, 10001)
y = sol.sol(theta)
D = y[2] + p * (2 - p) * y[0]
T = -(1 - p) * y[1]
Q2 = 0.75 * D**2 + 3 * T**2
rhs_errors = [abs(submitted.rhs_maker(p, n)(t, yi)[3]
                  - reference.explicit_fourth(t, yi, n, 2 - p)[3])
              for t, yi in zip(theta[::100], y.T[::100])]
stt = (2 - p) * (1 - p)
srr = (2 - p) + curvature

output = dict(
    displayed_number_errors={"p_0.1562": abs(0.1562 - key_p),
                             "R_0.3733": abs(0.3733 - key_R),
                             "p_0.15618": abs(0.15618 - key_p),
                             "R_0.37326": abs(0.37326 - key_R)},
    shooting_tolerance_study=studies,
    angular_profile={"min_f_interior": float(np.min(y[0, 1:-1])),
                     "min_Q2": float(np.min(Q2)),
                     "theta_at_min_Q2": float(theta[np.argmin(Q2)]),
                     "max_RHS_difference_from_reference": float(max(rhs_errors))},
    bisector_stresses={"sigma_theta_theta": stt, "sigma_rr": srr,
                      "sigma_zz": (stt + srr) / 2},
    strain_radial_power=-n * p,
    displacement_radial_power=1 - n * p,
    energy_density_radial_power=-(n + 1) * p,
    energy_integral_radial_power=2 - (n + 1) * p,
    analytic_elastic_wedge=reference.linear_exact(beta))
serialized = json.dumps(output, indent=2)
(HERE / "grading_checks.json").write_text(serialized + "\n")
print(serialized)
