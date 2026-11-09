"""Read-only mechanics/tolerance audit of the copied submission."""
import json
from pathlib import Path

import numpy as np
import solve as S
from scipy.optimize import brentq


reference = json.loads(Path('../reference_results.json').read_text())
key = reference['answer']
M, _ = S.M_from_tip()
L = S.L_for_M(M)
delta, coef = S.solve_delta(L, M)
computed = {
    'M_each_arm_N_m': float(M),
    'L_mm': float(L * 1e3),
    'delta_at_L_over_2_um': float(delta(L / 2)[0] * 1e6),
}
reported = dict(zip(key, [0.2126258, 10.80010, 38.21890]))
table_reported = dict(zip(key, [0.2126, 10.80, 38.22]))
errors = {
    source: {name: abs(values[name] / key[name] - 1) * 100 for name in key}
    for source, values in [('reported', reported), ('table_reported', table_reported),
                           ('rerun', computed)]
}

x = np.linspace(0, L, 1001)
A, B, C, D = coef[:4]
curvature = S.k**2 * (A * np.exp(S.k * x) + B * np.exp(-S.k * x)
                     - C * np.cos(S.k * x) - D * np.sin(S.k * x))
slope = S.k * (A * np.exp(S.k * x) - B * np.exp(-S.k * x)
               - C * np.sin(S.k * x) + D * np.cos(S.k * x))
shear = S.k**3 * (A * np.exp(S.k * x) - B * np.exp(-S.k * x)
                 + C * np.sin(S.k * x) - D * np.cos(S.k * x))
gap = delta(x)
work = S.sb * gap - S.sb * gap**2 / (2 * S.dc)
invariant = S.Dtil / S.b * (0.5 * curvature**2 - slope * shear) - work

a0_original = S.a0
a0_checks = []
for a0 in [0.005, 0.020, 0.100]:
    S.a0 = a0
    length = S.L_for_M(M)
    profile, _ = S.solve_delta(length, M)
    a0_checks.append({'a0_mm': a0 * 1e3, 'L_mm': length * 1e3,
                      'delta_mid_um': float(profile(length / 2)[0] * 1e6)})
S.a0 = a0_original

# A narrowly bracketed shooting root checks the first event, in contrast to
# verify.py's bracket extending to 50 mm and spanning nonphysical branches.
import verify as V
L_shoot = brentq(lambda length: float(V.fd_delta_at(length, M, [length])[0]) - S.dc,
                0.0105, 0.0110, xtol=1e-12)

result = {
    'reference': key,
    'submission_reported': reported,
    'submission_rerun': computed,
    'relative_error_percent': errors,
    'profile_diagnostics': {
        'delta_min_um': float(gap.min() * 1e6),
        'delta_max_um': float(gap.max() * 1e6),
        'slope_min_rad': float(slope.min()),
        'slope_max_rad': float(slope.max()),
        'curvature_at_front_inverse_m': float(curvature[0]),
        'required_curvature_at_front_inverse_m': float(np.sqrt(2 * S.b * S.Gam0 / S.Dtil)),
        'front_energy_J_m2': float(S.Dtil * curvature[0]**2 / (2 * S.b)),
        'rear_shear_third_derivative': float(shear[-1]),
        'rear_moment_N_m': float(S.Dtil * curvature[-1]),
        'maximum_first_integral_error_J_m2': float(np.max(np.abs(invariant - S.Gam0))),
    },
    'a0_sensitivity': a0_checks,
    'first_event_shooting_length_mm': float(L_shoot * 1e3),
}
Path('audit_results.json').write_text(json.dumps(result, indent=2) + '\n')
print(json.dumps(result, indent=2))
