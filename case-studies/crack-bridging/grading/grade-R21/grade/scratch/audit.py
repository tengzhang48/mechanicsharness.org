import contextlib
import io
import json
from pathlib import Path
import runpy

import numpy as np
from scipy.optimize import brentq

BASE = Path(__file__).resolve().parent.parent
with contextlib.redirect_stdout(io.StringIO()):
    sub = runpy.run_path(str(BASE / 'scratch/submission/verify.py'))
ref = runpy.run_path(str(BASE / 'scratch/reference/reference.py'))
key = json.loads((BASE / 'reference_results.json').read_text())['answer']

names = ['M_each_arm_N_m', 'L_mm', 'delta_at_L_over_2_um']
table = [0.2126, 10.80, 38.22]
extra_digits = [0.21262579, 10.8001, 38.2189]
rerun = [sub['M'], 1e3 * sub['ell'], 1e6 * sub['d'](sub['ell'] / 2)]
errors = {}
for i, name in enumerate(names):
    errors[name] = {
        'reference': key[name],
        'displayed_value': table[i],
        'displayed_relative_error_percent': abs(table[i] / key[name] - 1) * 100,
        'extra_digits_value': extra_digits[i],
        'extra_digits_relative_error_percent': abs(extra_digits[i] / key[name] - 1) * 100,
        'rerun_value': rerun[i],
        'rerun_relative_error_percent': abs(rerun[i] / key[name] - 1) * 100,
    }

ell = sub['ell']
C = sub['b'] * sub['R']
grid = np.linspace(0, ell, 2001)
gamma = sub['G0'] / (sub['sig'] * sub['dc'])
exact_data, exact_profile = ref['exact_softening_solution'](gamma)
mirrored_profile = sub['d'](ell - grid)
reference_profile = sub['dc'] * exact_profile(grid / ref['ELL'])
actual_front_curvature = float(sub['d2'](ell))
required_front_curvature = float(np.sqrt(2 * C * sub['G0']))
implied_front_work = actual_front_curvature**2 / (2 * C)
first_integral = (0.5 * sub['d2'](grid)**2 - sub['d1'](grid) * sub['d3'](grid)) / C
connector_work_to_gap = sub['sig'] * (sub['d'](grid) - sub['d'](grid)**2 / (2 * sub['dc']))
local_energy_error = first_integral - sub['G0'] - connector_work_to_gap

branches = []
scan = np.linspace(1e-4, 0.06, 2000)
indices = np.where(np.diff(np.sign([sub['resid'](p) for p in scan])))[0]
for i in indices:
    length = brentq(sub['resid'], scan[i], scan[i + 1])
    u = sub['beta'] * length
    B = -sub['A'] * (np.sinh(u) + np.sin(u)) / (np.cosh(u) + np.cos(u))
    x = np.linspace(0, length, 2001)
    z = sub['beta'] * x
    d = sub['dc'] + sub['A'] * (np.cosh(z) - np.cos(z)) + B * (np.sinh(z) + np.sin(z))
    front_curvature = sub['beta']**2 * (sub['A'] * (np.cosh(u) + np.cos(u)) + B * (np.sinh(u) - np.sin(u)))
    branches.append({'L_mm': float(length * 1e3), 'minimum_gap_um': float(d.min() * 1e6),
                     'maximum_gap_um': float(d.max() * 1e6),
                     'bonded_front_curvature_inverse_m': float(front_curvature)})

output = {
    'numeric_comparison': errors,
    'coordinate_equivalence': {
        'submission_coordinate': 's = L - x, from starter-cut root toward bonded front',
        'actual_front_curvature_inverse_m': actual_front_curvature,
        'required_front_curvature_inverse_m': required_front_curvature,
        'implied_intrinsic_work_J_m2': implied_front_work,
        'maximum_mirrored_profile_error_um': float(np.max(np.abs(mirrored_profile - reference_profile)) * 1e6),
        'minimum_gap_um': float(sub['d'](grid).min() * 1e6),
        'maximum_gap_um': float(sub['d'](grid).max() * 1e6),
        'minimum_slope_in_problem_coordinate': float(np.min(-sub['d1'](grid))),
        'maximum_slope_in_problem_coordinate': float(np.max(-sub['d1'](grid))),
        'maximum_local_energy_residual_J_m2': float(np.max(np.abs(local_energy_error))),
    },
    'scalar_equation_roots_in_scan': branches,
    'claimed_scalar_residuals_actual': {'at_8mm': float(sub['resid'](0.008)),
                                      'at_12mm': float(sub['resid'](0.012))},
}
out = BASE / 'scratch/audit_results.json'
out.write_text(json.dumps(output, indent=2) + '\n')
print(json.dumps(output, indent=2))
