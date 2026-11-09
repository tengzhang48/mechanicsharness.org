"""Local grading audit; evaluates submitted mechanics without its optional checker."""
import contextlib
import io
import json
import runpy
from pathlib import Path

import numpy as np

with contextlib.redirect_stdout(io.StringIO()):
    sub = runpy.run_path(str(Path(__file__).with_name('submission_core.py')))
reference = json.loads(Path('scratch/reference_rerun.json').read_text())['route_A_closed_reduction']

# Use the values identifiable in the submitted results table, including its rounding.
reported = {
    'T_MPa': 4816.47,
    'C_rr_60deg_sqrtmm': 0.0039420,
    'C_thetatheta_60deg_sqrtmm': 0.00063318,
    'C_engineering_rtheta_60deg_sqrtmm': 0.0057311,
    'C_y_90deg_sqrtmm': 0.0091390,
    'r_min_mm': 0.01000,
    'r_max_mm': 0.06133,
}
comparisons = {}
for key, value in reported.items():
    tolerance = 0.01 if key == 'T_MPa' else 0.02
    error = abs(value - reference[key]) / abs(reference[key])
    comparisons[key] = {
        'reported': value, 'reference': reference[key],
        'relative_error_percent': error * 100,
        'tolerance_percent': tolerance * 100, 'passes': bool(error <= tolerance),
    }

t, nu_eff, k = sub['T'], sub['nus'], sub['K']
strain_coefficients = k / (t * np.sqrt(2 * np.pi)) * np.array([
    1 - nu_eff, 1 - nu_eff, -2 * nu_eff,
])
j_checks = {str(r): float(sub['J_circle'](r, k, t, nu_eff))
            for r in (0.02, 0.05, 0.1, 0.2)}

full_law_tests = {}
for r in (0.01, 0.2):
    stresses = sub['sig_K'](r, 0.0, k)
    ezz = sub['full_strain'](*stresses, t)[3]
    coefficient = float(np.sqrt(r) * ezz)
    full_law_tests[str(r)] = {
        'sqrt_r_ezz_from_leading_stress': coefficient,
        'difference_from_leading_coefficient_percent':
            abs(coefficient + 0.005) / 0.005 * 100,
    }

assert all(item['passes'] for item in comparisons.values())
assert all(abs(value / 0.8 - 1) < 1e-8 for value in j_checks.values())
assert np.argmax(np.abs(strain_coefficients)) == 2
assert sub['r_lo'] <= min(sub['r_hi'], 0.2)
assert abs((1 - 2 * nu_eff) / t - (1 - 2 * 0.3) / 70000) < 1e-18

out = {
    'tolerance_comparisons': comparisons,
    'three_principal_strain_coefficients_on_theta_zero': strain_coefficients.tolist(),
    'direct_J_checks_N_per_mm': j_checks,
    'leading_stress_full_law_tests': full_law_tests,
    'scope': 'The full-law tests substitute leading stresses; they are not a finite-r boundary-value solution.',
}
Path('scratch/audit_results.json').write_text(json.dumps(out, indent=2) + '\n')
print(json.dumps(out, indent=2))
