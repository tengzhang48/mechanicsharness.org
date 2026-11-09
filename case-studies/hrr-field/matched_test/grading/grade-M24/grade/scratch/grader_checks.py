"""Reproducible grading checks; original submission and reference are unchanged."""

import contextlib
import io
import json
import runpy

import numpy as np

with open('reference_run.json') as handle:
    reference = json.load(handle)['route_A_closed_reduction']

reported = {
    'T_MPa': 4816.0,
    'C_rr_60deg_sqrtmm': 0.003942,
    'C_thetatheta_60deg_sqrtmm': 0.0006332,
    'C_engineering_rtheta_60deg_sqrtmm': 0.005731,
    'C_y_90deg_sqrtmm': 0.009139,
    'r_min_mm': 0.0100,
    'r_max_mm': 0.0613,
}
tolerances = {key: 0.01 if key == 'T_MPa' else 0.02 for key in reported}
comparisons = {
    key: {
        'reported': value,
        'reference': reference[key],
        'relative_error_percent': 100 * abs(value / reference[key] - 1),
        'allowed_percent': 100 * tolerances[key],
        'within_tolerance': abs(value / reference[key] - 1) <= tolerances[key],
    }
    for key, value in reported.items()
}
assert all(row['within_tolerance'] for row in comparisons.values())

with contextlib.redirect_stdout(io.StringIO()):
    submission = runpy.run_path('submission_solve.py')

t, nut, k = (submission[key] for key in ('T', 'nup', 'K'))
actual_outputs = {
    'T_MPa': t,
    'C_rr_60deg_sqrtmm': submission['err'],
    'C_thetatheta_60deg_sqrtmm': submission['ett'],
    'C_engineering_rtheta_60deg_sqrtmm': 2 * submission['ert'],
    'C_y_90deg_sqrtmm': submission['uy'],
    'r_min_mm': submission['r_lo'],
    'r_max_mm': submission['r_hi'],
}
code_relative_errors = {
    key: abs(value / reference[key] - 1) for key, value in actual_outputs.items()
}
assert max(code_relative_errors.values()) < 1e-9

js = {str(r): float(submission['J_integral'](r)) for r in (0.001, 0.01, 0.1)}
assert all(abs(value / 0.8 - 1) < 1e-8 for value in js.values())

e, nu, sy = 70000.0, 0.30, 10.0
a = 1 / t - 1 / e
bulk_original = e / (3 * (1 - 2 * nu))
bulk_reduced = t / (3 * (1 - 2 * nut))
assert abs(bulk_reduced / bulk_original - 1) < 1e-12

ahead_coefficients = k / (np.sqrt(2 * np.pi) * t) * np.array([1 - nut, 1 - nut, -2 * nut])
principal_at_lower = ahead_coefficients / np.sqrt(submission['r_lo'])
assert np.argmax(abs(principal_at_lower)) == 2
assert abs(np.max(abs(principal_at_lower)) - 0.05) < 1e-10

# The displayed calibration equation in ANSWER.md line 30 has sqrt(T), not T,
# in its denominator. Compare its literal evaluation with the correct equation.
literal_czz_equation = -2 * nut * np.sqrt(0.8 * t) / (np.sqrt(2 * np.pi) * np.sqrt(t))
correct_czz_equation = -2 * nut * np.sqrt(0.8 * t) / (np.sqrt(2 * np.pi) * t)

# Complementary energy Phi(sigma) generates strain. Work W equals sigma:epsilon-Phi.
uniaxial_stress = 30.0
phi = uniaxial_stress**2 / (2 * e) + a / 2 * (uniaxial_stress - sy)**2
strain = uniaxial_stress / e + a * (uniaxial_stress - sy)
work = uniaxial_stress * strain - phi

out = {
    'reported_output_tolerances': comparisons,
    'maximum_code_output_relative_error': max(code_relative_errors.values()),
    'submission_J_integral_N_per_mm': js,
    'bulk_modulus_original_MPa': bulk_original,
    'bulk_modulus_reduced_MPa': bulk_reduced,
    'all_three_principal_strains_at_lower_endpoint': principal_at_lower.tolist(),
    'literal_displayed_Czz_equation': literal_czz_equation,
    'correct_Czz_equation': correct_czz_equation,
    'complementary_energy_at_30MPa': phi,
    'strain_work_at_30MPa': work,
}
with open('grader_checks.json', 'w') as handle:
    json.dump(out, handle, indent=2)
print(json.dumps(out, indent=2))
