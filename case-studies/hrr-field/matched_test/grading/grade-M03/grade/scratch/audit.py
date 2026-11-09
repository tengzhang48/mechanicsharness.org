import contextlib
import io
import json
import runpy

import numpy as np

with contextlib.redirect_stdout(io.StringIO()):
    submitted = runpy.run_path('solve.py')
    flux_code = runpy.run_path('verify_J.py')

with open('../reference_output.json') as handle:
    reference = json.load(handle)['route_A_closed_reduction']

reported = {
    'T_MPa': 4816.47,
    'C_rr_60deg_sqrtmm': 0.0039420,
    'C_thetatheta_60deg_sqrtmm': 0.00063318,
    'C_engineering_rtheta_60deg_sqrtmm': 0.0057311,
    'C_y_90deg_sqrtmm': 0.0091390,
    'r_min_mm': 0.01000,
    'r_max_mm': 0.061325,
}
tolerances = {key: 0.01 if key == 'T_MPa' else 0.02 for key in reported}
comparison = {}
for key, value in reported.items():
    error = abs(value - reference[key]) / abs(reference[key])
    comparison[key] = {
        'reported': value,
        'reference': reference[key],
        'relative_error_percent': 100 * error,
        'tolerance_percent': 100 * tolerances[key],
        'within_tolerance': error <= tolerances[key],
    }
    assert error <= tolerances[key]

# Exercise the submitted flux function with the parameters and radii claimed
# in the answer. This is a grader rerun, not an added submission verification.
direct_flux = {
    str(radius): float(flux_code['J_circle'](
        submitted['K'], radius, submitted['E_eff'], submitted['nu_eff']))
    for radius in (0.001, 0.015, 0.04, 0.05)
}

# Illustrate why the 10% bracket screen is not a total-strain error bound.
# Evaluate the complete law on the leading stress, not a finite-r BVP solution.
t = submitted['T']
sigma = np.diag([100.0, 100.0, 0.0])
trace = np.trace(sigma)
deviator = sigma - np.eye(3) * trace / 3.0
equivalent = np.sqrt(1.5 * np.sum(deviator ** 2))
elastic = 1.3 / 70000.0 * sigma - 0.3 / 70000.0 * trace * np.eye(3)
deviatoric = 1.5 * (1.0 / t - 1.0 / 70000.0) * deviator
leading_strain = elastic + deviatoric
complete_strain = elastic + max(1 - 10 / equivalent, 0) * deviatoric

output = {
    'reported_number_comparison': comparison,
    'submitted_flux_function_at_calibrated_parameters': direct_flux,
    'displacement_check': {
        'formula_C_y': float(submitted['C_y']),
        'strain_reconstruction_C_y': float(submitted['uy_factor']),
        'absolute_difference': float(abs(submitted['C_y'] - submitted['uy_factor'])),
    },
    'theta_zero_constitutive_comparison_at_sigma_e_100_MPa': {
        'leading_principal_strains': np.diag(leading_strain).tolist(),
        'complete_law_at_same_stress_principal_strains': np.diag(complete_strain).tolist(),
        'zz_difference_relative_to_complete_law_percent': float(
            100 * abs(leading_strain[2, 2] - complete_strain[2, 2])
            / abs(complete_strain[2, 2])),
    },
}
print(json.dumps(output, indent=2))
