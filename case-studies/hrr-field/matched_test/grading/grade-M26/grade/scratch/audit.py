"""Grader audit, run from M26; do not modify the submitted calculation."""
from pathlib import Path
import contextlib
import io
import json
import runpy
import numpy as np

with contextlib.redirect_stdout(io.StringIO()):
    scope = runpy.run_path('scratch/solve_core.py')
reference = json.loads(Path('reference_output.json').read_text())['route_A_closed_reduction']
reported = {
    'T_MPa': 4816.47,
    'C_rr_60deg_sqrtmm': 3.9420e-3,
    'C_thetatheta_60deg_sqrtmm': 6.3318e-4,
    'C_engineering_rtheta_60deg_sqrtmm': 5.7311e-3,
    'C_y_90deg_sqrtmm': 9.1390e-3,
    'r_min_mm': 0.010000,
    'r_max_mm': 0.061325,
}
checks = {}
for key, value in reported.items():
    relative_error = abs(value / reference[key] - 1)
    tolerance = 0.01 if key == 'T_MPa' else 0.02
    checks[key] = {
        'reported': value,
        'reference': reference[key],
        'relative_error_percent': 100*relative_error,
        'tolerance_percent': 100*tolerance,
        'within_tolerance': relative_error <= tolerance,
    }
checks['r_max_coarsely_rounded'] = {
    'reported': 0.0613,
    'relative_error_percent': 100*abs(0.0613/reference['r_max_mm'] - 1),
    'within_tolerance': abs(0.0613/reference['r_max_mm'] - 1) <= 0.02,
}
print('Reported-value comparisons:')
print(json.dumps(checks, indent=2))

# Correcting this diagnostic here establishes the source of the error. It is
# not counted as a verification performed by the submitter.
T, nue, K = scope['T'], scope['nue'], scope['K']
u, f = scope['u_closed'], scope['f_cart']
audit = []
for r, theta in [(0.05, 0.6), (0.137, 2.1), (0.02, np.pi/2)]:
    x, y = r*np.cos(theta), r*np.sin(theta)
    h = 1e-6
    dx = (np.asarray(u(x+h,y)) - np.asarray(u(x-h,y)))/(2*h)
    dy = (np.asarray(u(x,y+h)) - np.asarray(u(x,y-h)))/(2*h)
    a, b, c = f(theta)
    sc = K/np.sqrt(2*np.pi*r)
    strain = np.array([
        [(a-nue*b)*sc/T, (1+nue)*c*sc/T],
        [(1+nue)*c*sc/T, (b-nue*a)*sc/T],
    ])
    grad = np.column_stack([dx,dy])
    reconstructed = (grad + grad.T)/2
    audit.append({
        'r_mm': r,
        'theta_rad': theta,
        'correct_max_strain_residual': float(np.max(np.abs(reconstructed-strain))),
        'submission_shear_diagnostic_residual': float(0.5*(dy[0]+dy[1])-strain[0,1]),
    })
print('Shear diagnostic audit:')
print(json.dumps(audit, indent=2))

# Czz(T)=-sqrt(J/(2pi))*(T^-1/2-b*T^1/2), b>0, is strictly
# increasing over the admissible interval, so the submitted root is unique.
b = (1-2*scope['nu'])/scope['E']
c = np.sqrt(scope['Jg']/(2*np.pi))
dczz = c*(0.5*T**(-1.5)+0.5*b*T**(-0.5))
print('Calibrated T is admissible:', 0 < T < scope['E'])
print('Calibration derivative dCzz/dT at solution:', dczz)
Path('scratch/audit_results.json').write_text(json.dumps({
    'numerical_comparisons': checks,
    'shear_derivative_audit': audit,
    'admissible_T': bool(0<T<scope['E']),
    'calibration_derivative_at_T': float(dczz),
    'reproduced_shear_compatibility_residual': float(scope['res_shear']),
    'reproduced_Cy_integrated': float(scope['Cy_num']),
    'reproduced_Cy_closed_form': float(scope['Cy_cf']),
}, indent=2)+'\n')
