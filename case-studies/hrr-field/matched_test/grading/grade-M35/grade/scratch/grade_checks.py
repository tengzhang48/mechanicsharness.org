import ast
import json
import math
from pathlib import Path

import numpy as np

root = Path(__file__).resolve().parent.parent
ref = json.loads((root / 'reference_output.json').read_text())['route_A_closed_reduction']
reported = {
    'T_MPa': 4816.5,
    'C_rr_60deg_sqrtmm': 0.0039420,
    'C_thetatheta_60deg_sqrtmm': 0.00063318,
    'C_engineering_rtheta_60deg_sqrtmm': 0.0057311,
    'C_y_90deg_sqrtmm': 0.0091390,
    'r_min_mm': 0.01000,
    'r_max_mm': 0.06133,
}
checks = {}
for key, value in reported.items():
    err = abs(value - ref[key]) / abs(ref[key])
    tol = 0.01 if key == 'T_MPa' else 0.02
    checks[key] = dict(reported=value, reference=ref[key], relative_error_percent=100*err,
                       tolerance_percent=100*tol, passes=bool(err <= tol))

# Load only the functions and imports from the submitted code, skipping its
# expensive top-level reruns. The constants are the already-reproduced solution.
def load_functions(filename, corrected_derivative=False):
    src = (root / 'submission' / filename).read_text()
    if corrected_derivative:
        src = src.replace('uypx = (uyp[1] - uym[1]) / (2 * h)',
                          'uypx = (uxp[1] - uxm[1]) / (2 * h)')
    tree = ast.parse(src)
    tree.body = [node for node in tree.body
                 if isinstance(node, (ast.Import, ast.ImportFrom, ast.FunctionDef))]
    env = dict(E=70000., nu=0.30, sy=10., J=0.800, Czz_meas=-0.005,
               Tsol=ref['T_MPa'], Et=ref['T_MPa'], nt=ref['nu_effective'],
               K=math.sqrt(0.8*ref['T_MPa']))
    env['phi0'] = env['K'] / (env['Et']*math.sqrt(2*math.pi))
    env['trapz'] = getattr(np, 'trapezoid', None) or np.trapz
    exec(compile(tree, str(root / 'submission' / filename), 'exec'), env)
    return env

s = load_functions('solve.py')
v = load_functions('verify3.py')
s_fixed = load_functions('solve.py', corrected_derivative=True)
angles = np.linspace(-math.pi, math.pi, 401)
old_displacement = np.asarray(s['lead_disp'](0.003, angles, ref['T_MPa']))
new_displacement = v['phi0']*math.sqrt(0.003)*np.asarray(v['Fcart'](angles))
identity_error = float(np.max(np.abs(old_displacement-new_displacement)))
strain_errors = []
for th in np.radians([0, 30, 45, 60, 90, 120]):
    F = np.array(v['Fcart'](th))
    Fp = np.array(v['Fcart_d'](th))
    scale = v['phi0']/math.sqrt(0.003)
    grad = scale*np.column_stack((0.5*math.cos(th)*F-math.sin(th)*Fp,
                                 0.5*math.sin(th)*F+math.cos(th)*Fp))
    strain = (grad+grad.T)/2
    g = v['ang'](th)
    law = scale*np.array([[g['gxx'],g['gxy']],[g['gxy'],g['gyy']]])
    strain_errors.append(float(np.max(np.abs(strain-law))))

bad_J = float(s['J_contour'](0.01, ref['T_MPa'], 'lin'))
fixed_J = float(s_fixed['J_contour'](0.01, ref['T_MPa'], 'lin'))
leading_J = float(v['J_contour'](0.01, 'lin'))
full_J = {str(r): float(v['J_contour'](r, 'full')) for r in (0.1, 0.01, 0.001)}
error_ratio = (0.8-full_J['0.1'])/(0.8-full_J['0.01'])

out = dict(tolerance_checks=checks,
           max_abs_old_vs_reconstructed_displacement=identity_error,
           max_abs_verify3_strain_error=float(max(strain_errors)),
           solve_original_J_N_per_mm=bad_J,
           solve_J_with_only_derivative_corrected_N_per_mm=fixed_J,
           verify3_leading_J_N_per_mm=leading_J,
           verify3_full_energy_with_leading_field_J_N_per_mm=full_J,
           full_energy_J_error_ratio_for_tenfold_radius_decrease=float(error_ratio),
           inferred_error_power=float(math.log10(error_ratio)))
(root/'scratch'/'grade_checks.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
assert all(c['passes'] for c in checks.values())
assert identity_error < 1e-15
assert max(strain_errors) < 1e-9
assert abs(fixed_J - 0.8) < 1e-8
assert abs(leading_J - 0.8) < 1e-8
