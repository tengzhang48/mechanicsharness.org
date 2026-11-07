"""Read-only audit of the delivered numerical checks and limiting claims."""
import ast
import json
from pathlib import Path

import numpy as np

tree = ast.parse(Path('submission/solve.py').read_text())
functions = [node for node in tree.body
             if isinstance(node, ast.FunctionDef) and node.name == 'sigma_kk']
namespace = {'np': np}
exec(compile(ast.Module(body=functions, type_ignores=[]), 'submission/solve.py', 'exec'), namespace)
reported_stress = namespace['sigma_kk']


def correct_stress(x, y):
    z = np.asarray(x) + 1j * y
    w = np.sqrt(z - 1.0) * np.sqrt(z + 1.0)
    return 2.0 * (np.real(z / w) - 1.0)


print('Reference-stress symmetry and decay:')
for x, y in [(0.5, 0.1), (-0.5, 0.1), (100.0, 0.1), (-100.0, 0.1)]:
    print(f'x={x:g}, y={y:g}: delivered={reported_stress(x,y):.12g}, '
          f'correct_branch={correct_stress(x,y):.12g}')

print('\nRefined face-integral claims at x endpoints +/-0.9999:')
xs = np.linspace(-0.9999, 0.9999, 200001)
for y in [1e-6, 1e-4, 1e-2]:
    delivered = np.trapezoid(reported_stress(xs, y), xs)
    fixed = np.trapezoid(correct_stress(xs, y), xs)
    z_end = 0.9999 + 1j*y
    exact_fixed = 4.0 * (np.real(np.sqrt(z_end-1.0)*np.sqrt(z_end+1.0)) - 0.9999)
    print(f'y={y:.0e}: delivered={delivered:.9f}, correct_branch={fixed:.9f}, '
          f'correct_analytic={exact_fixed:.9f}')

reported = json.loads(Path('submission/results.json').read_text())
nu, nu_u, alpha = 0.2, 0.4, 0.89
eta = alpha*(1-2*nu)/(2*(1-nu))
key = {'K1_0': 0.807593, 'K1_inf': 1.0,
       'V1_0': np.pi*(1-nu_u), 'V1_inf': np.pi*(1-nu),
       'V1_short_coeff': 2.440220,
       'p_faces_0': 2*(nu_u-nu)/(alpha*(1-2*nu)),
       'K2_0': 0.0, 'K2_inf': -eta,
       'V2_inf': -alpha*np.pi*(1-2*nu)/2}
print('\nCorrected-key numerical comparison:')
for name, target in key.items():
    value = reported[name]
    error = abs(value-target)
    tolerance = 0.01 if target == 0 else 0.01*abs(target)
    relative = 'n/a' if target == 0 else f'{100*error/abs(target):.6f}%'
    print(f'{name}: submission={value:.12g}, key={target:.12g}, '
          f'relative_error={relative}, in_band={error <= tolerance}')
print('Corrected faces contribution:', reported['V1_short_coeff']*(1-nu_u)/(1-nu))
print('Corrected tips contribution:', 2.440220-reported['V1_short_coeff']*(1-nu_u)/(1-nu))
print('Loading-2 normalized volume/K relation:', np.pi*(1-nu)*key['K2_inf'])
