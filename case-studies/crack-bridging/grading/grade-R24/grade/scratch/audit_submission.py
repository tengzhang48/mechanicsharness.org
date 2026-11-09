"""Grade supplied routes A/B without the unavailable spectral dependency.

Only existing imports, constants, function bodies, and root-finding statements
are extracted. No replacement spectral solver or mechanical formula is added.
Extra diagnostics below are grader checks, not claimed submission work.
"""
import ast
import json
from pathlib import Path

folder = Path(__file__).resolve().parent
source = (folder / 'solve.py').read_text()
tree = ast.parse(source)
constant_names = {'E', 'nu', 'Ep', 'a', 't', 'dc', 'G0', 'K0', 'LA', 'LB'}
function_names = {'Sigma_of_c', 'delta_t_A', 'delta_t_B'}
selected = []
for node in tree.body:
    if isinstance(node, ast.Import) and any(n.name == 'numpy' for n in node.names):
        selected.append(node)
    elif isinstance(node, ast.ImportFrom) and node.module in {'scipy.optimize', 'scipy.integrate'}:
        selected.append(node)
    elif isinstance(node, ast.Assign) and all(isinstance(t, ast.Name) and t.id in constant_names for t in node.targets):
        selected.append(node)
    elif isinstance(node, ast.FunctionDef) and node.name in function_names:
        selected.append(node)
extracted = '\n\n'.join(ast.get_source_segment(source, node) for node in selected) + '\n'
(folder / 'extracted_routes.py').write_text(extracted)
ns = {}
exec(compile(ast.Module(body=selected, type_ignores=[]), str(folder / 'solve.py'), 'exec'), ns)
np = ns['np']
Ep, a, t, dc, K0 = (ns[k] for k in ('Ep', 'a', 't', 'dc', 'K0'))
LA, LB = ns['LA'], ns['LB']
Sigma = ns['Sigma_of_c']
dA, dB = ns['delta_t_A'], ns['delta_t_B']
c = a + LA
stress = Sigma(c)
reference = json.loads((folder.parent / 'reference_results.json').read_text())
answer = reference['answer']
reported = {
    'stress_133.73_MPa': {'value': 133.73, 'target': answer['Sigma_MPa']},
    'stress_133.730_MPa': {'value': 133.730, 'target': answer['Sigma_MPa']},
    'length_4.086_mm': {'value': 4.086, 'target': answer['L_mm_per_end']},
    'length_4.0856_mm': {'value': 4.0856, 'target': answer['L_mm_per_end']}
}
for item in reported.values():
    item['absolute_relative_error_percent'] = 100 * abs(item['value'] / item['target'] - 1)
    item['within_0.3_percent'] = item['absolute_relative_error_percent'] <= 0.3
branch = []
for length in (0.5e-3, 2e-3, LA, 8e-3):
    branch.append({'L_mm': length * 1e3, 'opening_A_um': dA(length) * 1e6,
                   'opening_B_um': dB(length) * 1e6,
                   'relative_opening_difference': abs(dA(length) / dB(length) - 1)})
ell = Ep * dc / t
alpha = a / ell
S47 = t * (2 / np.pi) * np.arccos(np.exp(-np.pi / (8 * alpha)))
L47 = a * (np.exp(np.pi * Ep * dc / (8 * t * a)) - 1)
small_L = 1e-5
K_bridge = np.sqrt(np.pi * c) * (2 * t / np.pi) * np.arccos(a / c)
k = K0 / (t * np.sqrt(np.pi * a))
u_end = np.log(c / a)
us = np.linspace(1e-9, u_end, 10001)
slopes = 2 / (np.pi * np.sqrt(np.expm1(2 * us))) - 0.5 * k * np.exp(-us / 2)
result = {
    'extraction': 'Unmodified supplied A/B definitions and root statements; no spectral code',
    'reported_answer_errors': reported,
    'route_A': {'Sigma_MPa': stress / 1e6, 'L_mm_per_end': LA * 1e3, 'opening_um': dA(LA) * 1e6},
    'route_B': {'Sigma_MPa': Sigma(a + LB) / 1e6, 'L_mm_per_end': LB * 1e3, 'opening_um': dB(LB) * 1e6},
    'A_B_relative_L_difference': abs(LA / LB - 1),
    'A_B_relative_stress_difference': abs(stress / Sigma(a + LB) - 1),
    'A_B_along_branch': branch,
    'tip_residual_relative': abs((stress * np.sqrt(np.pi * c) - K_bridge) / K0 - 1),
    'minimum_sampled_load_slope_before_event': float(slopes.min()),
    'event_load_slope': float(slopes[-1]),
    'zero_tip_limit': {'Sigma_MPa': S47 / 1e6, 'L_mm_per_end': L47 * 1e3},
    'decomposition': {'K_part_stress_MPa': K0 / np.sqrt(np.pi * c) / 1e6,
                      'bridging_part_stress_MPa': (2 * t / np.pi) * np.arccos(a / c) / 1e6,
                      'K_part_opening_um': (4 * K0 / (Ep * np.sqrt(np.pi * c))) * np.sqrt(c*c-a*a) * 1e6,
                      'bridging_part_opening_um': 8 * t * a / (np.pi * Ep) * np.log(c/a) * 1e6,
                      'K_applied_MPa_sqrt_m': stress * np.sqrt(np.pi*c) / 1e6,
                      'K_bridging_MPa_sqrt_m': K_bridge / 1e6},
    'small_zone_check_as_supplied': {
        'L_over_a': small_L / a,
        'bridge_K_ratio': (np.sqrt(np.pi*(a+small_L))*(2*t/np.pi)*np.arccos(a/(a+small_L)))/(t*np.sqrt(8*small_L/np.pi)),
        'printed_opening_ratio_claimed_to_approach_1': dA(small_L)*Ep/(8*t*small_L),
        'zero_tip_correct_opening_ratio_to_8tL_over_piEp': a*np.log1p(small_L/a)/small_L
    }
}
(folder / 'audit_results.json').write_text(json.dumps(result, indent=2) + '\n')
print(json.dumps(result, indent=2))
