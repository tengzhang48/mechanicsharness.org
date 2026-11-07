"""Independent grading diagnostics; does not change submitted files."""
import ast
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
tree = ast.parse((ROOT / 'submission/verify.py').read_text())
fn = next(node for node in tree.body if isinstance(node, ast.FunctionDef) and node.name == 'crack_solve')
env = {'np': np, 'G': 1.0, 'L': 1.0, 'nu': 0.2, 'kappa': 2.2}
exec(compile(ast.Module(body=[fn], type_ignores=[]), 'submission/verify.py', 'exec'), env)
solver = env['crack_solve']

Nc, Nm = 241, 160
theta = (np.arange(Nc) + 0.5) * np.pi / Nc
q = 1 + 0.5 * (4 * np.cos(theta)**2 - 1)
th, x, a, w, reported_V, K, n = solver(q, Nc, Nm)
true_V = -np.pi * env['L'] * a[0] / 2
th_mid = (np.arange(10000) + 0.5) * np.pi / 10000
integrated_V = np.sum(w(th_mid) * env['L'] * np.sin(th_mid)) * np.pi / 10000
env['L'] = 2.0
_, _, _, _, _, K_L2, _ = solver(np.ones(Nc), Nc, Nm)

nu, nu_u, alpha, G, permeability = .2, .4, .89, 1.0, 1.0
gamma = alpha * (1 - 2 * nu) / (1 - nu)
storage_H = alpha**2 * (1 - 2 * nu)**2 * (1 - nu_u) / ((nu_u - nu) * (1 - nu))
c_from_balance = 2 * G * permeability / storage_H
c_given = 2 * permeability * G * (1 - nu) * (nu_u - nu) / (alpha**2 * (1 - 2 * nu)**2 * (1 - nu_u))

dt = 5e-6
neighbor = 1.0
erroneous_neumann_boundary = -dt * neighbor / (1 - dt)
h = 8 / 2048

report = {
    'manufactured_load': 'q(x)=1+0.5*U_2(x/L)',
    'coded_volume': float(reported_V),
    'volume_from_opening_integral': float(integrated_V),
    'correct_volume_from_a1': float(true_V),
    'uniform_load_K_L2_coded_normalized': float(K_L2 / np.sqrt(np.pi * 2)),
    'uniform_load_K_L2_expected_normalized': 1.0,
    'neumann_bc_neighbor': neighbor,
    'neumann_bc_code_imposes_boundary': erroneous_neumann_boundary,
    'neumann_bc_expected_boundary': neighbor,
    'mirroring_negative_x_minus_h_maps_to_x': 0.0,
    'mirroring_negative_x_minus_h_should_map_to_x': h,
    'strip_p_y_at_x2L_L1': float(2 / (3 * np.pi)),
    'uncracked_plane_p_y_expected': 0.0,
    'c_conditional_from_balance': c_from_balance,
    'c_given': c_given,
    'literal_3D_extension_diffusivity_at_given_parameters': 2 * permeability * G * (nu_u - nu) / (alpha**2 * (1 - 2 * nu) * (1 - 2 * nu_u)),
    'general_pressure_equation': 'p_t = c*laplacian(p) - alpha*(1-2*nu)/storage_H * H_t, H=s+gamma*p',
}
fft_convergence = []
for box, count in [(4, 256), (4, 512), (4, 1024), (4, 2048), (8, 1024), (16, 2048)]:
    xx = (np.arange(count) - count // 2) * 2 * box / count
    X, Y = np.meshgrid(xx, xx, indexing='xy')
    P = (np.arctan2(X + 1, np.abs(Y)) + np.arctan2(1 - X, np.abs(Y))) / np.pi
    k = 2 * np.pi * np.fft.fftfreq(count, d=2 * box / count)
    KX, KY = np.meshgrid(k, k, indexing='xy')
    K2 = KX**2 + KY**2
    kernel = np.zeros_like(K2)
    np.divide(KX**2, K2, out=kernel, where=K2 > 0)
    stress = np.fft.ifft2(-gamma * kernel * np.fft.fft2(P)).real
    sample = float(stress[count // 2, np.argmin(np.abs(xx - 0.5))])
    fft_convergence.append({'box': box, 'N': count, 'trace_at_x05': sample})
report['fft_convergence'] = fft_convergence
(ROOT / 'scratch/audit_results.json').write_text(json.dumps(report, indent=2) + '\n')
print(json.dumps(report, indent=2))
