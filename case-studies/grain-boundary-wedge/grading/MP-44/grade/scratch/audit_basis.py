"""Grader audit; imports only the copied submission and changes no submitted files."""
import json
import os
from pathlib import Path
import sys

import numpy as np
from scipy.linalg import eig

sys.path.insert(0, str(Path(__file__).parent / "reproduction"))
import solve as submitted

print("Recorded basis-audit PID:", os.getpid(), flush=True)
N = 64
Zc = 0.5 * (1 - np.cos(np.arange(1, N + 1) * np.pi / (N + 1)))
R = np.array([[submitted.Rf(k, z) for k in range(N)] for z in Zc])
M = np.array([[submitted.fk(k, z) for k in range(N)] for z in Zc])
parity = (-1.0) ** np.arange(N)

def select(lam):
    return sorted((v.real, i) for i, v in enumerate(lam)
                  if v.real > 1e-6 and abs(v.imag) < 1e-6 * max(1, abs(v.real)))

lam, vec = eig(R, M)
candidates = select(lam)
lambda_bad, i = candidates[0]
c = vec[:, i].real
z_test = np.linspace(0.05, 0.95, 19)
Mt = np.array([[submitted.fk(k, z) for k in range(N)] for z in z_test])
Rt = np.array([[submitted.Rf(k, z) for k in range(N)] for z in z_test])
f_bad = Mt @ c
same_matrix_ratio = (Rt @ c) / f_bad
physical_ratio = (Rt @ (parity * c)) / f_bad
print("Submitted eigenvalues (first 3 positive real):", [v for v, _ in candidates[:3]])
print("Negative-real eigenvalues:", sorted(v.real for v in lam if abs(v.imag) < 1e-6 * max(1, abs(v.real)) and v.real < 0))
print("Complex eigenvalues with positive real part below reported lambda1:",
      [str(v) for v in lam if 0 < v.real < lambda_bad and abs(v.imag) >= 1e-6 * max(1, abs(v.real))])
print("Points:", z_test.tolist())
print("Submitted residual max:", float(np.max(np.abs(same_matrix_ratio - lambda_bad))))
print("Physical R[f]/f (differentiate submitted f correctly):", physical_ratio.tolist())
print("Physical residual max:", float(np.max(np.abs(physical_ratio - lambda_bad))))

# Since submitted.fk(k,z)''' = (-1)^k * g_k(z), flip only odd M columns.
fixed_lam, fixed_vec = eig(R, M * parity)
fixed_candidates = select(fixed_lam)
lambda_fixed, fi = fixed_candidates[0]
fc = fixed_vec[:, fi].real
f_at_one = np.array([submitted.fk(k, 1.0) for k in range(N)]) * parity
f_at_half = np.array([submitted.fk(k, 0.5) for k in range(N)]) * parity
corrected = {
    "lambda1": lambda_fixed,
    "lambda2": fixed_candidates[1][0],
    "lambda3": fixed_candidates[2][0],
    "mode1_shape_half": float((f_at_half @ fc) / (f_at_one @ fc)),
    "relaxation_time_t0": 1 / lambda_fixed,
}
print("Corrected odd-column signs:", json.dumps(corrected, indent=2))
print("Basis k=1 at z=0.25: claimed g=", submitted.chebT(1, -0.5) / np.sqrt(0.25 * 0.75),
      "; actual fk'''=", -submitted.chebT(1, -0.5) / np.sqrt(0.25 * 0.75))
