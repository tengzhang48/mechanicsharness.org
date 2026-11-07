import numpy as np
from scipy.linalg import expm
from solver2 import assemble, lobatto_grid

def build_P(N, kappa=1e6):
    x, V, Vinv, W, D1v, D2v, D3v, A, B = assemble(N)
    P = np.zeros((N+1, N+1))
    for i in range(1, N):
        P[i, :] = A[i, :] @ D3v        # evolution row: L[f] = K[f''']
    # pin BCs stiffly: rows 0 (x=1): f'(1)->0 ; N (x=0): f(0)->0 ; N-1: f''(1)->0
    P[0, :] = -kappa * D1v[0, :]
    P[N, :] = -kappa * np.eye(N+1)[N, :]
    P[N-1, :] = -kappa * D2v[0, :]
    return x, P

N = 240
x, P = build_P(N)
sig0 = np.ones(N+1)
i05 = np.argmin(np.abs(x - 0.5)); i025 = np.argmin(np.abs(x - 0.25))
ss = np.array([0.001, 0.003, 0.01, 0.02, 0.04, 0.07, 0.1, 0.15, 0.2, 0.3, 0.45, 0.6])
s05, s025 = [], []
for s in ss:
    v = expm(P*s) @ sig0
    s05.append(v[i05]); s025.append(v[i025])
s05 = np.array(s05); s025 = np.array(s025)
for s, a, b in zip(ss, s025, s05):
    print(f"s={s:6.3f}  sigma(0.25)={a:+.8e}  sigma(0.5)={b:+.8e}")
lam_eff = -np.diff(np.log(np.abs(s05)))/np.diff(ss)
print("lambda_eff(0.5) between consecutive s:", lam_eff)
