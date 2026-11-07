import numpy as np
from scipy.linalg import expm
from solver2 import assemble, dcoef_matrix, cheb_V
import warnings; warnings.filterwarnings('ignore')

def build_P(N, kappa=1e4):
    x, V, Vinv, W, D1v, D2v, D3v, A, B = assemble(N)
    D = dcoef_matrix(N)
    P = np.zeros((N+1, N+1))
    P[1:N, :] = A[1:N, :]          # interior rows in coefficient space
    P[0, :] = -kappa * (V @ D)[0, :]       # pin f'(1)=0
    P[N, :] = -kappa * V[N, :]             # pin f(0)=0
    P[N-1, :] = -kappa * (V @ D @ D)[0, :] # pin f''(1)=0
    return x, V, P

N = 200
x, V, P = build_P(N)
c0 = np.linalg.solve(V, np.ones(N+1))   # project sigma0=1
ss = np.array([0.002, 0.005, 0.01, 0.02, 0.035, 0.05, 0.08, 0.12, 0.18, 0.25, 0.35, 0.5])
i05 = np.argmin(np.abs(x-0.5)); i03 = np.argmin(np.abs(x-0.3)); i08 = np.argmin(np.abs(x-0.8))
rows = []
for s in ss:
    c = expm(P*s) @ c0
    v = V @ c
    rows.append((s, v[i03], v[i05], v[i08]))
for r in rows:
    print(f"s={r[0]:6.3f}  sig(0.3)={r[1]:+.8e}  sig(0.5)={r[2]:+.8e}  sig(0.8)={r[3]:+.8e}")
s_ = np.array([r[0] for r in rows]); y = np.array([r[2] for r in rows])
print("lam_eff(0.5):", -np.diff(np.log(np.abs(y)))/np.diff(s_))
