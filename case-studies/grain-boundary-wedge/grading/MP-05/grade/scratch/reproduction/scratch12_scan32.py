import numpy as np
from solver7 import RhoOps32
from solver5 import cheb_T
import warnings; warnings.filterwarnings('ignore')

def smin32(ops, lam, mult=3, x0c=1e-7):
    N = ops.N
    Nc = mult*(N+1)
    ii = np.arange(1, Nc+1)
    xc = (1.0 - np.cos((2*ii-1)*np.pi/(2*Nc)))/2.0
    hp = ops.h_pieces(xc, kmax=3)[3]
    rho = (1.0-xc)[:, None]**1.5 * cheb_T(xc, N)
    M = np.zeros((Nc+3, N+1))
    M[:Nc, :] = hp + lam*rho
    M[Nc, :]   = ops.h_pieces(np.array([x0c]), 0)[0]
    h1 = ops.h_pieces(np.array([1.0]), 2)
    M[Nc+1, :] = h1[1, 0]
    M[Nc+2, :] = h1[2, 0]
    return np.linalg.svd(M, compute_uv=False)

N = 24
ops = RhoOps32(N)
lams = np.logspace(np.log10(3), np.log10(5000), 80)
vals = np.array([smin32(ops, l)[-1] for l in lams])
# print local minima
for i in range(1, len(lams)-1):
    if vals[i] < vals[i-1] and vals[i] < vals[i+1]:
        print(f"local min: lam={lams[i]:.4f} smin={vals[i]:.3e}")
print("global min:", lams[np.argmin(vals)], vals.min())
print("smin range:", vals.min(), vals.max())
