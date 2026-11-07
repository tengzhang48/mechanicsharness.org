import numpy as np
from eig import eigen_s
from gbcore import lagrange_eval, bary_weights
import scipy.linalg as sla

N = 72
xg, lam, V, lams = eigen_s(N=N, verbose=False)
Lmat = None
# rebuild to get eigenvector: recompute pencil
from gbcore import cheb_lobatto_01, diff_mat
from eig import build_Q_general
xg = cheb_lobatto_01(N); D = diff_mat(xg); D2 = D@D; D3 = D2@D
Q = build_Q_general(xg[1:N-1], xg)
G = -(Q @ D3)
Lm = np.zeros((N+1,N+1)); Bm = np.zeros((N+1,N+1))
Lm[1:N-1,:] = G; Bm[1:N-1,1:N-1] = np.eye(N-2)
Lm[0,0]=1.0; Lm[N-1,:]=D2[N,:]; Lm[N,:]=D[N,:]
lam, V = sla.eig(Lm, Bm)
fin = lam[np.isfinite(lam)]
lams_all = lam
# find eigenvalue closest to 11.75
i0 = np.argmin(np.abs(lam - 11.75))
f = V[:, i0].real
f = f / f[-1]
print("lambda1:", lam[i0])
# endpoint slopes: |f| vs x near 0
for j in [1,2,3,4,5,8,12]:
    print(f"x={xg[j]:.3e}  f={f[j]:+.4e}")
lx = np.log(xg[1:8]); lf = np.log(np.abs(f[1:8]))
s0 = np.polyfit(lx, lf, 1)[0]
print("fitted slope near x=0:", s0)
lx2 = np.log(1-xg[-2:-8:-1]); lf2 = np.log(np.abs(f[-2:-8:-1]))
s1 = np.polyfit(lx2, lf2, 1)[0]
print("fitted slope near x=1 (of f itself):", s1, " f(1)=",f[-1])
# also f''' behavior
f3 = D3 @ f
for j in [-1,-2,-3,-5,-8]:
    print(f"1-x={1-xg[j]:.3e} f'''={f3[j]:+.3e}")
