import numpy as np
import scipy.linalg as sla
from gbcore import cheb_lobatto_01, diff_mat, bary_weights, lagrange_eval, Rfun
from eig import build_Q_general

def d3_mode(N=96):
    xg = cheb_lobatto_01(N); D = diff_mat(xg); D2=D@D; D3=D2@D
    Q = build_Q_general(xg[1:N-1], xg)
    Lm = np.zeros((N+1,N+1)); Bm = np.zeros((N+1,N+1))
    Lm[1:N-1,:] = -(Q@D3); Bm[1:N-1,1:N-1]=np.eye(N-2)
    Lm[0,0]=1.0; Lm[N-1,:]=D2[N,:]; Lm[N,:]=D[N,:]
    lam, V = sla.eig(Lm,Bm)
    i0 = np.argmin(np.abs(lam-11.7427)); f = np.real(V[:,i0]); f/=f[-1]
    w = bary_weights(xg)
    def fev(x):
        return (lagrange_eval(xg,w,np.array([x]))@f)[0]
    return lam[i0].real, fev(0.5)/fev(1.0), [fev(x) for x in [0.1,0.25,0.5,0.75,0.9]]

lam1, r_half, prof = d3_mode(96)
print(f"D3 pencil N=96: lambda1={lam1:.4f} (extrapolates to 11.7427), f(0.5)/f(1)={r_half:.5f}")
print("profile:", [f'{p:.5f}' for p in prof])
print("g-form values:      ['0.21316','0.48341','0.80196','0.96326','0.99619']")
