import numpy as np
import scipy.linalg as sla
from gbcore import cheb_lobatto_01, diff_mat, bary_weights, lagrange_eval, Rfun
from eig import build_Q_general
from numpy.polynomial.legendre import leggauss

N = 96
xg = cheb_lobatto_01(N); D = diff_mat(xg); D2=D@D; D3=D2@D
Q = build_Q_general(xg[1:N-1], xg)
G = -(Q @ D3)
Lm = np.zeros((N+1,N+1)); Bm = np.zeros((N+1,N+1))
Lm[1:N-1,:] = G; Bm[1:N-1,1:N-1]=np.eye(N-2)
Lm[0,0]=1.0; Lm[N-1,:]=D2[N,:]; Lm[N,:]=D[N,:]
lam, V = sla.eig(Lm,Bm)
i0 = np.argmin(np.abs(lam-11.75)); f = np.real(V[:,i0]); f/=f[-1]
f3 = D3@f
# residual on independent dense grid using barycentric interp of f and f3
w = bary_weights(xg)
xr = np.array([0.01,0.05,0.1,0.3,0.5,0.7,0.9,0.95,0.99,0.995])
Lr = lagrange_eval(xg, w, xr)
fr = Lr@f; f3r = Lr@f3
xi, wi = leggauss(4000); xi=0.5*(xi+1); wi=0.5*wi
Lq = lagrange_eval(xg, w, xi)
gq = Lq@f3   # f''' at quadrature nodes
for i,x in enumerate(xr):
    d = Lr[i]  # not needed
    # A[f3](x): PV split
    gx = (lagrange_eval(xg,w,np.array([x]))@f3)[0]
    t1 = np.sum(wi*(gq-gx)/(x-xi))
    t2 = gx*np.log(x/(1-x))
    t3 = np.sum(wi*Rfun(x,xi)*gq)
    lhs = -lam[i0].real*fr[i]
    rhs = t1+t2+t3
    print(f"x={x:6.3f}  -lam*f={lhs:+.5f}  A[f3]={rhs:+.5f}  resid={lhs-rhs:+.2e}")
# endpoint behavior of f near 0: ratio f/x
print("\nf/x near 0:")
for j in [1,2,3,4,6,10,20]:
    print(f"  x={xg[j]:.3e} f/x={f[j]/xg[j]:.5f}")
print("\nf''' near 0 (log-log slope):")
j = np.arange(1,8)
lx = np.log(xg[j]); lf3 = np.log(np.abs(f3[j])+1e-30)
print("  slopes:", np.polyfit(lx,lf3,1)[0], " f3 values:", f3[1:8])
