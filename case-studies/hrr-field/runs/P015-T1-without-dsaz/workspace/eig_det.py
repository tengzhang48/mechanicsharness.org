"""Eigenvalue condition: for even solutions, parametrized by (phi(0),phi''(0)),
the free-surface conditions A(pi)=0, C(pi)=0 are a 2x2 linear homogeneous
system. Its determinant must vanish -> eigenvalue p.
Basis: sol1 with (phi0,phi20)=(1,0); sol2 with (0,1). (phi'(0)=phi'''(0)=0.)
Compute D(p) = A1*C2 - A2*C1 and find its zeros."""
import numpy as np, sys
import hrr_clean as H
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

n = H.N
def even_sols(p):
    k=2.0-p
    def shoot(phi0,phi20):
        s=solve_ivp(lambda th,Y: H.rhs(th,Y,n=n,p=p,k=k),[0.0,np.pi],
                    [phi0,0.0,phi20,0.0], rtol=1e-10, atol=1e-13,
                    method='Radau', max_step=0.01)
        phi,d1,d2,d3=s.y[:,-1]
        return k*phi+d2, -(k-1.0)*d1
    A1,C1=shoot(1.0,0.0)
    A2,C2=shoot(0.0,1.0)
    return A1*C2-A2*C1, (A1,C1,A2,C2)

print(" p        D(p)         A1      C1      A2      C2")
ps=[0.5,0.4,0.3,0.25,0.2,0.15,0.1333,0.12,0.1,0.08,0.05,0.02]
prev=None
for p in ps:
    D,vals=even_sols(p)
    print(f" {p:.5f}  {D:+.8f}  {vals[0]:+.5f} {vals[1]:+.5f} {vals[2]:+.5f} {vals[3]:+.5f}")
    sys.stdout.flush()