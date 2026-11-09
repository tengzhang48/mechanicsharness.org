"""Fast p-scan: for each p, solve A(pi;t)=0 for t, report C(pi). Flushed."""
import numpy as np, sys
import hrr_clean as H
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

n = H.N
def AC(t,p):
    k=2.0-p
    s=solve_ivp(lambda th,Y: H.rhs(th,Y,n=n,p=p,k=k),[0.0,np.pi],
                [1.0,0.0,t,0.0], rtol=1e-9, atol=1e-11, method='Radau',
                max_step=0.01)
    phi,d1,d2,d3=s.y[:,-1]
    return k*phi+d2, -(k-1.0)*d1, d1

def solve_t_for_A(p, lo=-2, hi=2, N=41):
    ts=np.linspace(lo,hi,N); vals=[AC(t,p)[0] for t in ts]; roots=[]
    for i in range(N-1):
        if vals[i]*vals[i+1]<0:
            r=brentq(lambda t: AC(t,p)[0], ts[i], ts[i+1], xtol=1e-10)
            roots.append(r)
    return roots

print(" p        t*         A(pi)      C(pi)      phi'(pi)")
for p in [0.6,0.5,0.4,0.35,0.3,0.25,0.2,0.1667,0.15,0.14,0.1333,0.12,0.1,0.08,0.05]:
    roots=solve_t_for_A(p)
    if not roots:
        print(f" {p:.5f}   none"); sys.stdout.flush(); continue
    for r in roots:
        A,C,d1=AC(r,p)
        print(f" {p:.5f}  {r:+.6f}  {A:+.2e}  {C:+.6f}  {d1:+.6f}"); sys.stdout.flush()