"""Scan p: for each p, solve A(pi;t)=0 for t (shoot with phi(0)=1,phi'(0)=0,
phi'''(0)=0), then evaluate C(pi). Find p where C(pi)=0."""
import numpy as np
import hrr_clean as H
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

n = H.N

def shoot(t, p):
    k = 2.0-p
    Y0=[1.0,0.0,t,0.0]
    s=solve_ivp(lambda th,Y: H.rhs(th,Y,n=n,p=p,k=k),
                [0.0,np.pi], Y0, rtol=1e-11, atol=1e-13,
                method='Radau', max_step=0.005)
    return s

def A_and_C(t, p):
    s=shoot(t,p)
    phi,d1,d2,d3=s.y[:,-1]
    k=2.0-p
    return k*phi+d2, -(k-1.0)*d1, d1

def solve_t_for_A(p):
    # find t with A(pi)=0
    ts=np.linspace(-1.5,1.5,31)
    vals=[A_and_C(t,p)[0] for t in ts]
    roots=[]
    for i in range(len(ts)-1):
        if vals[i]*vals[i+1]<0:
            r=brentq(lambda t: A_and_C(t,p)[0], ts[i], ts[i+1], xtol=1e-12)
            roots.append(r)
    return roots

print(" p      t*        A(pi)      C(pi)=-(k-1)phi'(pi)")
for p in [0.6,0.5,0.4,0.3,0.25,0.2,0.15,0.13,0.1333,0.12,0.1,0.05]:
    roots=solve_t_for_A(p)
    if not roots:
        print(f" {p:.4f}   (no A-root in [-1.5,1.5])"); continue
    for r in roots:
        A,C,d1=A_and_C(r,p)
        print(f" {p:.4f}  {r:+.5f}   {A:+.2e}   {C:+.6f}")