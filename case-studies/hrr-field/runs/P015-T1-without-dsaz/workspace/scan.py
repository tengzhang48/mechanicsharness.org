"""Shooting scan for n=6.5: integrate from theta=0 with phi(0)=1,phi'(0)=0,
phi'''(0)=0 and varying phi''(0)=t. Record A(pi) and C(pi)."""
import numpy as np
import hrr_clean as H
from scipy.integrate import solve_ivp

n=H.N; p=H.P; k=H.K
print("n=",n," p=",p," k=",k)

def shoot(t, nnt=n):
    Y0=[1.0,0.0,t,0.0]
    s=solve_ivp(lambda th,Y: H.rhs(th,Y,n=nnt,p=1.0/(nnt+1),k=2-1.0/(nnt+1)),
                [0.0, np.pi], Y0, rtol=1e-12, atol=1e-14,
                method='Radau', dense_output=True, max_step=0.005)
    return s

print(" t      A(pi)      C(pi)     phi'(pi)   phi(pi)")
for t in [-2,-1,-0.5,0,0.5,1,2,3,5]:
    s=shoot(t)
    if not s.success:
        print(f" {t:+.2f}  FAIL {s.message}"); continue
    phi,d1,d2,d3=s.y[:,-1]
    A=k*phi+d2; C=-(k-1)*d1
    print(f" {t:+.2f}  {A:+.6f}  {C:+.6f}  {d1:+.6f}  {phi:+.6f}")