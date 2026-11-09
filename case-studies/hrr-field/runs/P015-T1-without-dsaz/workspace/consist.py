"""Verify constitutive consistency of the two derivative choices."""
import numpy as np
from scipy.integrate import solve_ivp
import hrr_clean as H
n=6.5; p=1.0/(n+1.0); k=2.0-p
s=solve_ivp(lambda t,y: np.array(H.rhs(t,y)),[0.0,np.pi],[1.0,0.0,-0.9142244740,0.0328463401],
            rtol=1e-13,atol=1e-15,method='Radau',dense_output=True,max_step=0.002)
sol=s.sol
N=1000001; tt=np.linspace(1e-9,np.pi,N); th=tt[1]-tt[0]
Yv=sol(tt); ph,d1,d2,d3=Yv
A=k*ph+d2; B=k*(k-1.0)*ph; C=-(k-1.0)*d1
se=np.sqrt(0.75*(A-B)**2+3.0*C**2); q=se**(n-1.0)
Ap=k*d1+d3; Bp=k*(k-1.0)*d1; Cp=-(k-1.0)*d2
sep=(0.75*(A-B)*(Ap-Bp)+3.0*C*Cp)/se; qp=(n-1.0)*se**(n-2.0)*sep
Ur=0.75*q*(A-B)/p
dUr=(0.75/p)*(qp*(A-B)+q*(Ap-Bp))
Ut=(3.0*q*C-dUr)/(p-1.0)
dUt=0.75*q*(B-A)-Ur
# check constitutive residuals:
print("eq(I)  pUr-0.75q(A-B):",np.max(np.abs(p*Ur-0.75*q*(A-B))))
print("eq(II) Ur+dUt-0.75q(B-A):",np.max(np.abs(Ur+dUt-0.75*q*(B-A))))
print("eq(III) dUr+(p-1)Ut-3qC:",np.max(np.abs(dUr+(p-1.0)*Ut-3.0*q*C)))
# compare convergence value
def J(mode):
    if mode=='fd':
        dUr2=np.gradient(Ur,th); Ut2=(3.0*q*C-dUr2)/(p-1.0); dUt2=np.gradient(Ut2,th)
    else:
        dUr2,Ut2,dUt2=dUr,Ut,dUt
    ur_x=np.cos(tt)*p*Ur-np.sin(tt)*(dUr2-Ut2)
    ut_x=np.cos(tt)*p*Ut2-np.sin(tt)*(dUt2+Ur)
    integ=(n/(n+1.0))*se**(n+1)*np.cos(tt)-(A*ur_x+C*ut_x)
    return 2*np.trapezoid(integ,tt)
print("I_n analytic=%.6f  fd=%.6f"%(J('a'),J('fd')))