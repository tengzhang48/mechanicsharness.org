"""Convergence study for I_n: analytic vs FD derivatives, refine grid."""
import numpy as np
from scipy.integrate import solve_ivp
import hrr_clean as H
n=6.5; p=1.0/(n+1.0); k=2.0-p
s=solve_ivp(lambda t,y: np.array(H.rhs(t,y)),[0.0,np.pi],[1.0,0.0,-0.9142244740,0.0328463401],
            rtol=1e-13,atol=1e-15,method='Radau',dense_output=True,max_step=0.002)
sol=s.sol
def run(N, mode):
    tt=np.linspace(0,np.pi,N); th=tt[1]-tt[0]
    Yv=sol(tt); ph,d1,d2,d3=Yv
    A=k*ph+d2; B=k*(k-1.0)*ph; C=-(k-1.0)*d1
    se=np.sqrt(0.75*(A-B)**2+3.0*C**2); q=se**(n-1.0)
    Ap=k*d1+d3; Bp=k*(k-1.0)*d1; Cp=-(k-1.0)*d2
    sep=(0.75*(A-B)*(Ap-Bp)+3.0*C*Cp)/se; qp=(n-1.0)*se**(n-2.0)*sep
    Ur=0.75*q*(A-B)/p
    if mode=='analytic':
        dUr=(0.75/p)*(qp*(A-B)+q*(Ap-Bp))
        Ut=(3.0*q*C-dUr)/(p-1.0)
        dUt=0.75*q*(B-A)-Ur
    else:
        dUr=np.gradient(Ur,th)
        Ut=(3.0*q*C-dUr)/(p-1.0)
        dUt=np.gradient(Ut,th)
    ur_x=np.cos(tt)*p*Ur-np.sin(tt)*(dUr-Ut)
    ut_x=np.cos(tt)*p*Ut-np.sin(tt)*(dUt+Ur)
    integ=(n/(n+1.0))*se**(n+1)*np.cos(tt)-(A*ur_x+C*ut_x)
    return np.trapezoid(integ,tt)
for mode in ['analytic','fd']:
    print("mode",mode)
    for N in [20001,100001,500001,2000001]:
        print("  N=%8d I_n=%.8f"%(N,run(N,mode)))