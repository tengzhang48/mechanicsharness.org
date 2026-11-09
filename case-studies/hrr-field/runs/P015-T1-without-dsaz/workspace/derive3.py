"""Numerically verify: the compat ODE residual is identically zero for the
Williams family phi = a cos(3t/2)+b cos(t/2) at n=1, using the SAME formula
used in the ODE (Ut from (III), Utp_II from (II), Utp_diff = d/dtheta Ut)."""
import numpy as np

p = 0.5
k = 1.5

def build(t, phi):
    # evaluate A,B,C and derivatives numerically
    def d(m):
        if m==0: return phi(t)
        if m==1: return -1.5*a*np.sin(1.5*t)-0.5*b*np.sin(0.5*t)
        if m==2: return -2.25*a*np.cos(1.5*t)-0.25*b*np.cos(0.5*t)
        if m==3: return 3.375*a*np.sin(1.5*t)+0.125*b*np.sin(0.5*t)
    pp = d(0); d1=d(1); d2=d(2); d3=d(3)
    A=k*pp+d2; B=k*(k-1)*pp; C=-(k-1)*d1
    Ap=k*d1+d3; Bp=k*(k-1)*d1; Cp=-(k-1)*d2
    se=np.sqrt(0.75*(A-B)**2+3*C**2)
    return A,B,C,se,Ap,Bp,Cp

a,b=1.3,0.7
for t in [0.3,1.0,2.0,2.9]:
    A,B,C,se,Ap,Bp,Cp=build(t,lambda x:x)
    q=se**0  # n=1 -> q=1
    UmB=A-B; UmB_p=Ap-Bp
    Ur=q*UmB/(2*p)
    Urp=(0+0)/1  # q'=0 -> Urp = q*UmB_p/(2p)
    Urp=q*UmB_p/(2*p)
    Ut=(2*q*C-Urp)/(p-1)
    Utp_II=q*(B-A)/2-Ur
    # Utp_diff = d/dtheta[(2qC-Urp)/(p-1)]. q=1 const. derivative:
    # = [2(q' C+q C') - Urpp]/(p-1), q'=0 -> [2 C' - Urpp]/(p-1)
    C_p=Cp
    # Urpp = q*UmB_p'/(2p) = q*(Ap-Bp)'/(2p); UmB_pp=(A-B)''=App-Bpp
    # simpler: use finite difference for Urp' from Urp across t
    dt=1e-5
    A2,B2,C2,se2,Ap2,Bp2,Cp2=build(t+dt,lambda x:x)
    Urp2=(Ap2-Bp2)/(2*p)
    Urpp_fd=(Urp2-Urp)/dt
    Utp_diff=(2*C_p-Urpp_fd)/(p-1)
    print(f"t={t}: Utp_II={Utp_II:.6f} Utp_diff={Utp_diff:.6f} diff={Utp_II-Utp_diff:.2e}")