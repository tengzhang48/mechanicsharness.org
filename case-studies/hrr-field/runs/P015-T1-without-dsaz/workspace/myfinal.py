"""FINAL HRR n=6.5 angular functions + I_n. Corrected derivation (u~r^p,
eps~r^{p-1}, 2eps_rt = r^{p-1}(Ur' + (p-1)Ut)).  Corner solution via 2 free
ICs (phi''(0),phi'''(0)) and 2 free-surface conditions."""
import numpy as np
import hrr_clean as H
from scipy.integrate import solve_ivp
from scipy.optimize import fsolve

n=H.N; p=H.P; k=H.K
print(f"n={n} p={p:.10f} k={k:.10f}")

def resid(u):
    s=solve_ivp(lambda th,Y: H.rhs(th,Y,n=n,p=p,k=k),[0.0,np.pi],
                [1.0,0.0,u[0],u[1]], rtol=1e-12, atol=1e-14,
                method='Radau', max_step=0.004)
    phi,d1,d2,d3=s.y[:,-1]
    return np.array([k*phi+d2, -(k-1.0)*d1])

u=np.array([-0.9142244740, 0.0328463401])
# quick 2-point polish
for _ in range(3):
    f0=resid(u)
    J=np.zeros((2,2)); h=1e-6
    for j in range(2):
        up=u.copy(); up[j]+=h
        J[:,j]=(resid(up)-f0)/h
    du=np.linalg.solve(J,-f0)
    u=u+du
    print("polish u=",u," |resid|=",np.linalg.norm(resid(u)))
print("phi''(0)=%.10f  phi'''(0)=%.10f"%(u[0],u[1]))
s=solve_ivp(lambda th,Y: H.rhs(th,Y,n=n,p=p,k=k),[0.0,np.pi],
            [1.0,0.0,u[0],u[1]], rtol=1e-13, atol=1e-15,
            method='Radau', dense_output=True, max_step=0.002)

def st(th): return s.sol(th)
def ABC(th):
    phi,d1,d2,d3=st(th)
    return k*phi+d2, k*(k-1.0)*phi, -(k-1.0)*d1

A0,B0,C0=ABC(0.0)
se0=np.sqrt(0.75*(A0-B0)**2+3*C0**2)
Sc=1.0/se0
print("RAW se(0)=%.10f"%se0)
print("NORM se(0)=1 -> A(0)=%.8f B(0)=%.8f C(0)=%.8f"%(Sc*A0,Sc*B0,Sc*C0))

# displacement angular fns (constitutive constant c absorbed; use c=1 relation)
def disp(th):
    phi,d1,d2,d3=st(th)
    A=k*phi+d2; B=k*(k-1.0)*phi; C=-(k-1.0)*d1
    se=np.sqrt(0.75*(A-B)**2+3*C**2)
    q=se**(n-1.0)
    Ap=k*d1+d3; Bp=k*(k-1.0)*d1; Cp=-(k-1.0)*d2
    sep=(0.75*(A-B)*(Ap-Bp)+3*C*Cp)/se
    qp=(n-1.0)*se**(n-2.0)*sep
    Ur=q*(A-B)/(2*p)
    Urp=(qp*(A-B)+q*(Ap-Bp))/(2*p)
    Ut=(Urp-2*q*C)/(p+1.0)  # from (III): Ur'+(p-1)Ut=2qC -> Ut=(2qC-Ur')/(p-1)
    return Ur,Ut,se
# NOTE: (III) gives Ur'+(p-1)Ut=2qC => (p-1)Ut=2qC-Ur' => Ut=(2qC-Ur')/(p-1).
# Equivalent: multiply num/den by -1: Ut=(Ur'-2qC)/(1-p). Different from above!
# Fix: use Ut=(2qC-Ur')/(p-1).
def disp2(th):
    phi,d1,d2,d3=st(th)
    A=k*phi+d2; B=k*(k-1.0)*phi; C=-(k-1.0)*d1
    se=np.sqrt(0.75*(A-B)**2+3*C**2)
    q=se**(n-1.0)
    Ap=k*d1+d3; Bp=k*(k-1.0)*d1; Cp=-(k-1.0)*d2
    sep=(0.75*(A-B)*(Ap-Bp)+3*C*Cp)/se
    qp=(n-1.0)*se**(n-2.0)*sep
    Ur=q*(A-B)/(2*p)
    Urp=(qp*(A-B)+q*(Ap-Bp))/(2*p)
    Ut=(2.0*q*C-Urp)/(p-1.0)
    return Ur,Ut,se
print("check disp vs disp2 at th=1 (III consistency):")
Ur,Ut,se=disp(1.0); Ur2,Ut2,se2=disp2(1.0)
print("  Ut(disp)=%.6f  Ut(disp2)=%.6f"%(Ut,Ut2))
np.savez("my_hrr_sol.npz", u=u, n=n,p=p,k=k, Sc=Sc)

# report shapes
print("\n th    A       B       C       se      Ur      Ut")
for th in np.linspace(0,np.pi,13):
    A,B,C=ABC(th); se=np.sqrt(0.75*(A-B)**2+3*C**2)
    Ur,Ut,_=disp2(th)
    print(" %.3f %+.5f %+.5f %+.5f %+.5f %+.5f %+.5f"%(th,Sc*A,Sc*B,Sc*C,Sc*se,Sc*Ur,Sc*Ut))