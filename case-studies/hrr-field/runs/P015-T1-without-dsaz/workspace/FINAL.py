"""
DEFINITIVE computation for n=6.5 using the VERIFIED ODE (hrr_clean.rhs).
ODE verified: reproduces n=1 Williams biharmonic exactly AND
phi''(0)=-0.914224474, phi'''(0)=+0.032846340 gives A(pi)=C(pi)=0.
Scalings: sigma~r^{-p}, eps~r^{p-1}, u~r^p, W~r^{-1}. Normalize max se = 1.
"""
import numpy as np
from scipy.integrate import solve_ivp
import hrr_clean as H

n=6.5; p=1.0/(n+1.0); k=2.0-p
Y0=[1.0,0.0,-0.9142244740,0.0328463401]
s=solve_ivp(lambda t,y: np.array(H.rhs(t,y)),[0.0,np.pi],Y0,
            rtol=1e-13,atol=1e-15,method='Radau',dense_output=True,max_step=0.002)
sol=s.sol
tt=np.linspace(0,np.pi,2000001)
Yv=sol(tt); ph,d1,d2,d3=Yv
A=k*ph+d2; B=k*(k-1.0)*ph; C=-(k-1.0)*d1
se=np.sqrt(0.75*(A-B)**2+3.0*C**2)          # RAW (unscaled) tilde-se
q=se**(n-1.0)
Wt=(n/(n+1.0))*se**(n+1)                     # raw tilde-W
print("RAW: se(0)=%.6f se_max=%.6f  A(0)=%.6f B(0)=%.6f"%(se[0],se.max(),A[0],B[0]))
print("free surface: A(pi)=%.2e C(pi)=%.2e"%(A[-1],C[-1]))

# displacement: (I) p Ur = 3/4 q(A-B) ; (III) Ur'+(p-1)Ut = 3qC
Ur=0.75*q*(A-B)/p
Ap=k*d1+d3; Bp=k*(k-1.0)*d1; Cp=-(k-1.0)*d2
sep=(0.75*(A-B)*(Ap-Bp)+3.0*C*Cp)/se; qp=(n-1.0)*se**(n-2.0)*sep
dUr=(0.75/p)*(qp*(A-B)+q*(Ap-Bp))
Ut=(3.0*q*C-dUr)/(p-1.0)
dUt=0.75*q*(B-A)-Ur
print("constitutive residuals: %.2e %.2e %.2e"%(
    np.max(np.abs(p*Ur-0.75*q*(A-B))),
    np.max(np.abs(Ur+dUt-0.75*q*(B-A))),
    np.max(np.abs(dUr+(p-1.0)*Ut-3.0*q*C))))

# ---- angular J integrand (correct polar->cartesian rotation) ----
# u_r=L r^p Ur, u_t=L r^p Ut ; du_i/dx = cos*d/dr - sin/r * d/dth
# u_{r,x}=L r^{p-1}[ cos*p*Ur - sin*(Ur' - Ut) ]
# u_{t,x}=L r^{p-1}[ cos*p*Ut - sin*(Ut' + Ur) ]
ur_x=np.cos(tt)*p*Ur-np.sin(tt)*(dUr-Ut)
ut_x=np.cos(tt)*p*Ut-np.sin(tt)*(dUt+Ur)
# t_i=(sig_rr,sig_rt); t_i u_{i,x} = A*ur_x + C*ut_x
tterm=A*ur_x+C*ut_x
Wcx=Wt*np.cos(tt)
I_n=2.0*np.trapezoid(Wcx-tterm,tt)   # raw I_n (J = alpha*eps0*sig0*K^{n+1} * I_n)
Wbar=2.0*np.trapezoid(Wt,tt)         # raw  (U = alpha*eps0*sig0*K^{n+1} * Wbar * R)
Braw=B[0]
print("I_n = %.8f   Wbar = %.8f   B(0)_raw = %.8f"%(I_n,Wbar,Braw))
# H, C are scale-invariant; compute from raw normalization phi(0)=1:
Hv=Braw/I_n**p
Cv=Wbar/I_n
print("p = %.10f"%p)
print("H = %.6f"%Hv)
print("C = %.6f"%Cv)
# Alternative: Wterm and traction term separately
print("[split] intWcos=%.6f  intTraction=%.6f"%(np.trapezoid(Wcx*0+Wcx,tt)*0 + 2*np.trapezoid(Wcx,tt), 2*np.trapezoid(tterm,tt)))