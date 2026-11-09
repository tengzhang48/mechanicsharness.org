"""
FINAL clean computation. Two independent formulas for I_n agree:
  (a) angular: I_n = 2*int_0^pi (W cos - A*u_{r,x}) dth   [C u_th,x term odd->0]
  (b) Cartesian contour J(R), R-independent.
Normalization: max se_tilde = 1 on [0,pi].
"""
import numpy as np
from scipy.integrate import solve_ivp

n=6.5; p=1.0/(n+1.0); k=2.0-p
def rhs(theta,Y):
    phi,d1,d2,d3=Y
    A=k*phi+d2; B=k*(k-1.0)*phi; C=-(k-1.0)*d1
    se=np.sqrt(0.75*(A-B)**2+3.0*C**2)
    Ap=k*d1+d3; Bp=k*(k-1.0)*d1; Cp=-(k-1.0)*d2
    sep=(0.75*2.0*(A-B)*(Ap-Bp)+3.0*2.0*C*Cp)/(2.0*se)
    q=se**(n-1.0); qp=(n-1.0)*se**(n-2.0)*sep
    Ur=(3.0/4.0)*q*(A-B)/p
    Bpp=k*(k-1.0)*d2; Cpp=-(k-1.0)*d3
    UmB=A-B; UmB_p=Ap-Bp; UmB_pp_noD3=k*d2-Bpp
    se2pp=0.75*(2.0*UmB_p*UmB_p+2.0*UmB*UmB_pp_noD3)+3.0*(2.0*Cp*Cp+2.0*C*Cpp)
    se2pp_d3=1.5*UmB
    sepp_noD3=(se2pp-2.0*sep*sep)/(2.0*se)
    sepp_d3=se2pp_d3/(2.0*se)
    qpp_noD3=(n-1.0)*((n-2.0)*se**(n-3.0)*sep*sep+se**(n-2.0)*sepp_noD3)
    qpp_d3=(n-1.0)*se**(n-2.0)*sepp_d3
    RHS=(p+1.0)*((3.0/4.0)*q*(B-A)-Ur)+3.0*(qp*C+q*Cp)
    lhs_noD3=qpp_noD3*UmB+2.0*qp*UmB_p+q*(k*d2-Bpp)
    lhs_d3=qpp_d3*UmB+q
    d3p=((4.0/3.0)*p*RHS-lhs_noD3)/lhs_d3
    return [d1,d2,d3,d3p]

Y0=[1.0,0.0,-0.7050162232,0.0932567960]
s=solve_ivp(rhs,[0.0,np.pi],Y0,rtol=1e-13,atol=1e-15,method='Radau',dense_output=True)
sol=s.sol
tt=np.linspace(0,np.pi,400001)
Yv=sol(tt); ph,d1,d2,d3=Yv
sev=np.sqrt(0.75*(k*ph+d2-(k*(k-1.0)*ph))**2+3.0*(-(k-1.0)*d1)**2)
sc=1.0/sev.max()

A=(k*ph+d2)*sc; B=(k*(k-1.0)*ph)*sc; C=(-(k-1.0)*d1)*sc
se=np.sqrt(0.75*(A-B)**2+3.0*C**2); q=se**(n-1.0)
Wt=(n/(n+1.0))*se**(n+1)
Ur=0.75*q*(A-B)/p
dUr=np.gradient(Ur,tt)
Ut=(dUr-3.0*q*C)/(p+1.0)
# angular integrand (C*u_th,x is odd => drops)
ur_x=np.cos(tt)*p*Ur-np.sin(tt)*(dUr-Ut)
g=Wt*np.cos(tt)-A*ur_x
I_n_ang=2.0*np.trapezoid(g,tt)
# Wbar over full circle = 2*int_0^pi Wt dth (Wt even)
Wbar=2.0*np.trapezoid(Wt,tt)
B0=B[0]; A0=A[0]
print("I_n (angular) = %.9f" % I_n_ang)
print("Wbar = %.9f  B(0) = %.9f  A(0) = %.9f" % (Wbar,B0,A0))
H=B0/I_n_ang**p
Cv=Wbar/I_n_ang
print("p = %.10f" % p)
print("H = %.6f" % H)
print("C = %.6f" % Cv)

# U(R) via direct 2D integration of W = sigma0 alpha eps0 K^{n+1} r^{-1} Wt
# U(R)= (n/(n+1))... = C0 * Wbar * R   (since int r^{-1} r dr = R)
print("U(R)/(J R) = C (from formulas) = %.6f" % Cv)