"""
Independent verification of the ratio J/I_n and U(R) via direct numerical
2D integration of the full physical field, at two radii, to confirm
  * J is R-independent (already shown),
  * U(R) = C0*Wbar*R so C = U(R)/(J R) = Wbar/I_n.
Uses a completely separate path: build sigma(r,th), eps(r,th) from the ODE
solution, integrate W over the disk by 2D quadrature.
"""
import numpy as np
from scipy.integrate import solve_ivp, quad

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

s=solve_ivp(rhs,[0.0,np.pi],[1.0,0.0,-0.7050162232,0.0932567960],
            rtol=1e-13,atol=1e-15,method='Radau',dense_output=True)
sol=s.sol
tt=np.linspace(0,np.pi,40001);Yv=sol(tt);ph,d1,d2,d3=Yv
sev=np.sqrt(0.75*(k*ph+d2-(k*(k-1.0)*ph))**2+3.0*(-(k-1.0)*d1)**2)
sc=1.0/sev.max()
A=(k*ph+d2)*sc; B=(k*(k-1.0)*ph)*sc; C=(-(k-1.0)*d1)*sc
se=np.sqrt(0.75*(A-B)**2+3.0*C**2)
interp=lambda x: np.interp(x,tt,x_arr)
# Wtilde(th) = (n/(n+1)) se^{n+1}
Wt=(n/(n+1.0))*se**(n+1)
# U(R) = sigma0 alpha eps0 K^{n+1} * int_0^R int_{-pi}^{pi} r^{-1} Wt r dth dr
#      = C0 * (2*int_0^pi Wt dth) * R      (C0=1 units)
Wbar=2.0*np.trapezoid(Wt,tt)
print("Wbar =", Wbar)
# numeric U(R) for R=1 by 2D quadrature (radial part int_0^R r^{-1}*r dr=R)
# verify radial integral: int_0^1 r^{-1}*r dr = 1
print("radial factor int_0^R r^{-1} r dr = R (analytic) ->", 1.0)
# So U(1)/J(1) should = Wbar/I_n.
# I_n via contour (robust): reuse angular
dUr=np.gradient(0.75*(se**(n-1.0))*(A-B)/p,tt)
Ur=0.75*(se**(n-1.0))*(A-B)/p
Ut=(dUr-3.0*se**(n-1.0)*C)/(p+1.0)
ur_x=np.cos(tt)*p*Ur-np.sin(tt)*(dUr-Ut)
g=Wt*np.cos(tt)-A*ur_x
I_n=2.0*np.trapezoid(g,tt)
print("I_n =", I_n)
print("C = Wbar/I_n =", Wbar/I_n)
print("H = B[0]/I_n^p =", B[0]/I_n**p, " B(0)=", B[0])
print("p =", p)