"""Compare angular-form integrand vs direct Cartesian finite-difference
integrand, pointwise, using identical source fields (no parity assumptions)."""
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
tt=np.linspace(0,np.pi,20001)
Yv=sol(tt); ph,d1,d2,d3=Yv
sev=np.sqrt(0.75*(k*ph+d2-(k*(k-1.0)*ph))**2+3.0*(-(k-1.0)*d1)**2)
sc=1.0/sev.max()
Sig0=1.0; C0=1.0  # work in normalized units; J in these units -> ratio to I_n

def pol(thd):
    Y=sol(abs(thd)); ph,d1,d2,d3=Y
    A=(k*ph+d2)*sc; B=(k*(k-1.0)*ph)*sc; C=(-(k-1.0)*d1)*sc
    se=np.sqrt(0.75*(A-B)**2+3.0*C**2); q=se**(n-1.0)
    Ur=0.75*q*(A-B)/p
    # dUr via central diff on 0..pi using sol
    h=1e-6
    Yh=sol(abs(thd)+h); phb=h and (k*Yh[0]+Yh[2])*sc; Bb=(k*(k-1.0)*Yh[0])*sc
    seb=np.sqrt(0.75*(phb-Bb)**2+3.0*((-(k-1.0)*Yh[1])*sc)**2); qb=seb**(n-1.0)
    Urb=0.75*qb*(phb-Bb)/p
    dUr=(Urb-Ur)/h if abs(thd)+h<=np.pi else 0.0
    Ut=(dUr-3.0*q*C)/(p+1.0)
    return A,B,C,se,q,Ur,Ut,dUr

def t_u_x(thd, R):
    """Return (W n_x, A*ur_x, C*uth_x) angular integrand values at physical R."""
    A,B,C,se,q,Ur,Ut,dUr=pol(thd)
    W=(n/(n+1.0))*se**(n+1)
    # finite-diff in r for dUr/dr, dUt/dr would need r-dep; but Ur,Ut are th-only.
    # du_r/dx = cos*r^{p-1}*(p Ur) ... use analytic with the -u_theta terms.
    c=np.cos(thd); sn=np.sin(thd)
    ur_x=c*p*Ur-sn*(dUr-Ut)      # radial comp of du/dx
    uth_x=c*p*Ut-sn*(dUr and 0)  # placeholder
    return W*c, A*ur_x

thg=np.linspace(0,np.pi,20001)
gvals=[]
for t in thg:
    A,B,C,se,q,Ur,Ut,dUr=pol(t)
    c=np.cos(t); s_=np.sin(t)
    W=(n/(n+1.0))*se**(n+1)
    ur_x=c*p*Ur-s_*(dUr-Ut)
    gvals.append(W*c-A*ur_x)
gvals=np.array(gvals)
I_n=2.0*np.trapezoid(gvals,thg)
print("ANGULAR I_n (with correct -u_theta term, symmetry) =", I_n)
print("  (contour method gave 5.179081)")