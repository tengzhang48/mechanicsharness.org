"""Scan shooting ICs to find all roots with A(pi)=C(pi)=0 and classify by I_n>0."""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import fsolve

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

def endres(u):
    s=solve_ivp(rhs,[0.0,np.pi],[1.0,0.0,u[0],u[1]],rtol=1e-10,atol=1e-12,
                method='Radau',max_step=0.02)
    if not s.success: return np.array([1e6,1e6])
    phi,d1,d2,d3=s.y[:,-1]
    return np.array([k*phi+d2, -(k-1.0)*d1])

def classify(u):
    s=solve_ivp(rhs,[0.0,np.pi],[1.0,0.0,u[0],u[1]],rtol=1e-12,atol=1e-14,
                method='Radau',dense_output=True)
    th=s.t; phi,d1,d2,d3=s.y
    A=k*phi+d2; B=k*(k-1.0)*phi; C=-(k-1.0)*d1
    se=np.sqrt(0.75*(A-B)**2+3.0*C**2); sc=1.0/se.max()
    A,B,C,se=A*sc,B*sc,C*sc,se*sc
    q=se**(n-1.0); Wt=(n/(n+1.0))*se**(n+1)
    Ur=0.75*q*(A-B)/p; dUr=np.gradient(Ur,th); Ut=(dUr-3.0*q*C)/(p+1.0)
    ur_x=np.cos(th)*p*Ur-np.sin(th)*(dUr-Ut)
    I_n=2.0*np.trapezoid(Wt*np.cos(th)-A*ur_x,th)
    return I_n, B[0], Ut[-1]

roots=[]
for a in np.linspace(-2.0,2.0,7):
    for b in np.linspace(-1.0,1.0,7):
        try:
            sol,info,ier,msg=fsolve(endres,[a,b],xtol=1e-11,full_output=True)
            r=endres(sol)
            if np.max(np.abs(r))<1e-6:
                key=(round(sol[0],5),round(sol[1],5))
                if not any(abs(sol[0]-x[0])<1e-4 and abs(sol[1]-x[1])<1e-4 for x in roots):
                    I0,B0,Utpi=classify(sol)
                    roots.append((sol[0],sol[1],I0,B0,Utpi))
        except Exception: pass
print("distinct roots (phi2, phi3, I_n, B0, Ut(pi)):")
for r in roots:
    print("  phi2=%+.6f phi3=%+.6f  I_n=%+.5f  B(0)=%.5f  Ut(pi)=%+.5f" % r)