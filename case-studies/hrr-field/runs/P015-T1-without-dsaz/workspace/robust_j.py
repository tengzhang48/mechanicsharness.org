"""
ROBUST contour J-integral. Build the full physical field on the circle by
using Mode I parity explicitly (fields defined for th in [0,pi], extended):
  sigma_rr, sigma_tt, sigma_rt, u_r  are EVEN in th (for th in [0,pi] base)
  Actually: define everything from |th| and apply vector rotation properly.
We compute sigma in Cartesian directly and u in Cartesian directly with parity.
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
tt=np.linspace(0,np.pi,40001)
Yv=sol(tt); ph,d1,d2,d3=Yv
sev=np.sqrt(0.75*(k*ph+d2-(k*(k-1.0)*ph))**2+3.0*(-(k-1.0)*d1)**2)
sc=1.0/sev.max()

def at(th):
    """fields at signed angle th (radians, any value). Uses |th| for polar
    quantities (all even in th for Mode I since they're components in the
    polar basis: sig_rr,sig_tt,sig_rt,u_r are even, u_th is odd)."""
    a=abs(th); Y=sol(a); ph,d1,d2,d3=Y
    A=(k*ph+d2)*sc; B=(k*(k-1.0)*ph)*sc; C=(-(k-1.0)*d1)*sc
    se=np.sqrt(0.75*(A-B)**2+3.0*C**2); q=se**(n-1.0)
    Ur=0.75*q*(A-B)/p
    Ur=Ur; Ut=None
    # u_th from eq3 using dUr/dth (even in th for u_r -> dUr/dth odd; but we
    # only need |th| magnitude); compute derivative wrt a
    h=1e-7
    def Ur_of(x):
        Yb=sol(min(x,np.pi)); pb,db1,db2,db3=Yb
        Ab=(k*pb+db2)*sc; Bb=(k*(k-1.0)*pb)*sc; Cb=(-(k-1.0)*db1)*sc
        seb=np.sqrt(0.75*(Ab-Bb)**2+3.0*Cb**2); qb=seb**(n-1.0)
        return 0.75*qb*(Ab-Bb)/p, qb*Cb
    Urm,_=Ur_of(a-h); Urp,_=Ur_of(a+h)
    dUrda=(Urp-Urm)/(2*h)
    qC=Ur_of(a)[1]
    Ut=(dUrda-3.0*qC)/(p+1.0)
    # sign of u_th: odd in th
    sgn=np.sign(th) if th!=0 else 0.0
    return A,B,C,se,Ur,sgn*Ut

def sig_cart(r,th):
    A,B,C,se,Ur,Ut=at(th)
    R=r**(-p)  # sigma0*K=1, alpha eps0 K^n=1 units
    srr=A*R; stt=B*R; srt=C*R
    c=np.cos(th); sn=np.sin(th)
    sxx=c*c*srr-2*c*sn*srt+sn*sn*stt
    syy=sn*sn*srr+2*c*sn*srt+c*c*stt
    sxy=c*sn*(srr-stt)+(c*c-sn*sn)*srt
    return sxx,syy,sxy

def u_cart(r,th):
    A,B,C,se,Ur,Ut=at(th)
    ur=r**p*Ur; ut=r**p*Ut
    c=np.cos(th); sn=np.sin(th)
    ux=c*ur-sn*ut
    uy=sn*ur+c*ut
    return ux,uy

def W(r,th):
    A,B,C,se,Ur,Ut=at(th)
    return r**(-1.0)*(n/(n+1.0))*se**(n+1)

def J(R,N=2001):
    th=np.linspace(-np.pi,np.pi,N)
    dth=th[1]-th[0]
    eps=1e-6*R
    total=0.0
    for t in th:
        sxx,syy,sxy=sig_cart(R,t); c=np.cos(t); sn=np.sin(t)
        tx=sxx*c+sxy*sn; ty=sxy*c+syy*sn
        x=R*c; y=R*sn
        uxp,uyp=u_cart(np.hypot(x+eps,y),np.arctan2(y,x+eps))
        uxm,uym=u_cart(np.hypot(x-eps,y),np.arctan2(y,x-eps))
        dux=(uxp-uxm)/(2*eps); duy=(uyp-uym)/(2*eps)
        total+=(W(R,t)*c-(tx*dux+ty*duy))*R*dth
    return total

print("Convergence of J(R) in R (sigma0*alpha*eps0*K^{n+1}=1 units):")
for R in [1e-3,1e-2,0.1,1.0]:
    print("  R=%.1e J=%.6f" % (R, J(R,N=3001)))
print("Grid convergence at R=1e-2:")
for N in [501,1001,2001,4001,8001]:
    print("  N=%5d J=%.6f" % (N, J(1e-2,N=N)))