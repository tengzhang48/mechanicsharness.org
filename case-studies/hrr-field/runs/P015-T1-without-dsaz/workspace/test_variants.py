"""Test candidate ICs against candidate ODE variants at n=6.5."""
import numpy as np
from scipy.integrate import solve_ivp

def make(nn, eq3sign):
    # eq3sign: 'A' -> Ur'+(p-1)Ut=3qC   (TASK kinematics, teammate)
    #          'B' -> Ur'-(p+1)Ut=3qC   (my old wrong)
    p=1.0/(nn+1.0); k=2.0-p
    def rhs(theta,Y):
        phi,d1,d2,d3=Y
        A=k*phi+d2; B=k*(k-1.0)*phi; C=-(k-1.0)*d1
        se=np.sqrt(0.75*(A-B)**2+3.0*C**2)
        Ap=k*d1+d3; Bp=k*(k-1.0)*d1; Cp=-(k-1.0)*d2
        sep=(0.75*2.0*(A-B)*(Ap-Bp)+3.0*2.0*C*Cp)/(2.0*se)
        q=se**(nn-1.0); qp=(nn-1.0)*se**(nn-2.0)*sep
        Ur=(3.0/4.0)*q*(A-B)/p
        if eq3sign=='A':
            # Ur' = 3qC - (p-1)Ut = 3qC+(1-p)Ut ; Ut=(Ur'-3qC)/(1-p)... use
            # eq2: Ut'=(3/4)q(B-A)-Ur ; eq3 diff: Ur''=3(q'C+qC')+(1-p)Ut'
            Utp=(3.0/4.0)*q*(B-A)-Ur
            Urpp=3.0*(qp*C+q*Cp)+(1.0-p)*Utp
        else:
            Utp=(3.0/4.0)*q*(B-A)-Ur
            Urpp=3.0*(qp*C+q*Cp)-(p+1.0)*Utp  # d/dth[Ur'-(p+1)Ut=3qC]
        Bpp=k*(k-1.0)*d2; Cpp=-(k-1.0)*d3
        UmB=A-B; UmB_p=Ap-Bp; UmB_pp_noD3=k*d2-Bpp
        se2pp=0.75*(2.0*UmB_p*UmB_p+2.0*UmB*UmB_pp_noD3)+3.0*(2.0*Cp*Cp+2.0*C*Cpp)
        se2pp_d3=1.5*UmB
        sepp_noD3=(se2pp-2.0*sep*sep)/(2.0*se); sepp_d3=se2pp_d3/(2.0*se)
        qpp_noD3=(nn-1.0)*((nn-2.0)*se**(nn-3.0)*sep*sep+se**(nn-2.0)*sepp_noD3)
        qpp_d3=(nn-1.0)*se**(nn-2.0)*sepp_d3
        lhs_noD3=qpp_noD3*UmB+2.0*qp*UmB_p+q*(k*d2-Bpp)
        lhs_d3=qpp_d3*UmB+q
        d3p=((4.0/3.0)*p*Urpp-lhs_noD3)/lhs_d3
        return [d1,d2,d3,d3p]
    def endres(u,rtol=1e-9):
        s=solve_ivp(rhs,[0.0,np.pi],[1.0,0.0,u[0],u[1]],rtol=rtol,atol=1e-11,
                    method='Radau',max_step=0.02)
        if not s.success: return np.array([1e6,1e6])
        ph,d1,d2,d3=s.y[:,-1]
        return np.array([k*ph+d2,-(k-1.0)*d1])
    return p,k,rhs,endres

ics=[[-0.9142244740,0.0328463401],[-0.7050162232,0.0932567960],[-0.4684958,-108.4]]
for name in ['A','B']:
    p,k,rhs,endres=make(6.5,name)
    print("=== ODE variant %s ==="%name)
    for ic in ics:
        r=endres(ic)
        print("  ic=%s -> A(pi)=%.4e C(pi)=%.4e"%(ic,r[0],r[1]))