"""
CORRECTED full computation (fixed epsilon_rtheta sign) for n=6.5.
ODE verified at n=1 against Williams biharmonic (residual 4e-16).
Continuation in n from 1 -> 6.5.
"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import fsolve

N=6.5

def make(nn):
    p=1.0/(nn+1.0); k=2.0-p
    def rhs(theta,Y):
        phi,d1,d2,d3=Y
        A=k*phi+d2; B=k*(k-1.0)*phi; C=-(k-1.0)*d1
        se=np.sqrt(0.75*(A-B)**2+3.0*C**2)
        Ap=k*d1+d3; Bp=k*(k-1.0)*d1; Cp=-(k-1.0)*d2
        sep=(0.75*2.0*(A-B)*(Ap-Bp)+3.0*2.0*C*Cp)/(2.0*se)
        q=se**(nn-1.0); qp=(nn-1.0)*se**(nn-2.0)*sep
        Ur=(3.0/4.0)*q*(A-B)/p
        # eq3: Ur'+(p-1)Ut=3qC ; eq2: Ut'=(3/4)q(B-A)-Ur
        # -> Ur'' = 3(q'C+qC') - (p-1)[(3/4)q(B-A)-Ur]
        Urpp = 3.0*(qp*C+q*Cp) - (p-1.0)*((3.0/4.0)*q*(B-A)-Ur)
        Bpp=k*(k-1.0)*d2; Cpp=-(k-1.0)*d3
        UmB=A-B; UmB_p=Ap-Bp; UmB_pp_noD3=k*d2-Bpp
        se2pp=0.75*(2.0*UmB_p*UmB_p+2.0*UmB*UmB_pp_noD3)+3.0*(2.0*Cp*Cp+2.0*C*Cpp)
        se2pp_d3=1.5*UmB
        sepp_noD3=(se2pp-2.0*sep*sep)/(2.0*se); sepp_d3=se2pp_d3/(2.0*se)
        qpp_noD3=(nn-1.0)*((nn-2.0)*se**(nn-3.0)*sep*sep+se**(nn-2.0)*sepp_noD3)
        qpp_d3=(nn-1.0)*se**(nn-2.0)*sepp_d3
        lhs_noD3=qpp_noD3*UmB+2.0*qp*UmB_p+q*(k*d2-Bpp)
        lhs_d3=qpp_d3*UmB+q
        # (3/4)/p [lhs_noD3 + lhs_d3*d3p] = Urpp
        d3p=((4.0/3.0)*p*Urpp-lhs_noD3)/lhs_d3
        return [d1,d2,d3,d3p]
    def endres(u, rtol=1e-10):
        s=solve_ivp(rhs,[0.0,np.pi],[1.0,0.0,u[0],u[1]],rtol=rtol,atol=1e-12,
                    method='Radau',max_step=0.03)
        ph,d1,d2,d3=s.y[:,-1]
        return np.array([k*ph+d2, -(k-1.0)*d1])
    return p,k,rhs,endres

# start at n=1 with Williams ICs: phi''(0)=-0.75, phi'''(0)=0
ns=np.linspace(1.0,6.5,23)
u=np.array([-0.75,0.0])
for nn in ns:
    p,k,rhs,endres=make(nn)
    u,info,ier,msg=fsolve(lambda x:endres(x),u,xtol=1e-10,full_output=True)
    r=endres(u,1e-12)
    print("n=%5.2f  phi''(0)=%.8f phi'''(0)=%.8f resid=%.2e"%(nn,u[0],u[1],np.max(np.abs(r))))

# final at n=6.5
nn=6.5; p,k,rhs,endres=make(nn)
u,info,ier,msg=fsolve(lambda x:endres(x),u,xtol=1e-12,full_output=True)
print("\nFINAL n=6.5: phi''(0)=%.10f phi'''(0)=%.10f resid=%.2e"%(u[0],u[1],np.max(np.abs(endres(u,1e-12)))))
np.savez('sol6p5_correct.npz',u=u)