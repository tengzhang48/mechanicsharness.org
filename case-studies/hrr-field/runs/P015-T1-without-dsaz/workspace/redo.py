"""
CORRECTED derivation (fixing epsilon_rtheta sign error).
Kinematics (TASK.md): 
  eps_rr = u_{r,r}                                  = L r^{p-1} p Ur
  eps_tt = (u_r + u_{t,t})/r                        = L r^{p-1} (Ur + Ut')
  2 eps_rt = u_{r,t}/r + u_{t,r} - u_t/r            = L r^{p-1} (Ur' + p Ut - Ut)
                                                      = L r^{p-1}(Ur' + (p-1)Ut)
Constitutive:  eps_ij = (3/2) alpha eps0 se^{n-1} s_ij   (L = alpha eps0 K^n)
  with s_rr=(A-B)/2, s_tt=(B-A)/2, s_rt=C
=> p Ur = (3/4) q (A-B)
   Ur + Ut' = (3/4) q (B-A)
   Ur' + (p-1) Ut = 3 q C
   q = se^{n-1}
Eliminate Ur,Ut -> 4th order ODE for phi.
"""
import numpy as np
n=6.5; p=1.0/(n+1.0); k=2.0-p

def rhs(theta,Y):
    phi,d1,d2,d3=Y
    A=k*phi+d2; B=k*(k-1.0)*phi; C=-(k-1.0)*d1
    se=np.sqrt(0.75*(A-B)**2+3.0*C**2)
    Ap=k*d1+d3; Bp=k*(k-1.0)*d1; Cp=-(k-1.0)*d2
    sep=(0.75*2.0*(A-B)*(Ap-Bp)+3.0*2.0*C*Cp)/(2.0*se)
    q=se**(n-1.0); qp=(n-1.0)*se**(n-2.0)*sep
    # Ur = (3/4) q (A-B)/p
    Ur=(3.0/4.0)*q*(A-B)/p
    # eq3: Ur' + (p-1)Ut = 3qC  -> Ut = (3qC - Ur')/(p-1)
    # Ut' = (3(qp C + q Cp) - Ur'')/(p-1)
    # eq2: Ut' = (3/4)q(B-A) - Ur
    # => (3(qpC+qCp) - Ur'')/(p-1) = (3/4)q(B-A) - Ur
    # => Ur'' = 3(qpC+qCp) - (p-1)[(3/4)q(B-A) - Ur]
    Urpp = 3.0*(qp*C+q*Cp) - (p-1.0)*((3.0/4.0)*q*(B-A) - Ur)
    # Ur'' = (3/4)/p [ q'' (A-B) + 2 q' (A'-B') + q (A''-B'') ]
    Bpp=k*(k-1.0)*d2; Cpp=-(k-1.0)*d3
    UmB=A-B; UmB_p=Ap-Bp; UmB_pp_noD3=k*d2-Bpp
    se2pp=0.75*(2.0*UmB_p*UmB_p+2.0*UmB*UmB_pp_noD3)+3.0*(2.0*Cp*Cp+2.0*C*Cpp)
    se2pp_d3=1.5*UmB
    sepp_noD3=(se2pp-2.0*sep*sep)/(2.0*se); sepp_d3=se2pp_d3/(2.0*se)
    qpp_noD3=(n-1.0)*((n-2.0)*se**(n-3.0)*sep*sep+se**(n-2.0)*sepp_noD3)
    qpp_d3=(n-1.0)*se**(n-2.0)*sepp_d3
    lhs_noD3=qpp_noD3*UmB+2.0*qp*UmB_p+q*(k*d2-Bpp)
    lhs_d3=qpp_d3*UmB+q
    # (3/4)/p [lhs_noD3 + lhs_d3*d3p] = Urpp
    d3p=((4.0/3.0)*p*Urpp - lhs_noD3)/lhs_d3
    return [d1,d2,d3,d3p]

if __name__=="__main__":
    from scipy.integrate import solve_ivp
    from scipy.optimize import fsolve
    def resid(u):
        s=solve_ivp(rhs,[0.0,np.pi],[1.0,0.0,u[0],u[1]],rtol=1e-8,atol=1e-10,method='Radau')
        phi,d1,d2,d3=s.y[:,-1]
        return [k*phi+d2, -(k-1.0)*d1]
    sol,info,ier,msg=fsolve(resid,[-0.9,0.03],xtol=1e-12,full_output=True)
    print("root: phi''(0)=%.9f phi'''(0)=%.9f resid=%s"%(sol[0],sol[1],resid(sol)))