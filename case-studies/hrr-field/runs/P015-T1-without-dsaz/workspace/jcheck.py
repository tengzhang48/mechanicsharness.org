"""Independent Cartesian J at several R, using physical fields, check R-independence.
J = int_{-pi}^{pi}[ W cos - t_i u_{i,x} ] r dtheta.
Build physical fields sigma_{ij}(r,th), u_i(r,th) with amplitude L=1, then
t_i u_{i,x} computed by cartesian differentiation with polar-basis rotation.
"""
import numpy as np
from scipy.integrate import solve_ivp
import hrr_clean as H
n=6.5; p=1.0/(n+1.0); k=2.0-p
s=solve_ivp(lambda t,y: np.array(H.rhs(t,y)),[0.0,np.pi],[1.0,0.0,-0.9142244740,0.0328463401],
            rtol=1e-13,atol=1e-15,method='Radau',dense_output=True,max_step=0.002)
sol=s.sol
L=1.0
def fields(R,th):
    Y=sol(th); ph,d1,d2,d3=Y
    A=k*ph+d2; B=k*(k-1)*ph; C=-(k-1)*d1
    se=np.sqrt(0.75*(A-B)**2+3*C**2)
    q=se**(n-1)
    # stress (scale-free tilde times L, times r^-p)
    rp=R**(-p)
    srr=A*L*rp; stt=B*L*rp; srt=C*L*rp
    Ur=0.75*q*(A-B)/p * L
    dUr=np.array([ (0.75*q*(A-B)/p) ])  # placeholder
    # displacement and derivatives wrt r, th via analytic (u~r^p)
    # u_r = L r^p Ur(th): u_{r,r}=L p r^{p-1} Ur; u_{r,th}=L r^p Ur'
    # u_t = L r^p Ut(th)
    return A,B,C,se,q
# Instead integrate in th with fine grid and finite differences for derivatives
tt=np.linspace(0,np.pi,2000001)
Yv=sol(tt); ph,d1,d2,d3=Yv
A=k*ph+d2; B=k*(k-1)*ph; C=-(k-1)*d1
se=np.sqrt(0.75*(A-B)**2+3*C**2); q=se**(n-1)
Ur=0.75*q*(A-B)/p*L
# analytic derivatives (avoid FD noise from q' at axis)
Ap=k*d1+d3; Bp=k*(k-1)*d1; Cp=-(k-1)*d2
dUr=(0.75/p)*(np.gradient(q,tt)*(A-B)+q*(Ap-Bp))
Ut=(3.0*q*C-dUr)/(p-1.0)
dUt=0.75*q*(B-A)-Ur        # from eq (II)
for R in [0.02,0.05,0.1,0.3,0.7]:
    rp=R**(-p)
    srr=A*rp; stt=B*rp; srt=C*rp
    err=p*R**(p-1)*Ur
    ett=R**(p-1)*(Ur+dUt)
    ert=0.5*R**(p-1)*(dUr+(p-1)*Ut)
    W=(n/(n+1.0))*se** (n+1)*R**(-1)
    ur_x=np.cos(tt)*R**(p-1)*(p*Ur)-np.sin(tt)*R**(p-1)*(dUr-Ut)
    ut_x=np.cos(tt)*R**(p-1)*(p*Ut)-np.sin(tt)*R**(p-1)*(dUt+Ur)
    integ=(W*np.cos(tt)-(srr*ur_x+srt*ut_x))*R
    J=2.0*np.trapezoid(integ,tt)
    print("R=%.3f  J=%.8f"%(R,J))
print("I_n analytic full: 1.67653262")