"""Verify disp2 satisfies (I),(II),(III) and compute I_n in schemes.
Scheme A (standard HRR, the one the user quoted):
  I_n = int_{-pi}^{pi} [ (n/(n+1)) se^{n+1} cos(th)
                         - (A ux + C uy) ] dth
where ux,uy are the angular parts of du_x/dx, du_y/dx at the self-similar
level, using displacements consistent with the SAME constitutive normalization.
"""
import numpy as np
import hrr_clean as H
from scipy.integrate import solve_ivp

n=H.N; p=H.P; k=H.K
u=np.array([-0.9142244740,0.0328463401])
s=solve_ivp(lambda th,Y: H.rhs(th,Y,n=n,p=p,k=k),[0.0,np.pi],
            [1.0,0.0,u[0],u[1]], rtol=1e-13, atol=1e-15,
            method='Radau', dense_output=True, max_step=0.002)
def st(th): return s.sol(th)
def qty(th):
    phi,d1,d2,d3=st(th)
    A=k*phi+d2; B=k*(k-1.0)*phi; C=-(k-1.0)*d1
    se=np.sqrt(0.75*(A-B)**2+3*C**2)
    return A,B,C,se,d1,d2,d3

# normalize se(0)=1
A0,B0,C0,se0,_,_,_=qty(0.0); Sc=1.0/se0

def UrUt(th):
    A,B,C,se,d1,d2,d3=qty(th)
    q=se**(n-1.0)
    Ap=k*d1+d3; Bp=k*(k-1.0)*d1; Cp=-(k-1.0)*d2
    sep=(0.75*(A-B)*(Ap-Bp)+3*C*Cp)/se
    qp=(n-1.0)*se**(n-2.0)*sep
    Ur=q*(A-B)/(2*p)
    Urp=(qp*(A-B)+q*(Ap-Bp))/(2*p)
    Ut=(2.0*q*C-Urp)/(p-1.0)     # from (III)
    return Ur,Ut,Urp

# --- verify (II): Ur + Ut' = q(B-A)/2 ---
print("verify (II) and (III):")
h=1e-6
for th in [0.5,1.2,2.0,2.8]:
    A,B,C,se,d1,d2,d3=qty(th)
    q=se**(n-1.0)
    Ur,Ut,Urp=UrUt(th)
    Ur1,Ut1,_=UrUt(th+h); Ur2,Ut2,_=UrUt(th-h)
    Utp=(Ut1-Ut2)/(2*h)
    lhs2=Ur+Utp; rhs2=q*(B-A)/2
    lhs3=Urp+(p-1.0)*Ut; rhs3=2.0*q*C
    print(f"  th={th}: (II) {lhs2:.4f} vs {rhs2:.4f} | (III) {lhs3:.4f} vs {rhs3:.4f}")

# --- I_n ---
# J = int [ W cos - t_i du_i/dx ] r dth.  Angular parts (r powers cancel):
# sigma = sigma0 K r^{-p}(A,B,C); u = C0 r^{p}(Ur,Ut) with C0 = alpha eps0 K^n.
# W = n/(n+1) sigma:eps.  sigma:eps = sigma0 K r^{-p} * C0 r^{p-1} [A e_rr+...]
#   e_rr = p Ur ; e_tt=Ur+Ut'; e_rt=(Ur'+(p-1)Ut)/2 (=qC)
#   sigma:eps angular = (sigma0 K)(C0) r^{-1} [ A pUr + B(Ur+Ut') + 2C e_rt ]
#   but A pUr = q(A-B)A/2 ... ; the full contraction simplifies. Use se relation:
#   sigma:eps = sigma_ij eps_ij = (C0/sigma0 r^{-1}) sigma_ij s_ij se^{n-1} *? 
# Simpler: sigma:eps = se^{n+1} (dimensionless, since eps=(3/2)alpha eps0(...)s/sigma0
#   and sigma:s = se^2 ... => sigma:eps=(1/sigma0)*C_n*se^{n-1}*sigma_ij s_ij...).
# We adopt the standard normalized definition where Wt= n/(n+1) se^{n+1}.
def integrand_standard(th):
    A,B,C,se,d1,d2,d3=qty(th)
    a=Sc*A; b=Sc*B; c=Sc*C; sen=Sc*se
    Ur,Ut,Urp=UrUt(th)
    urn=Sc*Ur; utn=Sc*Ut; urnp=Sc*Urp
    # du_r/dx, du_t/dx at self-similar level:
    # u_r=r^p Ur -> du_r/dx = r^{p-1}(p Ur cos - Ur' sin)
    # angular: p urn cos - urnp sin ; similarly for u_t
    dUr_dx = p*urn*np.cos(th) - urnp*np.sin(th)
    Ur1,Ut1,_=UrUt(th+h); Ur2,Ut2,_=UrUt(th-h)
    utnp=Sc*(Ut1-Ut2)/(2*h)
    dUt_dx = p*utn*np.cos(th) - utnp*np.sin(th)
    dux_dx = dUr_dx*np.cos(th) - dUt_dx*np.sin(th)
    duy_dx = dUr_dx*np.sin(th) + dUt_dx*np.cos(th)
    tr=a*np.cos(th)+c*np.sin(th)     # t_r = A cos + C sin
    tt=c*np.cos(th)+b*np.sin(th)     # t_t = C cos + B sin
    Wt=(n/(n+1.0))*sen**(n+1.0)*np.cos(th)
    return Wt - (tr*dux_dx + tt*duy_dx)

from scipy.integrate import quad
vals=[]
for th in np.linspace(0.001,np.pi-0.001,400):
    vals.append(integrand_standard(th))
vals=np.array(vals)
ths=np.linspace(0.001,np.pi-0.001,400)
I_half=np.trapezoid(vals,ths)
print("I_n (int over [0,pi] doubled) =",2*I_half)
print("I_n (int over [0,pi]) =",I_half)