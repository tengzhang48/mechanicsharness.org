"""
Fully independent cross-check of I_n, B(0) by building the PHYSICAL field
sigma_ij(r,th) and u_i(r,th) on a polar grid, then computing the contour
J-integral at radius R numerically from the Cartesian definition:
   J = int_G ( W n_x - t_i du_i/dx ) ds ,
with G a circle of radius R, n = radial unit vector, t = sigma.n.
This does NOT reuse the angular-formula algebra.
Also compute U(R) = int_0^R int W r dth dr to get C = U/(J R).
"""
import numpy as np
from scipy.integrate import solve_ivp, quad

n = 6.5
p = 1.0/(n+1.0)
k = 2.0-p
sigma0=400e6; E=200e9; alpha=0.02; eps0=sigma0/E; K=1.0


def rhs(theta, Y):
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
# normalize so max se = 1
tt=np.linspace(0,np.pi,20001)
Yv=sol(tt); ph,d1,d2,d3=Yv
Av=k*ph+d2; Bv=k*(k-1.0)*ph; Cv=-(k-1.0)*d1
sev=np.sqrt(0.75*(Av-Bv)**2+3.0*Cv**2)
sc=1.0/sev.max()

def angular(thd):
    """return tilde fields at angle thd in [0,pi], scaled."""
    Y=sol(thd); ph,d1,d2,d3=Y
    A=(k*ph+d2)*sc; B=(k*(k-1.0)*ph)*sc; C=(-(k-1.0)*d1)*sc
    se=np.sqrt(0.75*(A-B)**2+3.0*C**2)
    q=se**(n-1.0)
    Ur=0.75*q*(A-B)/p
    # Ut from eq3: (p+1)Ut = Ur' - 3 q C ; Ur' needs derivative -> finite diff
    h=1e-6
    Yp=sol(thd+h); php,d1p,d2p,d3p=Yp
    Ap=(k*php+d2p)*sc; Bp=(k*(k-1.0)*php)*sc; Cp=(-(k-1.0)*d1p)*sc
    sep=np.sqrt(0.75*(Ap-Bp)**2+3.0*Cp**2); qp=sep**(n-1.0)
    Urp=0.75*qp*(Ap-Bp)/p
    dUr=(Urp-Ur)/h
    Ut=(dUr-3.0*q*C)/(p+1.0)
    return A,B,C,se,q,Ur,Ut

# Build Cartesian polar field. angle th measured from x-axis.
# For theta in [0,pi] use mirror (Mode I symmetric): sigma_xx(-th)=sigma_xx(th), etc.
def sigma_cart(r,th):
    # get A,B,C for |th|
    a=abs(th)
    A,B,C,se,q,Ur,Ut=angular(a)
    R=sigma0*K*r**(-p)
    srr=A*R; stt=B*R; srt=C*R
    c=np.cos(th); sn=np.sin(th)
    sxx=c*c*srr-2*c*sn*srt+sn*sn*stt
    syy=sn*sn*srr+2*c*sn*srt+c*c*stt
    sxy=c*sn*(srr-stt)+(c*c-sn*sn)*srt
    # note: for negative th, sxy flips sign (odd), syy even -> handled by using th
    if th<0:
        sxy=-sxy
    return sxx,syy,sxy

def Wdens(r,th):
    A,B,C,se,q,Ur,Ut=angular(abs(th))
    return sigma0*alpha*eps0*K**(n+1)*r**(-1.0)*(n/(n+1.0))*se**(n+1)

def disp_cart(r,th):
    a=abs(th)
    A,B,C,se,q,Ur,Ut=angular(a)
    ur=alpha*eps0*K**n*r**p*Ur
    ut=alpha*eps0*K**n*r**p*Ut
    c=np.cos(th); sn=np.sin(th)
    ux=c*ur-sn*ut
    uy=sn*ur+c*ut
    if th<0:
        uy=-uy
    return ux,uy

def J_contour(R, N=4000):
    th=np.linspace(-np.pi,np.pi,N)
    dth=th[1]-th[0]
    integ=0.0
    # J = int (W n_x - t_i du_i/dx) ds ; ds = R dth ; n=(cos,sin)
    # t_i = sigma_ij n_j = (sxx c + sxy s , sxy c + syy s)
    # du_i/dx = d/dx u_i ; need derivatives. Use finite diff in x.
    eps=1e-5*R
    for t in th:
        sxx,syy,sxy=sigma_cart(R,t)
        c=np.cos(t); sn=np.sin(t)
        W=Wdens(R,t)
        # t_i
        tx=sxx*c+sxy*sn; ty=sxy*c+syy*sn
        # du/dx via finite differences in Cartesian x
        x=R*c; y=R*sn
        def u_at(xx,yy):
            rr=np.hypot(xx,yy); tt=np.arctan2(yy,xx)
            return disp_cart(rr,tt)
        uxp,uvp=u_at(x+eps,y); uxm,uvm=u_at(x-eps,y)
        dux=(uxp-uxm)/(2*eps); duy=(uvp-uvm)/(2*eps)
        integ+=(W*c-(tx*dux+ty*duy))*R*dth
    return integ

for R in [1e-3,1e-2,1e-1,1.0]:
    print("R=%.1e  J_contour=%.6f" % (R, J_contour(R, N=3000)))

# U(R) = int_0^R int W r dth dr
def Wbar_ang():
    th=np.linspace(-np.pi,np.pi,8001)
    val=[]
    for t in th:
        A,B,C,se,q,Ur,Ut=angular(abs(t))
        val.append((n/(n+1.0))*se**(n+1))
    return np.trapezoid(val,th)
Wb=Wbar_ang()
print("Wbar (ang integral of Wtilde) =", Wb)
print("C_pred = U/(J R) = Wbar/I_n")
# I_n from contour*R^{2p}? No: J_contour*... = sigma0 alpha eps0 K^{n+1} I_n
Jn=J_contour(1e-2,N=4000)
I_n = Jn/(sigma0*alpha*eps0*K**(n+1))
print("I_n from contour =", I_n)
A0,B0,C0f,se0,q0,Ur0,Ut0 = angular(0.0)
print("scaled B(0)=%.8f  A(0)=%.8f" % (B0, A0))
print("H = B(0)/I_n^p = %.8f" % (B0/I_n**p))
print("C = Wbar/I_n = %.8f" % (Wb/I_n))
print("=== FINAL ===  p=%.8f  H=%.8f  C=%.8f" % (p, B0/I_n**p, Wb/I_n))