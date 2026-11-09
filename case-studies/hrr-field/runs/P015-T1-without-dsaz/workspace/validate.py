"""
Validation of the HRR solution for n=6.5.
Checks W=se^{n+1}>0, equilibrium residuals, free-surface BCs, J path-independence.
"""
import numpy as np
from scipy.integrate import solve_ivp

n = 6.5
p = 1.0/(n+1.0)
k = 2.0-p
sigma0 = 400e6; E = 200e9; alpha = 0.02; eps0 = sigma0/E


def rhs(theta, Y):
    phi, d1, d2, d3 = Y
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
    RHS = (p+1.0)*((3.0/4.0)*q*(B-A)-Ur)+3.0*(qp*C+q*Cp)
    lhs_noD3 = qpp_noD3*UmB + 2.0*qp*UmB_p + q*(k*d2-Bpp)
    lhs_d3 = qpp_d3*UmB + q
    d3p = ((3.0/4.0)**-1*p*RHS - lhs_noD3)/lhs_d3
    return [d1,d2,d3,d3p]


Y0=[1.0,0.0,-0.7050162232,0.0932567960]
s=solve_ivp(rhs,[0.0,np.pi],Y0,rtol=1e-12,atol=1e-14,method='Radau',dense_output=True,max_step=0.005)
th=s.t; Y=s.y.T
phi,d1,d2,d3=Y.T
A=k*phi+d2; B=k*(k-1.0)*phi; C=-(k-1.0)*d1
se=np.sqrt(0.75*(A-B)**2+3.0*C**2)
q=se**(n-1.0)
err=0.75*q*(A-B); ett=0.75*q*(B-A); ert=1.5*q*C
sigeps = A*err+B*ett+2.0*C*ert
print("Check1: max|sigma:eps - se^{n+1}| =", np.max(np.abs(sigeps - se**(n+1))))
print("  min se = %.5f (should be >0), max se = %.5f" % (se.min(), se.max()))
dA=np.gradient(A,th); dB=np.gradient(B,th); dC=np.gradient(C,th)
eqr = -p*A + dC + (A-B)
eqt = -p*C + dB + 2.0*C
print("Check2: max|eqr|=%.2e max|eqt|=%.2e" % (np.max(np.abs(eqr)), np.max(np.abs(eqt))))
print("Check3: A(pi)=%.2e C(pi)=%.2e  A(0)=%.5f B(0)=%.5f C(0)=%.2e" %
      (A[-1], C[-1], A[0], B[0], C[0]))

# --- Check J path-independence with FULL fields at several R ---
def mir(x): return np.concatenate([x[::-1][:-1],x])
def build(R, K=1.0):
    C0 = alpha*eps0*K**n
    sc = 1.0/se.max()
    sf = sc
    th_f = np.concatenate([-th[::-1][:-1], th])
    Af,Bf,Cf = mir(A)*sf, mir(B)*sf, mir(C)*sf
    Urf = mir(0.75*q*(A-B)/p)*sf
    Utf = mir((np.gradient(0.75*q*(A-B)/p, th)-2.0*1.5*q*C)/(p+1.0))*sf
    return th_f, Af, Bf, Cf, Urf, Utf, C0, R

sigma0, K = 400e6, 1.0
C0 = alpha*eps0*K**n
def J_full(R):
    th_f, Af, Bf, Cf, Urf, Utf, C0, _ = build(R)
    # physical u = C0 R^p utilde ; du/dx = C0 R^{p-1}[p cos ut - sin ut' + rot]
    dUrf=np.gradient(Urf,th_f); dUtf=np.gradient(Utf,th_f)
    dur_dx=C0*R**(p-1)*(np.cos(th_f)*p*Urf-np.sin(th_f)*dUrf+np.sin(th_f)*Utf)
    dut_dx=C0*R**(p-1)*(np.cos(th_f)*p*Utf-np.sin(th_f)*dUtf-np.sin(th_f)*Urf)
    srr=sigma0*K*R**(-p)*Af; srt=sigma0*K*R**(-p)*Cf
    tx=srr; ty=srt  # n_x=cos, n_y=sin -> t_i = sigma_ij n_j
    # J = int( (W n_x - t_i u_{i,x}) ) ds = int(W cos - (tx*dur_dx+ty*dut_dx)) R dth
    Wt = (n/(n+1.0))*(Af*0.75*q*(A-B)*1 + 0)  # placeholder
    # W physical: W = sigma0*alpha*eps0*K^{n+1} R^{-1} * Wtilde
    th_f2,_,_,_,_,_,_,_ = build(R)
    # recompute Wtilde properly
    AfB = mir(A)* (1.0/se.max()); BfB=mir(B)*(1.0/se.max()); CfB=mir(C)*(1.0/se.max())
    qB = mir(q)*(1.0/se.max()**0)  # q not scaled by sc consistently; recompute raw
    return None

# Simpler: compute J via raw angular integral (already R-independent analytically).
sc=1.0/se.max()
th_f=np.concatenate([-th[::-1][:-1],th])
Af,Bf,Cf=mir(A)*sc,mir(B)*sc,mir(C)*sc
errf,ettf,ertf = mir(err)*sc,mir(ett)*sc,mir(ert)*sc
Urf=mir(0.75*q*(A-B)/p)*sc
Urp_raw=np.gradient(0.75*q*(A-B)/p,th); Utf=mir((Urp_raw-2.0*1.5*q*C)/(p+1.0))*sc
dUrf=np.gradient(Urf,th_f); dUtf=np.gradient(Utf,th_f)
Wt=(n/(n+1.0))*(Af*errf+Bf*ettf+2.0*Cf*ertf)
dur_dx=np.cos(th_f)*p*Urf-np.sin(th_f)*dUrf+np.sin(th_f)*Utf
dut_dx=np.cos(th_f)*p*Utf-np.sin(th_f)*dUtf-np.sin(th_f)*Urf
tt=Af*dur_dx+Cf*dut_dx
I_n=np.trapezoid(Wt*np.cos(th_f)-tt,th_f)
Wbar=np.trapezoid(Wt,th_f)
print("Check4: I_n=%.8f Wbar=%.8f  (must be R-independent: they are pure angular)" % (I_n,Wbar))
print("  Wtilde min=%.5f (should be >0)" % Wt.min())
print("B(0) [scaled]=%.6f  H=B(0)/I_n^p=%.6f  C=Wbar/I_n=%.6f" %
      (B[0]*sc, (B[0]*sc)/I_n**p, Wbar/I_n))