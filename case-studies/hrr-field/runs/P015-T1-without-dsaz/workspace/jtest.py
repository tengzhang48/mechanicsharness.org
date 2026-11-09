"""
Direct test: is the contour J-integral at radius R R-independent for the leading
HRR self-similar field?
  J(R) = int_{-pi}^{pi} ( W cos(th) - t_i du_i/dx ) R dth.
Dimensional: sigma=sigma0 K R^{-p} sigtilde, u=alpha eps0 K^n R^{-p} utilde.
Predicted scaling J(R) ~ R^{-2p}.
"""
import numpy as np
import hrr2 as H

n = 6.5; p = 1.0/(n+1.0); k = 2.0-p
H.n, H.p, H.k = n, p, k
sigma0 = 400e6; E = 200e9; alpha = 0.02; eps0 = sigma0/E; K = 1.0

d = np.load("hrr_sol_6p5.npz")
th = d["th"]; Y = d["Y"]
phi, d1, d2, d3 = Y.T
A = k*phi + d2; B = k*(k-1.0)*phi; C = -(k-1.0)*d1
se = np.sqrt(0.75*(A-B)**2 + 3.0*C**2)
Ap = k*d1+d3; Bp = k*(k-1.0)*d1; Cp = -(k-1.0)*d2
sep = (0.75*2.0*(A-B)*(Ap-Bp)+3.0*2.0*C*Cp)/(2.0*se)
q = se**(n-1.0); qp=(n-1.0)*se**(n-2.0)*sep
Urp = -(1.0/(2.0*p))*(qp*(A-B)+q*(Ap-Bp))
Ur = -q*(A-B)/(2.0*p); Ut = (Urp-2.0*q*C)/(p+1.0)
sc = 1.0/se.max()
A,B,C,se = A*sc,B*sc,C*sc,se*sc
Ur,Ut,Urp = Ur*sc,Ut*sc,Urp*sc
def mir(x): return np.concatenate([x[::-1][:-1], x])
th_f = np.concatenate([-th[::-1][:-1], th])
Af,Bf,Cf = mir(A),mir(B),mir(C)
Urf,Utf,Urpf = mir(Ur),mir(Ut),mir(Urp)
dUtf = np.gradient(Utf, th_f)
err=-p*Urf; ett=Urf+dUtf; ert=0.5*(Urpf-p*Utf-Utf)
Wt = (n/(n+1.0))*(Af*err+Bf*ett+2.0*Cf*ert)
dur_dx = np.cos(th_f)*(-p*Urf)-np.sin(th_f)*Urpf+np.sin(th_f)*Utf
dut_dx = np.cos(th_f)*(-p*Utf)-np.sin(th_f)*dUtf-np.sin(th_f)*Urf
tt = Af*dur_dx+Cf*dut_dx
Iang = np.trapezoid(Wt*np.cos(th_f)-tt, th_f)

C0 = sigma0*alpha*eps0*K**(n+1)
for R in [1e-3, 1e-2, 1e-1]:
    J = C0*R**(-2*p)*Iang
    print("R=%.1e  J(R)=%.6e   J*R^{2p}=%.6e" % (R, J, J*R**(2*p)))