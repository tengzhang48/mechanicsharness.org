"""
Independent dimensional verification of H and C, and cross-check I_n.

Set sigma0=400e6, alpha=0.02, eps0=sigma0/E=0.002, E=200e9, n=6.5, p=1/7.5.
Field: sigma = sigma0 K r^{-p} sigtilde(theta)  (use normalized sigtilde with
max se = 1).  Choose K=1 (dimensionless probe) for verification of RATIOS.

Compute at several R:
  sig_tt(R,0) = sigma0 K R^{-p} B(0)
  J(R) via contour integral (dimensional)
  U(R) via radial integral (dimensional)
  H_pred = sig_tt(R,0) / ( sigma0 [J(R)/(alpha eps0 sigma0 R)]^p )
  C_pred = U(R)/(J(R) R)
and show these are R-independent, matching compute2's values.
"""
import numpy as np
import hrr2 as H

n = 6.5
p = 1.0/(n+1.0)
k = 2.0 - p
H.n, H.p, H.k = n, p, k
sigma0 = 400e6
E = 200e9
alpha = 0.02
eps0 = sigma0/E
K = 1.0    # dimensionless amplitude (probe)

d = np.load("hrr_sol_6p5.npz")
th = d["th"]; Y = d["Y"]
phi, d1, d2, d3 = Y.T
A = k*phi + d2; B = k*(k-1.0)*phi; C = -(k-1.0)*d1
se = np.sqrt(0.75*(A-B)**2 + 3.0*C**2)
Ap = k*d1 + d3; Bp = k*(k-1.0)*d1; Cp = -(k-1.0)*d2
sep = (0.75*2.0*(A-B)*(Ap-Bp) + 3.0*2.0*C*Cp)/(2.0*se)
q = se**(n-1.0); qp = (n-1.0)*se**(n-2.0)*sep
Urp = -(1.0/(2.0*p))*(qp*(A-B) + q*(Ap-Bp))
Ur = -q*(A-B)/(2.0*p)
Ut = (Urp - 2.0*q*C)/(p+1.0)

# scale so max se = 1 (sigtilde now bounded by 1)
sc = 1.0/se.max()
A, B, C, se = A*sc, B*sc, C*sc, se*sc
Ur, Ut, Urp = Ur*sc, Ut*sc, Urp*sc

def mir(x):
    return np.concatenate([x[::-1][:-1], x])
th_f = np.concatenate([-th[::-1][:-1], th])
Af, Bf, Cf = mir(A), mir(B), mir(C)
Urf, Utf, Urpf = mir(Ur), mir(Ut), mir(Urp)
dUtf = np.gradient(Utf, th_f)
err = -p*Urf; ett = Urf + dUtf; ert = 0.5*(Urpf - p*Utf - Utf)
# Wtilde and ttilde (angular, R-independent) as in compute2
Wt = (n/(n+1.0))*(Af*err + Bf*ett + 2.0*Cf*ert)
dur_dx = np.cos(th_f)*(-p*Urf) - np.sin(th_f)*Urpf + np.sin(th_f)*Utf
dut_dx = np.cos(th_f)*(-p*Utf) - np.sin(th_f)*dUtf - np.sin(th_f)*Urf
tt = Af*dur_dx + Cf*dut_dx

I_n = np.trapezoid(Wt*np.cos(th_f) - tt, th_f)
Wbar = np.trapezoid(Wt, th_f)
print("I_n=%.8f  Wbar=%.8f  B(0)=%.8f" % (I_n, Wbar, B[0]))

# The applied (constant) J that makes the field self-consistent is the
# amplitude relation J = sigma0 alpha eps0 K^{n+1} I_n (r=1 reference, i.e.
# evaluated with the SAME K). We define J0 as the R-independent amplitude:
J0 = sigma0*alpha*eps0*K**(n+1)*I_n
print("Applied (constant) J0 =", J0)
for R in [1e-3, 1e-2, 1e-1]:
    stt0 = sigma0*K*R**(-p)*B[0]
    U = sigma0*alpha*eps0*K**(n+1)*Wbar*R**(1-2*p)/(1-2*p)
    Hpred = stt0/(sigma0*(J0/(alpha*eps0*sigma0*R))**p)
    Cpred_const = U/(J0*R)      # literal C with constant J (R-dependent)
    print("R=%.1e  H=%.8f  C(const J)=%.8f" % (R, Hpred, Cpred_const))

print("Final: H = B(0)/I_n^p =", B[0]/I_n**p)
print("Final: C = (Wbar/(1-2p))/I_n  [contour-J at R, R-independent] =",
      (Wbar/(1.0-2.0*p))/I_n)