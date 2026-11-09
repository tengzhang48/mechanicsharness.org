"""
Clean computation of I_n, H, C.

Normalized self-similar field:
  sigma_ij = sigma0 * K * r^{-p} * sigtilde_ij(theta)
  u_i      = alpha*eps0*K^n * r^{-p} * utilde_i(theta)
  eps_ij   = alpha*eps0*K^n * r^{-p-1} * epstilde_ij(theta)

J-contour integral at radius R (per unit thickness):
  J = int_{-pi}^{pi} ( W cos(th) - t_i du_i/dx ) R dth.
With sigma ~ R^{-p}, eps ~ R^{-p-1}:
  W = n/(n+1) sigma:eps = (sigma0 alpha eps0 K^{n+1}) R^{-1-2p} * Wtilde
  t_i du_i/dx = (sigma0 alpha eps0 K^{n+1}) R^{-2p} * (ttilde)
  => J(R) = (sigma0 alpha eps0 K^{n+1}) R^{-2p} * I_n,
     I_n = int ( Wtilde cos - ttilde ) dth     [R-independent angular integral]

U(R) = int_0^R int W R dth dR
     = (sigma0 alpha eps0 K^{n+1}) * Wbar * R^{1-2p}/(1-2p),
     Wbar = int Wtilde dth.

=> C = U(R)/(J(R) R) = [Wbar/(1-2p)] / I_n.

H: sigma_tt(r,0) = sigma0 K r^{-p} B(0) ; J = sigma0 alpha eps0 K^{n+1} I_n
   sigma0 [J/(alpha eps0 sigma0 r)]^p = sigma0 (I_n K^{n+1} r^{-1})^p
       = sigma0 I_n^p K r^{-p}   (using p(n+1)=1)
   => H = B(0)/I_n^p.
"""
import numpy as np
import hrr2 as H

n = 6.5
p = 1.0/(n+1.0)
k = 2.0 - p
H.n, H.p, H.k = n, p, k

d = np.load("hrr_sol_6p5.npz")
th = d["th"]; Y = d["Y"]
phi, d1, d2, d3 = Y.T
A = k*phi + d2                 # sigtilde_rr
B = k*(k-1.0)*phi              # sigtilde_tt
C = -(k-1.0)*d1                # sigtilde_rt
se = np.sqrt(0.75*(A-B)**2 + 3.0*C**2)
Ap = k*d1 + d3
Bp = k*(k-1.0)*d1
Cp = -(k-1.0)*d2
sep = (0.75*2.0*(A-B)*(Ap-Bp) + 3.0*2.0*C*Cp)/(2.0*se)
q = se**(n-1.0)
qp = (n-1.0)*se**(n-2.0)*sep
Urp = -(1.0/(2.0*p))*(qp*(A-B) + q*(Ap-Bp))
Ur = -q*(A-B)/(2.0*p)
Ut = (Urp - 2.0*q*C)/(p+1.0)

# normalize so max se = 1 (SCALE only; H, C are scale-invariant)
smax = se.max()
scale = 1.0/smax
A, B, C, se = A*scale, B*scale, C*scale, se*scale
Ur, Ut, Urp = Ur*scale, Ut*scale, Urp*scale

# Full domain [-pi, pi]
th_f = np.concatenate([-th[::-1][:-1], th])
def mir(x):
    return np.concatenate([x[::-1][:-1], x])
Af, Bf, Cf, sef = mir(A), mir(B), mir(C), mir(se)
Urf, Utf, Urpf = mir(Ur), mir(Ut), mir(Urp)
# theta-derivatives on full domain
dUrf = Urpf
dUtf = np.gradient(Utf, th_f)

# strain angular functions (epstilde)
err = -p*Urf
ett = Urf + dUtf
ert = 0.5*(Urpf - p*Utf - Utf)
Wt = (n/(n+1.0))*(Af*err + Bf*ett + 2.0*Cf*ert)   # Wtilde

# Cartesian displacement + its derivative incl basis rotation
ux = Urf*np.cos(th_f) - Utf*np.sin(th_f)
uy = Urf*np.sin(th_f) + Utf*np.cos(th_f)
# d/dx = cos d/dR - (sin/R) d/dth, with basis-rotation corrections
dur_dx = np.cos(th_f)*(-p*Urf) - np.sin(th_f)*dUrf + np.sin(th_f)*Utf
dut_dx = np.cos(th_f)*(-p*Utf) - np.sin(th_f)*dUtf - np.sin(th_f)*Urf
# traction t = sigma . e_r = (A, C) in (r,t) basis ; t_i du_i/dx = A*dur_dx + C*dut_dx
tt = Af*dur_dx + Cf*dut_dx

I_n = np.trapezoid(Wt*np.cos(th_f) - tt, th_f)
Wbar = np.trapezoid(Wt, th_f)
print("se_max (before scaling was %.5f)" % smax)
print("A(0)=%.6f  B(0)=%.6f  se(0)=%.6f" % (A[0]/scale, B[0]/scale, se[0]/scale))
print("B(0) [normalized max se=1] =", B[0])
print("I_n =", I_n)
print("Wbar =", Wbar)
Cval = (Wbar/(1.0-2.0*p))/I_n
print("C = U(R)/(J(R)*R) =", Cval)
Hval = B[0]/I_n**p
print("H =", Hval)