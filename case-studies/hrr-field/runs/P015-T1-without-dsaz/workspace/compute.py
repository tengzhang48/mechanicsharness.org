"""
Compute the three numbers from the solved HRR angular functions (n=6.5).

Loads hrr_sol_6p5.npz produced by continuation.py, reconstructs
A=sigtilde_rr, B=sigtilde_tt, C=sigtilde_rt, se, and the displacement
angular functions Ur, Ut (with a common factor C0 = alpha eps0 K^n).

Then:
  I_n   = HRR normalization integral:
          I_n = int_{-pi}^{pi} [ (n/(n+1)) se^{n+1} cos(theta)
                                 - ( A * ux + C * uy ) ] dtheta
          where (ux,uy) are the Cartesian components of the dimensionless
          displacement angular function (Ur, Ut rotated).
  H     = B(0) / I_n^p
  C     = U(R)/(J(R) R)  (dimensionless, R-independent)
"""
import numpy as np
import hrr2 as H

n = 6.5
p = 1.0/(n+1.0)
k = 2.0 - p
H.n, H.p, H.k = n, p, k


def build(th, Y):
    phi, d1, d2, d3 = Y.T
    A = k*phi + d2
    B = k*(k-1.0)*phi
    C = -(k-1.0)*d1
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
    return A, B, C, se, Ur, Ut, Urp


d = np.load("hrr_sol_6p5.npz")
th = d["th"]; Y = d["Y"]
A, B, C, se, Ur, Ut, Urp = build(th, Y)
print("th range", th[0], th[-1])
print("A(0)=%.6f B(0)=%.6f C(0)=%.2e se(0)=%.6f" % (A[0], B[0], C[0], se[0]))
print("A(pi)=%.3e C(pi)=%.3e" % (A[-1], C[-1]))

# Cartesian displacement angular components:
# ux = Ur cos(th) - Ut sin(th);  uy = Ur sin(th) + Ut cos(th)
ux = Ur*np.cos(th) - Ut*np.sin(th)
uy = Ur*np.sin(th) + Ut*np.cos(th)

# Full symmetric theta in [-pi,pi]
th_full = np.concatenate([-th[::-1][:-1], th])
def mirror(x):
    return np.concatenate([x[::-1][:-1], x])
Af, Bf, Cf, sef = mirror(A), mirror(B), mirror(C), mirror(se)
uxf = mirror(ux); uyf = mirror(uy)

# I_n integral (HRR): use the standard form
integrand = (n/(n+1.0))*sef**(n+1.0)*np.cos(th_full) - (Af*uxf + Cf*uyf)
I_n = np.trapezoid(integrand, th_full)
print("I_n (with ux,uy in the -Au form) =", I_n)

# Also compute the J contour integral directly on radius R to confirm scaling
# J = int(W cos - t_i du_i/dx) R dtheta ; W = n/(n+1) sigma:eps
# In normalized terms, at radius R the factors: sigma ~ R^{-p}, eps ~ R^{-1-p}
# Compute the ANGULAR part of t_i du_i/dx and W.
# t_r = A (normalized), t_t = C.  du_i/dx involves rotating basis.
# We'll compute J_ang such that J(R) = R^{-2p} * J_ang(normalized).

# strain components normalized (divide out C0 R^{-p-1})  -- half domain [0,pi]
eps_rr = -p*Ur
eps_tt = Ur + np.gradient(Ut, th)
eps_rt = 0.5*(Urp - p*Ut - Ut)   # 2 eps_rt = Ur'-pUt-Ut  => eps_rt = half
# W_norm = n/(n+1) sigma:eps (with sigma normalized by R^{-p}, eps by R^{-p-1})
# use half-domain arrays, then mirror
Wnorm_half = (n/(n+1.0))*(A*eps_rr + B*eps_tt + 2.0*C*eps_rt)
Wnorm = mirror(Wnorm_half)
# J_ang = int(Wnorm cos - t_i du_i/dx_ang) dtheta
# du_i/dx: for the self-similar field, d/dx = (cos d/dr - sin/r d/dtheta)
# ux = R^{-p} uxbar ; dux/dR = -p R^{-p-1} uxbar
# term t_i du_i/dx with t=(A,C): t_r d ur/dx + t_t d ut/dx ; in (r,t) basis
# d ur/dx contains rotation: dur/dx = cos*ur_{,r} - sin/r*ur_{,theta} - (u_t rotation)
# Use the well-known polar form of du/dx with basis rotation:
#   dur/dx = cos ur_{,r} - (sin/r)(ur_{,theta} - u_t)
#   dut/dx = cos ut_{,r} - (sin/r)(ut_{,theta} + u_r)
# with ur_{,r} = -p R^{-p-1} Urbar, ur_{,theta}=R^{-p} Urbar'
# Cartesian derivative of displacement incl. rotation of polar basis.
# u_r = R^{-p} Urbar(th), u_t = R^{-p} Utbar(th).
# d u_r / dR = -p R^{-p-1} Urbar ;  d u_r/dth = R^{-p} Urbar'
# d u_t / dR = -p R^{-p-1} Utbar ;  d u_t/dth = R^{-p} Utbar'
# du_r/dx = cos*du_r/dR - (sin/R)*du_r/dth + (sin/R)*u_t
# du_t/dx = cos*du_t/dR - (sin/R)*du_t/dth - (sin/R)*u_r
# (the last terms come from the rotation of the basis e_r,e_t)
R = 1.0
Urb = mirror(Ur); Utb = mirror(Ut); Urpb = mirror(Urp)
dUrb = Urpb
# Utbar' : differentiate
dUtb = np.gradient(Utb, th_full)
dur_dx = np.cos(th_full)*(-p*Urb) - np.sin(th_full)*dUrb + np.sin(th_full)*Utb
dut_dx = np.cos(th_full)*(-p*Utb) - np.sin(th_full)*dUtb - np.sin(th_full)*Urb
# traction t_i = sigma_ij N_j with N = e_r: t = (A, C) in (r,t) basis
term_tt = Af*dur_dx + Cf*dut_dx
J_ang = np.trapezoid(Wnorm*np.cos(th_full) - term_tt, th_full)
W_ang = np.trapezoid(Wnorm, th_full)
print("J_ang (normalized contour integral) =", J_ang)
print("W_ang (normalized) =", W_ang)

# Check: contour J integral should be R-independent after the R^{-2p} scaling.
# U(R) = int_0^R int W R dth dR ; W_norm as defined => W = C0 sigma0 R^{-1-2p} Wnorm
# U(R) = C0 sigma0 * W_ang * R^{1-2p}/(1-2p)
# J(R) = C0 sigma0 * J_ang * R^{-2p}
# C = U(R)/(J(R) R) = [W_ang/(1-2p)] / J_ang
Cconst = (W_ang/(1.0-2.0*p))/J_ang
print("C = U(R)/(J(R)*R) =", Cconst)

# H = B(0)/I_n^p
Hval = B[0]/I_n**p
print("B(0) =", B[0], " I_n =", I_n, " H =", Hval)
print("p =", p)