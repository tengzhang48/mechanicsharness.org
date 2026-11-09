"""
CORRECTED HRR plane-strain Mode I, Airy form (factor 3/2 included).

Phi = r^k phi, k=2-p, p=1/(n+1).
  A = sigtilde_rr = k phi + phi''
  B = sigtilde_tt = k(k-1) phi
  C = sigtilde_rt = -(k-1) phi'
  se^2 = (3/4)(A-B)^2 + 3 C^2

u_i = r^p U_i(theta); strain ~ r^{p-1}.
Constitutive (eps=(3/2)alpha eps0 (se/sigma0)^{n-1} s_ij/sigma0), plane strain:
  (I)  p U_r            = (3/4) q (A-B)
  (II) U_r + U_t'       = (3/4) q (B-A)
  (III) U_r' - p U_t - U_t = 3 q C
with q = se^{n-1}.

From (I):  U_r = (3/4) q (A-B)/p
From (III): U_t = (U_r' - 3 q C)/(p+1)
From (II): U_t' = (3/4) q (B-A) - U_r
Equate U_t' from (II) and from differentiating (III):
   [U_r'' - 3(q' C + q C')]/(p+1) = (3/4)q(B-A) - U_r
   => U_r'' = 3(q' C + q C') + (p+1)[(3/4) q (B-A) - U_r]
And U_r'' = (3/4)[q''(A-B) + 2 q'(A-B)' + q (A-B)'']/p
with (A-B)'' = k p phi'' + phi''''.
"""
import numpy as np

N = 6.5
P = 1.0/(N+1.0)
K = 2.0 - P


def quantities(phi, d1, d2, k=K):
    A = k*phi + d2
    B = k*(k-1.0)*phi
    C = -(k-1.0)*d1
    se = np.sqrt(0.75*(A-B)**2 + 3.0*C**2)
    return A, B, C, se


def rhs(theta, Y, n=N, p=P, k=K):
    phi, d1, d2, d3 = Y
    A, B, C, se = quantities(phi, d1, d2, k)
    se = np.maximum(se, 1e-30)
    Ap = k*d1 + d3
    Bp = k*(k-1.0)*d1
    Cp = -(k-1.0)*d2
    UmB = A - B
    UmB_p = Ap - Bp
    sep = (0.75*UmB*UmB_p + 3.0*C*Cp)/se
    q = se**(n-1.0)
    qp = (n-1.0)*se**(n-2.0)*sep
    Ur = 0.75*q*UmB/p
    # (III): U_r' + (p-1)U_t = 3 q C  => U_t=(3qC-U_r')/(p-1)
    # (II):  U_t' = (3/4)q(B-A) - U_r
    # differentiate (III): (p-1)U_t' = 3(q'C+qC') - U_r''
    # equate: U_r'' = 3(q'C+qC') - (p-1)[(3/4)q(B-A) - U_r]
    Urpp = 3.0*(qp*C + q*Cp) - (p-1.0)*(0.75*q*(B-A) - Ur)
    Bpp = k*(k-1.0)*d2
    Cpp = -(k-1.0)*d3
    S0 = 1.5*(UmB_p**2 + UmB*(k*p*d2 - Bpp)) + 6.0*(Cp**2 + C*Cpp)
    sepp_0 = (S0 - 2.0*sep**2)/(2.0*se)
    sepp_d4 = (1.5*UmB)/(2.0*se)
    qpp_0 = (n-1.0)*(n-2.0)*se**(n-3.0)*sep**2 + (n-1.0)*se**(n-2.0)*sepp_0
    qpp_d4 = (n-1.0)*se**(n-2.0)*sepp_d4
    d4 = (p*Urpp - 0.75*(qpp_0*UmB + 2.0*qp*UmB_p + q*k*p*d2))/(0.75*(qpp_d4*UmB + q))
    return np.stack([np.broadcast_to(np.asarray(x, dtype=float), np.asarray(d1).shape)
                     for x in (d1, d2, d3, d4)])