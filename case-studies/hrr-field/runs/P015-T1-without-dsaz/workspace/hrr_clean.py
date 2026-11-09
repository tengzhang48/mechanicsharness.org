"""
HRR plane-strain Mode I, Airy form. Corrected derivation. Numeric RHS.

Phi = r^k phi(theta), k=2-p, p=1/(n+1).
  A = sigtilde_rr = k phi + phi''
  B = sigtilde_tt = k(k-1) phi
  C = sigtilde_rt = -(k-1) phi'

Plane strain, s_zz=0:
  s_rr=(A-B)/2, s_tt=(B-A)/2, s_rt=C
  se^2 = (3/4)(A-B)^2 + 3 C^2

Displacement u_i = r^{p} U_i(theta) (p=1/(n+1)), strain ~ r^{p-1}:
  e_rr  = p Ur
  e_tt  = Ur + Ut'
  e_rt  = (Ur' + (p-1) Ut)/2
Constitutive (scale free; q=se^{n-1}>0):
  p Ur            = q (A-B)/2        (I)
  Ur + Ut'        = q (B-A)/2        (II)
  Ur' + (p-1) Ut  = 2 q C            (III)

=> Ur  = q(A-B)/(2p)
   Ut  = (2 q C - Ur')/(p-1)               from (III)
   Ut' = q(B-A)/2 - Ur                     from (II)
   Ut' = d/dtheta[(2 q C - Ur')/(p-1)] = [2(q' C + q C') - Ur'']/(p-1)
Equate:  Ur'' = 2(q' C + q C') - (p-1)(q(B-A)/2 - Ur)
And Ur'' = [q'' (A-B) + 2 q' (A-B)' + q (A-B)'']/(2p), gives phi''''.
  (A-B)'' = (k - k(k-1)) phi'' + phi'''' = k p phi'' + phi''''  (since 2-k=p).
"""
import numpy as np

N = 6.5
P = 1.0 / (N + 1.0)
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
    Ur = q*UmB/(2.0*p)
    Urpp = 2.0*(qp*C + q*Cp) - (p-1.0)*(q*(B-A)/2.0 - Ur)
    Bpp = k*(k-1.0)*d2
    Cpp = -(k-1.0)*d3
    S0 = 1.5*(UmB_p**2 + UmB*(k*p*d2 - Bpp)) + 6.0*(Cp**2 + C*Cpp)
    sepp_0 = (S0 - 2.0*sep**2)/(2.0*se)
    sepp_d4 = (1.5*UmB)/(2.0*se)
    qpp_0 = (n-1.0)*(n-2.0)*se**(n-3.0)*sep**2 + (n-1.0)*se**(n-2.0)*sepp_0
    qpp_d4 = (n-1.0)*se**(n-2.0)*sepp_d4
    d4 = (2.0*p*Urpp - qpp_0*UmB - 2.0*qp*UmB_p - q*k*p*d2)/(qpp_d4*UmB + q)
    return np.stack([np.broadcast_to(np.asarray(x, dtype=float), np.asarray(d1).shape)
                     for x in (d1, d2, d3, d4)])


def bc(ya, yb):
    # theta=0 Mode I symmetry: C(0)=-(k-1)phi'(0)=0 => phi'(0)=0;
    #   phi'''(0)=0 (even symmetry of phi about 0).
    # normalization phi(0)=1.
    # theta=pi traction-free: A(pi)=0 (k phi+phi''=0).
    return np.array([ya[1], ya[3], ya[0]-1.0, K*yb[0]+yb[2]])