"""Spectral eigen-solver for  lambda f = -L[f''']  with
f(0)=0, f'(1)=f''(1)=0.

Basis (all satisfy the 3 BCs):
  phi_0(Z) = 1 - (1-Z)^(5/2)
  phi_j(Z) = Z (1-Z)^(5/2) T_{j-1}(2Z-1),  j = 1..N
Collocation at Z_i = (1-cos(pi i/(N+1)))/2, i=1..N+1 (Chebyshev-Lobatto on (0,1]).

Generalized problem:  B c = lam A c,  A_ij = phi_j(Z_i),  B_ij = -L[phi_j'''](Z_i).
"""
import numpy as np
from scipy.linalg import eig
from kernel import LOperator, chebT

POW = 2.5


def Pderivs(j, Zp):
    """P(Z)=Z*T_j(2Z-1) and its Z-derivatives up to 4th, at points Zp."""
    x = 2.0*Zp - 1.0
    T0 = chebT(j, x, 0); T1 = chebT(j, x, 1); T2 = chebT(j, x, 2)
    T3 = chebT(j, x, 3); T4 = chebT(j, x, 4)
    P = Zp*T0
    P1 = T0 + 2.0*Zp*T1
    P2 = 4.0*T1 + 4.0*Zp*T2
    P3 = 12.0*T2 + 8.0*Zp*T3
    P4 = 32.0*T3 + 16.0*Zp*T4
    return P, P1, P2, P3, P4


def make_Gfun(j):
    """Return G(y)=2y phi_j'''(1-y^2) and its derivative G'(y)."""
    if j == 0:
        # phi_0 = 1-(1-Z)^2.5 ; phi_0''' = (15/8)(1-Z)^(-1/2) ; G = 2y*(15/8)/y = 15/4
        return (lambda y: 3.75*np.ones_like(np.asarray(y, dtype=float)),
                lambda y: np.zeros_like(np.asarray(y, dtype=float)))
    def Gfun(y):
        y = np.asarray(y, dtype=float)
        Zp = 1.0 - y**2
        P, P1, P2, P3, P4 = Pderivs(j-1, Zp)
        return (2.0*y**6*P3 - 15.0*y**4*P2 + 22.5*y**2*P1 - 3.75*P)
    def Gpfun(y):
        y = np.asarray(y, dtype=float)
        Zp = 1.0 - y**2
        P, P1, P2, P3, P4 = Pderivs(j-1, Zp)
        return (-4.0*y**7*P4 + 42.0*y**5*P3 - 105.0*y**3*P2 + 52.5*y*P1)
    return Gfun, Gpfun


def phi_vals(j, Z):
    Z = np.asarray(Z, dtype=float)
    x = (1.0 - Z)**POW
    if j == 0:
        return 1.0 - x
    return Z*x*chebT(j-1, 2.0*Z-1.0, 0)


def build_matrices(N, M=512):
    Lop = LOperator(M)
    i = np.arange(1, N+2)
    Zc = 0.5*(1.0 - np.cos(np.pi*i/(N+1.0)))
    A = np.zeros((N+1, N+1))
    B = np.zeros((N+1, N+1))
    for j in range(N+1):
        A[:, j] = phi_vals(j, Zc)
        Gfun, Gpfun = make_Gfun(j)
        B[:, j] = -Lop.apply(Zc, Gfun, Gpfun)
    return A, B, Zc


def spectrum(N, M=512):
    A, B, Zc = build_matrices(N, M)
    lam, vec = eig(B, A)
    # filter finite, positive-real
    lam = np.asarray(lam)
    ok = np.isfinite(lam)
    lam = lam[ok]
    vec = vec[:, ok]
    # take eigenvalues with small imaginary part and positive real part
    real_ok = (np.abs(lam.imag) < 1e-6*np.maximum(1.0, np.abs(lam.real))) & (lam.real > 1e-9)
    lamr = lam[real_ok]; vecr = vec[:, real_ok]
    order = np.argsort(lamr.real)
    return lamr[order], vecr[:, order], A, B, Zc


def mode_value(vec, lam_order, Zpts):
    """f(Zpts) for the eigenvector at index lam_order."""
    c = vec[:, lam_order].real
    c = c/np.linalg.norm(c)
    N = len(c) - 1
    out = np.zeros(len(Zpts))
    for j in range(N+1):
        out += c[j]*phi_vals(j, Zpts)
    return out, c


if __name__ == "__main__":
    import sys
    Ns = [10, 14, 18, 22, 26, 30, 36, 44, 52]
    print(f"{'N':>4} {'lam1':>14} {'lam2':>14} {'f(0.5)/f(1)':>14}")
    for N in Ns:
        lam, vec, A, B, Zc = spectrum(N, M=384)
        l1, l2 = lam[0].real, lam[1].real
        # mode1 shape: evaluate at 0.5 and 1
        f05, c = mode_value(vec, 0, np.array([0.5]))
        f1 = c[0]  # phi_0(1)=1, others 0
        print(f"{N:4d} {l1:14.8f} {l2:14.8f} {(f05[0]/f1):14.8f}")
