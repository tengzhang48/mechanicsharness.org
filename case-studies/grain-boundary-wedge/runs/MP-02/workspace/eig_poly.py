"""Independent check: polynomial basis (no fractional tip power).
  phi_0(Z) = 1-(1-Z)^3  (=3Z-3Z^2+Z^3), phi_0''' = 6
  phi_j(Z) = Z(1-Z)^3 T_{j-1}(2Z-1), j = 1..N   (polynomial, smooth)
All satisfy phi(0)=0, phi'(1)=phi''(1)=0.
In x=2Z-1: Z(1-Z)^3 = (1+x)(1-x)^3/16 = (5T0-4T1-4T2+4T3-T4)/128.
Expect slower (algebraic) convergence to the same eigenvalues.
"""
import numpy as np
from numpy.polynomial import chebyshev as C
from scipy.linalg import eig
from kernel import LOperator

MULT = np.zeros(5)
MULT[0] = 5.0/128.0; MULT[1] = -4.0/128.0; MULT[2] = -4.0/128.0
MULT[3] = 4.0/128.0; MULT[4] = -1.0/128.0   # (1+x)(1-x)^3/16 in Chebyshev-x


def phi_cheb(j):
    """Chebyshev-x coefficients of phi_j(Z(x))."""
    if j == 0:
        # 1 - (1-x)^3/8 ; (1-x)^3 = 1-3x+3x^2-x^3
        # = T0 -3T1 +3(T2+T0)/2?  x^2=(T2+T0)/2? wait T2=2x^2-1 -> x^2=(T2+1)/2
        # x^3 = (T3+3T1)/4
        # (1-x)^3 = T0 -3 T1 + 3*(T2+T0)/2 - (T3+3T1)/4
        c = np.zeros(4)
        c[0] = 1.0 + 1.5          # T0 coeff: 1 + 3/2
        c[1] = -3.0 - 0.75        # -3 - 3/4
        c[2] = 1.5
        c[3] = -0.25
        phi = np.zeros(4)
        phi[0] = 1.0
        return phi - c/8.0
    base = np.zeros(j); base[j-1] = 1.0   # T_{j-1}
    return C.chebmul(MULT, base)


def make_Gfun_poly(j):
    ch = phi_cheb(j)
    c3 = ch.copy()
    for _ in range(3):
        c3 = C.chebder(c3)
    c3 = 8.0*c3                     # d^3/dZ^3 = 8 d^3/dx^3
    c4 = C.chebder(c3)
    c4 = 2.0*c4                     # d^4/dZ^4 = 16 d^4/dx^4; relative to c3: *2
    def Gfun(y):
        y = np.asarray(y, dtype=float)
        return 2.0*y*C.chebval(1.0-2.0*y**2, c3)
    def Gpfun(y):
        y = np.asarray(y, dtype=float)
        return 2.0*C.chebval(1.0-2.0*y**2, c3) - 4.0*y**2*C.chebval(1.0-2.0*y**2, c4)
    return Gfun, Gpfun


def phi_vals_poly(j, Z):
    x = 2.0*np.asarray(Z, dtype=float)-1.0
    return C.chebval(x, phi_cheb(j))


def spectrum_poly(N, M=384):
    Lop = LOperator(M)
    # interior Chebyshev roots on (0,1): exclude Z=1 because L[phi_j'''] can
    # diverge there for polynomial phi (phi'''(1) != 0).
    i = np.arange(1, N+2)
    Zc = 0.5*(1.0 - np.cos(np.pi*(2.0*i-1.0)/(2.0*(N+1.0))))
    A = np.zeros((N+1, N+1)); B = np.zeros((N+1, N+1))
    for j in range(N+1):
        A[:, j] = phi_vals_poly(j, Zc)
        Gf, Gp = make_Gfun_poly(j)
        B[:, j] = -Lop.apply(Zc, Gf, Gp)
    lam, vec = eig(B, A)
    lam = lam[np.isfinite(lam)]
    ok = (np.abs(lam.imag) < 1e-6*np.maximum(1.0, np.abs(lam.real))) & (lam.real > 1e-9)
    lam = lam[ok]
    return np.sort(lam.real), vec, A, B


if __name__ == "__main__":
    print(f"{'N':>4} {'lam1':>14} {'lam2':>14}")
    for N in [10, 16, 24, 32, 44, 60, 80]:
        lam, vec, A, B = spectrum_poly(N)
        print(f"{N:4d} {lam[0]:14.8f} {lam[1]:14.8f}")
