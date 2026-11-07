#!/usr/bin/env python3
"""
MP-18 solver: leading relaxation eigenvalues of
    d_tau s = L s ,   (L s)(x) = Int_0^1 S(x,y) s'''(y) dy ,
    S(x,y) = 1/(x-y) - 1/(x+y) - 2y(x-y)/(x+y)^3   (PV at y=x),
    s(0)=0 , s'(1)=s''(1)=0 , s(x,0)=1 ,
with t0 = 4 pi kT h^3 /(E* dgb Dgb Omega)  (so tau = t/t0 and the late-time decay is
exp(-lambda1 tau)).

Basis (satisfies the three BCs identically):  phi_n(x) = b(x) T_n(2x-1),
    b(x) = x^2 - 4/3 x^3 + 1/2 x^4      ( b(0)=0, b'(1)=b''(1)=0 ).

The operator is evaluated in the equivalent potential form
    (L s)(x) = d_x Int_0^1 K(x,y) s'''(y) dy ,
    K(x,y) = ln|(x-y)/(x+y)| + 2 x y/(x+y)^2 ,
computed with log-aware adaptive quadrature and a Chebyshev spectral x-derivative.
"""
import numpy as np
from numpy.polynomial import polynomial as PL
from scipy.integrate import quad

# ---------------------------------------------------------------- basis
def cheb_T(n, x):
    """T_n(2x-1) evaluated at array x via recurrence."""
    u = 2.0 * x - 1.0
    if n == 0:
        return np.ones_like(u)
    t0 = np.ones_like(u)
    t1 = u.copy()
    for _ in range(2, n + 1):
        t0, t1 = t1, 2.0 * u * t1 - t0
    return t1


def b_of_x(x):
    return x ** 2 - (4.0 / 3.0) * x ** 3 + 0.5 * x ** 4


def phi(n, x):
    return b_of_x(x) * cheb_T(n, x)


def phi_coeffs(n):
    """Polynomial coefficients (low->high) of phi_n."""
    tprev = np.array([1.0])
    if n == 0:
        tc = tprev
    else:
        tcur = np.array([-1.0, 2.0])
        for _ in range(1, n):
            tc = PL.polysub(PL.polymul(np.array([-2.0, 4.0]), tcur), tprev)
            tprev, tcur = tcur, tc
        tc = tcur
    b = np.array([0.0, 0.0, 1.0, -4.0 / 3.0, 0.5])
    return PL.polymul(b, tc)


# ---------------------------------------------------------------- kernel
def K(x, y):
    return np.log(np.abs((x - y) / (x + y))) + 2.0 * x * y / (x + y) ** 2


def mu_of_x(x, p3):
    """mu(x) = Int_0^1 K(x,y) g(y) dy ,  g = polynomial with coeffs p3 (low->high).
    Split at y=x so the logarithmic singularity sits at an interval endpoint, which
    scipy's adaptive quadrature resolves to full precision."""
    def integ(y):
        return K(x, y) * PL.polyval(y, p3)
    left = quad(integ, 0.0, x, limit=400, epsabs=1e-14, epsrel=1e-14)[0]
    right = quad(integ, x, 1.0, limit=400, epsabs=1e-14, epsrel=1e-14)[0]
    return left + right


def cheb_nodes(N):
    """Chebyshev-Gauss-Lobatto nodes on (0,1), interior only (i=1..N-1)."""
    i = np.arange(1, N)
    return 0.5 * (1.0 + np.cos(np.pi * i / N))


def cheb_D(N):
    """Chebyshev differentiation matrix (Lobatto) restricted to the interior nodes."""
    c = np.hstack([2.0, np.ones(N - 1), 2.0]) * ((-1.0) ** np.arange(N + 1))
    X = np.cos(np.pi * np.arange(N + 1) / N)
    dX = np.zeros((N + 1, N + 1))
    for i in range(N + 1):
        for j in range(N + 1):
            if i != j:
                dX[i, j] = c[i] / c[j] / (X[i] - X[j])
    dX = dX - np.diag(dX.sum(axis=1))
    Dx = 2.0 * dX                       # x = (1+xi)/2
    return Dx[1:N, 1:N]


def build(nb, N):
    """A (operator) and Phi (mass) matrices: nb basis functions, N-1 collocation nodes."""
    xc = cheb_nodes(N)
    Dx = cheb_D(N)
    nb = min(nb, N - 1)
    A = np.zeros((N - 1, nb))
    Phi = np.zeros((N - 1, nb))
    for n in range(nb):
        p3 = PL.polyder(phi_coeffs(n), 3)
        mu = np.array([mu_of_x(xv, p3) for xv in xc])
        A[:, n] = Dx @ mu
        Phi[:, n] = phi(n, xc)
    return xc, A, Phi


def spectrum(nb, N):
    xc, A, Phi = build(nb, N)
    M = np.linalg.lstsq(Phi, A, rcond=None)[0]
    ev, V = np.linalg.eig(M)
    order = np.argsort(-ev.real)
    return xc, A, Phi, ev[order], V[:, order]


if __name__ == "__main__":
    print(" nb  N   rightmost eigenvalues")
    for nb, N in [(8, 24), (12, 32), (16, 40), (20, 48), (24, 56)]:
        xc, A, Phi, ev, V = spectrum(nb, N)
        top = ev[:5]
        print(f"{nb:3d} {N:3d}  " + "  ".join(f"{v.real:+.5f}{v.imag:+.5f}i" for v in top))