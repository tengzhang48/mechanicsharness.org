#!/usr/bin/env python3
"""
MP-18 weak (Galerkin) solver.

    d_tau s = L s ,   (L s)(x) = Int_0^1 S(x,y) s'''(y) dy ,
    S(x,y) = 1/(x-y) - 1/(x+y) - 2y(x-y)/(x+y)^3 ,  (PV at y=x)
    s(0)=0 , s'(1)=s''(1)=0 ,  s(x,0)=1 ,  tau=t/t0 .

With S = d_x K ,  K(x,y) = ln|(x-y)/(x+y)| + 2xy/(x+y)^2 , and test functions that
vanish at x=1 (so the boundary term [psi K] drops),

    <psi_j , L phi_k> = - Int_0^1 Int_0^1 psi_j'(x) K(x,y) phi_k'''(y) dy dx .

Only a *logarithmic* singularity (at y=x) appears, so the double integral can be
evaluated to near machine precision by splitting each integration at the singular
line.  The Galerkin matrices are then

    G_jk = <psi_j, L phi_k> ,   M_jk = <psi_j, phi_k> ,
    and d_tau a = M^{-1} G a .

BC-satisfying trial functions:  phi_k(x) = b(x) T_k(2x-1) ,
test  functions (also vanish at 1):  psi_j(x) = (1-x) b(x) T_j(2x-1) .
"""
import numpy as np
from numpy.polynomial import polynomial as PL
from scipy.integrate import quad

B_C = np.array([0.0, 0.0, 1.0, -4.0 / 3.0, 0.5])


def cheb_coeffs(n):
    if n == 0:
        return np.array([1.0])
    tprev = np.array([1.0])
    tcur = np.array([-1.0, 2.0])
    for _ in range(1, n):
        tcur, tprev = PL.polysub(PL.polymul(np.array([-2.0, 4.0]), tcur), tprev), tcur
    return tcur


def phi_coeffs(k):
    return PL.polymul(B_C, cheb_coeffs(k))


def psi_coeffs(j):
    return PL.polymul(np.array([1.0, -1.0]), phi_coeffs(j))   # (1-x) b T_j


def K(x, y):
    return np.log(np.abs((x - y) / (x + y))) + 2.0 * x * y / (x + y) ** 2


def inner(psi_c, p3_c):
    """Int_0^1 Int_0^1 psi'(x) K(x,y) p3(y) dy dx  (log-singular, split)."""
    psip = PL.polyder(psi_c)                 # psi'(x)

    def gy(y):
        return PL.polyval(y, p3_c)

    def inner_y(x):
        # Int_0^1 K(x,y) gy(y) dy  split at y=x
        f = lambda y: K(x, y) * gy(y)
        return (quad(f, 0.0, x, limit=300, epsabs=1e-12, epsrel=1e-12)[0]
                + quad(f, x, 1.0, limit=300, epsabs=1e-12, epsrel=1e-12)[0])

    f = lambda x: PL.polyval(x, psip) * inner_y(x)
    return quad(f, 0.0, 1.0, limit=200, epsabs=1e-11, epsrel=1e-11)[0]


def mass(psi_c, phi_c):
    """Int_0^1 psi(x) phi(x) dx."""
    prod = PL.polymul(psi_c, phi_c)
    return quad(lambda x: PL.polyval(x, prod), 0.0, 1.0,
                limit=200, epsabs=1e-13, epsrel=1e-13)[0]


def assemble(nb):
    G = np.zeros((nb, nb))
    M = np.zeros((nb, nb))
    for j in range(nb):
        ps = psi_coeffs(j)
        for k in range(nb):
            p3 = PL.polyder(phi_coeffs(k), 3)
            G[j, k] = -inner(ps, p3)
            M[j, k] = mass(ps, phi_coeffs(k))
    return G, M


def run(nb):
    G, M = assemble(nb)
    A = np.linalg.solve(M, G)
    ev = np.linalg.eigvals(A)
    return ev[np.argsort(-ev.real)]


if __name__ == "__main__":
    for nb in [4, 6, 8, 10, 12, 14]:
        ev = run(nb)
        print(f"nb={nb:2d}  " + "  ".join(f"{v.real:+.6f}{v.imag:+.6f}i" for v in ev[:4]))