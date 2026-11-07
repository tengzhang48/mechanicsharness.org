"""Eigenanalysis via g = f''' formulation with built-in BCs and tip singularity.

f(x) = C1 x + C2 x^2 + int_0^x (x-t)^2/2 g(t) dt, with f'(1)=f''(1)=0 fixing
C1 = int t g, C2 = -int g /2.  Basis: g(xi) = xi (1-xi)^{-1/2} T_m(2xi-1).
Pencil: A[g] = -lambda f  ->  eig(A, -F).
All quadrature via t-substitution xi = 1 - t^2 (smooth integrands).
"""
import numpy as np
from numpy.polynomial.legendre import leggauss
from gbcore import Rfun


def tquad(a, b, M):
    """GL nodes/weights on [a,b] mapped to t-space."""
    t, w = leggauss(M)
    return 0.5 * (b - a) * t + 0.5 * (a + b), 0.5 * (b - a) * w


def eigen_g(Nv=48, Mt=400, verbose=True, nkeep=6):
    """Nv: number of Chebyshev basis functions for v (g = xi(1-xi)^-1/2 v)."""
    # basis functions evaluated at points xi: B_m(xi)
    def Tbasis(xi, mmax):
        return np.cos(np.outer(np.arccos(np.clip(2 * xi - 1, -1, 1)), np.arange(mmax + 1)))

    # collocation points (interior)
    k = np.arange(1, Nv + 2)
    xc = 0.5 * (1 - np.cos(np.pi * k / (Nv + 2)))  # Nv+1 points strictly in (0,1)

    A = np.zeros((Nv + 1, Nv + 1))
    F = np.zeros((Nv + 1, Nv + 1))

    for i, x in enumerate(xc):
        # ---- A[g](x) = PV int Shat g dxi + int R g dxi, in t-space xi=1-t^2
        t, wt = tquad(0.0, 1.0, Mt)
        xi = 1 - t**2
        Bq = Tbasis(xi, Nv)               # (Mt, Nv+1) basis v at xi
        Gq = xi[:, None] / t[:, None] * Bq  # g(xi) = xi/t * v  (t>0)
        Gx = x / np.sqrt(1 - x) * 1.0     # placeholder; recompute per basis below
        # PV part: 2*int [ t*(G(xi)-G(x)) / (x-1+t^2) ] dt ; G(x) = x/sqrt(1-x) * v(x)
        Bx = Tbasis(np.array([x]), Nv)[0]  # v(x)
        Gx = (x / np.sqrt(1 - x)) * Bx     # (Nv+1,)
        num = t[:, None] * Gq - t[:, None] * Gx[None, :]   # t*G(xi) - t*G(x)
        den = x - 1 + t**2
        I1 = 2 * (wt[:, None] * num / den[:, None]).sum(axis=0)
        I2 = Gx * np.log(x / (1 - x))
        R = Rfun(x, xi)
        I3 = 2 * (wt * R)[:, None] * (t[:, None] * Gq)
        I3 = I3.sum(axis=0)
        A[i] = I1 + I2 + I3

        # ---- f[g](x): C1 = int t g dt, C2 = -int g /2, plus double-int part
        # int_0^1 phi(xi) g(xi) dxi = 2 int_0^1 (1-t^2) v(1-t^2) phi(1-t^2) dt
        w1 = 2 * wt[:, None] * (1 - t**2)[:, None] * Bq      # for phi=1
        Ig = w1.sum(axis=0)                                   # int g
        Itg = (w1 * xi[:, None]).sum(axis=0)                  # int xi g
        C1 = Itg
        C2 = -0.5 * Ig
        # int_0^x (x-xi)^2/2 g dxi : t in [sqrt(1-x), 1]
        tx, wtx = tquad(np.sqrt(1 - x), 1.0, Mt)
        xix = 1 - tx**2
        Bqx = Tbasis(xix, Nv)
        Ix = 2 * (wtx[:, None] * (1 - tx**2)[:, None] * Bqx * ((x - xix)[:, None] ** 2 / 2))
        Ix = Ix.sum(axis=0)
        F[i] = C1 * x + C2 * x**2 + Ix

    import scipy.linalg as sla
    lam, V = sla.eig(A, -F)
    lam = lam[np.isfinite(lam)]
    sel = (np.abs(lam.imag) < 1e-8 * np.maximum(np.abs(lam.real), 1e-3)) & (lam.real > 1e-6)
    lams = np.sort(lam[sel].real)
    if verbose:
        print(f"Nv={Nv}: lambda = {lams[:nkeep]}")
    return lams, V, xc


def mode_shape(lam1, V, xc, Mt=800):
    """Reconstruct f-mode for eigenvalue closest to lam1; return f(x)/f(1) evaluator."""
    j = np.argmin(np.abs(lam1 - lam1))
    return None


if __name__ == "__main__":
    for Nv in [24, 32, 48, 64, 80]:
        eigen_g(Nv=Nv)
