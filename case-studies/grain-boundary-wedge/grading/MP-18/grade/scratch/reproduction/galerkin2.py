import numpy as np
import scipy.linalg as sla
from numpy.polynomial import chebyshev as C2
from numpy.polynomial.legendre import leggauss
from spectral_galerkin import cheb_basis_with_BCs


def Jm_all(xq, Mmax, ny=4000):
    gy, wy = leggauss(ny)
    yq = 0.5 * (gy + 1.0)
    wy = 0.5 * wy
    X = xq[:, None]
    Y = yq[None, :]
    S = 1.0 / (Y - X) + 1.0 / (X + Y) + 2.0 * Y * (X - Y) / (X + Y) ** 3
    S = np.where(np.abs(Y - X) < 1e-13, 0.0, S)
    Jmat = np.zeros((Mmax + 1, len(xq)))
    for m in range(Mmax + 1):
        Jmat[m] = (wy * Y ** m * S).sum(axis=1)
    return Jmat


def run(Nc, nx=600, ny=2000):
    ns = cheb_basis_with_BCs(Nc)
    nb = ns.shape[0]
    dx_, wx = leggauss(nx)
    xq = 0.5 * (dx_ + 1.0)
    wx = 0.5 * wx
    Jmat = Jm_all(xq, Nc, ny=ny)         # (Nc+1, nx)
    PsiX = np.array([C2.chebval(2 * xq - 1, ns[j]) for j in range(nb)])   # (nb,nx)
    # Build monomial-coefficient matrix of psi_i''' in x: shape (nb, deg)
    C3mono = np.zeros((nb, Nc + 1))
    for i in range(nb):
        c3 = C2.chebder(C2.chebder(C2.chebder(ns[i])))
        mono = C2.cheb2poly(c3)
        acc = np.zeros(1)
        for k in range(len(mono) - 1, -1, -1):
            acc = np.polymul(acc, [2.0, -1.0])
            acc = np.polyadd(acc, [mono[k]])
        C3mono[i, :len(acc)] = acc
    C3mono = C3mono[:, :Nc + 1]
    # inner_i(x) = sum_m C3mono[i,m] Jmat[m,x]  -> (nb, nx)
    Inner = C3mono @ Jmat[:Nc + 1, :]
    # G[i,j] = -(1/pi) sum_x wx * PsiX[j,x] * Inner[i,x]
    G = -(1.0 / np.pi) * (Inner @ (PsiX * wx).T)     # (nb, nb)
    Mass = (PsiX * wx) @ PsiX.T
    return ns, G, Mass


if __name__ == "__main__":
    for Nc in [10, 14, 18, 22, 26]:
        ns, G, Mass = run(Nc, nx=600, ny=2000)
        ev = sla.eig(G, Mass, right=False)
        ev = ev.real[np.argsort(-ev.real)]
        print("Nc=%2d nb=%2d top eigs: %s" % (Nc, ns.shape[0], np.round(ev[:5], 5)))