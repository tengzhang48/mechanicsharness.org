import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.linalg import eig

SRC = open('<round-runs>/runs/MP-40/explore7.py').read().split('for K in')[0]


def make(NMU):
    ns = {}
    exec(SRC.replace("NMU = 1600", "NMU = %d" % NMU), ns)
    return ns


def opmatrix(ns, K):
    XI, WQ = ns['XI'], ns['WQ']
    pint_S_g = ns['pint_S_g']
    cheb_basis_poly = ns['cheb_basis_poly']
    dcoef = ns['dcoef']
    Bq = np.zeros((len(XI), K)); Nq = np.zeros((len(XI), K))
    for k in range(K):
        c = cheb_basis_poly(k); c3 = dcoef(c, 3)
        Bq[:, k] = np.polynomial.polynomial.polyval(XI, c)
        g = (lambda cc: (lambda t: np.polynomial.polynomial.polyval(t, cc)))(c3)
        Nq[:, k] = np.array([pint_S_g(g, x) for x in XI])
    G = (Bq*WQ[:, None]).T @ Nq
    H = (Bq*WQ[:, None]).T @ Bq
    return G, H


for NMU in [300, 600, 1200, 2400, 4800]:
    ns = make(NMU)
    G, H = opmatrix(ns, 14)
    sym = np.linalg.norm(G-G.T)/np.linalg.norm(G)
    ev = np.sort(np.array([v.real for v in eig(G, -H, right=False)
                           if abs(v.imag) < 1e-6*max(1, abs(v.real)) and v.real > 0]))
    print("NMU=%5d  sym=%.3e  first eigs:" % (NMU, sym), np.round(ev[:4], 5))