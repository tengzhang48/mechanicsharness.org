import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.linalg import eig, null_space

NMU = 1400
gl_x, gl_w = leggauss(NMU)
XI = 0.5*(gl_x+1.0)
WQ = 0.5*gl_w


def Dfun(x):
    return -np.log(abs(1-x)) + np.log(1+x) + 6*x/(1+x) - 2*x**2/(1+x)**2 - 4.0


def pint_S_g(g, x):
    gx = g(np.array([x]))[0]
    gxi = g(XI)
    diff = x - XI
    close = np.abs(diff) < 1e-9
    ratio = np.empty_like(diff)
    safe = ~close
    ratio[safe] = (gxi[safe]-gx)/diff[safe]
    if np.any(close):
        hh = 1e-5
        ratio[close] = -(g(np.array([x+hh]))[0]-g(np.array([x-hh]))[0])/(2*hh)
    integ = 8.0*x*XI**2*ratio/(x+XI)**3
    integ[close] = 0.0
    return np.sum(WQ*integ) + gx*Dfun(x)


def basis_ops(ps):
    K = len(ps)
    Bq = np.zeros((len(XI), K)); Nq = np.zeros((len(XI), K))
    for i, p in enumerate(ps):
        c3 = p*(p-1)*(p-2)
        g = (lambda pp, cc: (lambda t: cc*t**(pp-3.0)))(p, c3)
        Bq[:, i] = XI**p
        Nq[:, i] = np.array([pint_S_g(g, x) for x in XI])
    G = (Bq*WQ[:, None]).T @ Nq
    H = (Bq*WQ[:, None]).T @ Bq
    C = np.zeros((2, K))
    for i, p in enumerate(ps):
        C[0, i] = p
        C[1, i] = p*(p-1)
    return G, H, C


def run(ps):
    G, H, C = basis_ops(ps)
    Z = null_space(C)
    G2 = Z.T @ G @ Z; H2 = Z.T @ H @ Z
    ev, V = eig(G2, -H2)
    re = np.array([v.real for v in ev if abs(v.imag) < 1e-6*max(1, abs(v.real)) and v.real > 0])
    idx = np.argsort(re)
    re = re[idx]
    # eigenvector of smallest
    order = [k for k in np.argsort(ev.real) if abs(ev[k].imag) < 1e-6*max(1, abs(ev[k].real)) and ev[k].real > 0]
    return re, Z, V, order


print("Integer basis p=1..N")
for N in [8, 10, 12, 14, 16, 18]:
    ps = list(np.arange(1.0, 1.0+N))
    re, Z, V, order = run(ps)
    print("N=%2d  lam1=%.5f lam2=%.4f lam3=%.3f lam4=%.2f" %
          (N, re[0], re[1], re[2], re[3]))

print()
print("pi^4/8 =", np.pi**4/8, "  pi^4 =", np.pi**4, "  (3*pi/4)^3*pi=", np.pi*(0.75*np.pi)**3)
print("(pi/2)^3*pi =", np.pi*(0.5*np.pi)**3)
print("(3pi/2)^3*pi =", np.pi*(1.5*np.pi)**3, " (5pi/2)^3*pi=", np.pi*(2.5*np.pi)**3)