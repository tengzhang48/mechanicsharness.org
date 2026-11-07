import numpy as np
from numpy.polynomial.legendre import leggauss

NMU = 2000
gl_x, gl_w = leggauss(NMU)
XII = 0.5*(gl_x+1.0); WQf = 0.5*gl_w


def Df(x):
    return -np.log(abs(1-x)) + np.log(1+x) + 6*x/(1+x) - 2*x**2/(1+x)**2 - 4.0


def applyN(f, fp, fp3, xv):
    # N f(x) = P int_0^1 S(x,xi) f3(xi) dxi, via subtraction with g=f3
    g = fp3(XII); gx = fp3(xv)
    diff = xv[:, None] - XII[None, :]
    close = np.abs(diff) < 1e-9
    ratio = np.where(close, -fp3p(xv)[:, None], (g[None, :]-gx[:, None])/np.where(close, 1.0, diff))
    integ = 8.0*xv[:, None]*XII[None, :]**2*ratio/(xv[:, None]+XII[None, :])**3
    return integ @ WQf + gx*Df(xv)


# ---------- build explore13-style eigenpairs ----------
import importlib.util
spec = importlib.util.spec_from_file_location("e13", "<round-runs>/runs/MP-40/explore13.py")
# instead replicate:
from numpy.polynomial import chebyshev as C


def chebT_poly(k):
    return C.cheb2poly([0]*k+[1.])


def to01(c):
    out = np.array([0.0])
    for i, ci in enumerate(c):
        term = np.array([ci])
        for _ in range(i):
            term = np.convolve(term, np.array([-1.0, 2.0]))
        if len(term) > len(out):
            out = np.concatenate([out, np.zeros(len(term)-len(out))])
        out[:len(term)] += term
    return out


def deriv(c, k):
    for _ in range(k):
        c = np.polynomial.polynomial.polyder(c)
    return c


from scipy.linalg import null_space

NMU2 = 1200
gl2, w2 = leggauss(NMU2)
XI = 0.5*(gl2+1.0); WQ = 0.5*w2


def pint_S_all(g, gp, xv):
    gx = g(xv); gxi = g(XI)
    diff = xv[:, None] - XI[None, :]
    close = np.abs(diff) < 1e-9
    num = gxi[None, :] - gx[:, None]
    ratio = np.where(close, -gp(xv)[:, None], num/np.where(close, 1.0, diff))
    integ = 8.0*xv[:, None]*XI[None, :]**2*ratio/(xv[:, None]+XI[None, :])**3
    return integ @ WQ + gx*Dfun_test(xv)


def Dfun_test(x):
    return -np.log(abs(1-x)) + np.log(1+x) + 6*x/(1+x) - 2*x**2/(1+x)**2 - 4.0


K = 20
polys = [to01(chebT_poly(k)) for k in range(K+1)]
Cm = np.zeros((3, K+1))
for k, p in enumerate(polys):
    Cm[0, k] = np.polynomial.polynomial.polyval(0.0, p)
    Cm[1, k] = np.polynomial.polynomial.polyval(1.0, deriv(p, 1))
    Cm[2, k] = np.polynomial.polynomial.polyval(1.0, deriv(p, 2))
Z = null_space(Cm)
Bmat = np.zeros((len(XI), K+1)); Nmat = np.zeros((len(XI), K+1))
for k, p in enumerate(polys):
    Bmat[:, k] = np.polynomial.polynomial.polyval(XI, p)
    c3 = deriv(p, 3)
    g = (lambda cc: (lambda t: np.polynomial.polynomial.polyval(t, cc)))(c3)
    gp = (lambda cc: (lambda t: np.polynomial.polynomial.polyval(t, deriv(cc, 1))))(c3)
    Nmat[:, k] = pint_S_all(g, gp, XI)
Bn = Bmat @ Z; Nn = Nmat @ Z
n = Bn.shape[1]; R = np.eye(n)
for j in range(n):
    v = Bn[:, j].copy(); r = R[:, j].copy()
    for i in range(j):
        c = np.sum(WQ*Bn[:, i]*v); v -= c*Bn[:, i]; r -= c*R[:, i]
    R[:, j] = r/np.sqrt(np.sum(WQ*v*v))
Q = Bn @ R; NQ = Nn @ R
G = (Q*WQ[:, None]).T @ NQ; H = (Q*WQ[:, None]).T @ Q
ev, V = np.linalg.eig(np.linalg.solve(H, G))
cands = [(-ev[k].real, k) for k in range(len(ev)) if abs(ev[k].imag) < 1e-6*max(1, abs(ev[k])) and -ev[k].real > 0]
cands.sort()
print("chebyshev-basis eigenvalues:", [round(c[0], 4) for c in cands[:4]])

# reconstruct sigma as polynomial in monomial basis on [0,1]
for lam, k in cands[:2]:
    a = np.real(V[:, k])
    coef = Z @ (R @ a)              # coefficients in [1,x,x^2,...]
    f = lambda t: np.polynomial.polynomial.polyval(t, coef)
    f1 = lambda t: np.polynomial.polynomial.polyval(t, deriv(coef, 1))
    f2 = lambda t: np.polynomial.polynomial.polyval(t, deriv(coef, 2))
    f3 = lambda t: np.polynomial.polynomial.polyval(t, deriv(coef, 3))
    global fp3p
    fp3p = f2  # derivative of g=f3 is f4; but for subtraction we need g'=f4
    f4 = lambda t: np.polynomial.polynomial.polyval(t, deriv(coef, 4))
    fp3p = f4
    xv = np.array([0.1, 0.25, 0.5, 0.75, 0.9])
    Ns = applyN(f, f1, f3, xv)
    sig = f(xv)
    print("lam=%.5f  Nsigma/-sigma ratio at x=%s:" % (lam, list(xv)))
    print("   ", np.round(Ns/-sig, 4), " resid(rel)=", np.max(np.abs(Ns+lam*sig))/np.max(np.abs(lam*sig)))