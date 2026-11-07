import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.linalg import eig, null_space

gl_x, gl_w = leggauss(800)
XI = 0.5*(gl_x+1.0)
WQ = 0.5*gl_w


def Dfun(x):
    return -np.log(abs(1-x)) + np.log(1+x) + 6*x/(1+x) - 2*x**2/(1+x)**2 - 4.0


def Sreg(x, xi):
    return -1.0/(x+xi) - 2.0*xi*(x-xi)/(x+xi)**3


def pint_S_g(g, x, gp=None):
    # P \int_0^1 S(x,xi) g(xi) dxi  for a general callable g (array-capable)
    gx = g(np.array([x]))[0]
    gxi = g(XI)
    diff = x - XI
    close = np.abs(diff) < 1e-7
    safe = ~close
    term1 = np.zeros_like(diff)
    term1[safe] = (gxi[safe]-gx)/diff[safe]
    if np.any(close) and gp is not None:
        term1[close] = -gp(x)
    elif np.any(close):
        # numerical derivative fallback
        hh = 1e-6
        term1[close] = -(g(np.array([x+hh]))[0]-g(np.array([x-hh]))[0])/(2*hh)
    bracket = term1 + Sreg(x, XI)*(gxi-gx)
    bracket[close] = term1[close]   # Sreg*(g-g)=0 there anyway
    return np.sum(WQ*bracket) + gx*Dfun(x)


# validation against exact polynomial result for g=xi^1, xi^2
from numpy.polynomial import polynomial as npoly
gl_x2, gl_w2 = leggauss(600)
XI2 = 0.5*(gl_x2+1.0); WQ2 = 0.5*gl_w2
def S_integral_poly(cg, x):
    gx = npoly.polyval(x, cg)
    n = len(cg); qc = np.zeros(n-1)
    for j in range(1, n):
        for m in range(j):
            qc[m] += cg[j]*x**(j-1-m)
    intq = sum(qc[m]/(m+1) for m in range(len(qc))) if len(qc) else 0.0
    gxi = npoly.polyval(XI2, cg)
    Sregv = -1.0/(x+XI2) - 2.0*XI2*(x-XI2)/(x+XI2)**3
    return -intq + gx*np.log(abs(x/(1-x))) + np.sum(WQ2*gxi*Sregv)

print("validation pint_S_g vs exact, g=xi^2:")
for x in [0.2, 0.5, 0.8]:
    a = pint_S_g(lambda t: t**2, x)
    b = S_integral_poly([0, 0, 1.0], x)
    print("  x=%.2f  general=%.9f  exact=%.9f  diff=%.2e" % (x, a, b, a-b))


XC = XI   # evaluate operator/basis at the GL nodes; fallback handles coincidence


def build(K, Kmax=None):
    ks = list(range(1, K+1))
    def sigk(k, x):
        return x**(k/2.0)
    def sigk3coef(k):
        p = k/2.0
        return p*(p-1)*(p-2)
    def sigk3(k, x):
        return sigk3coef(k)*x**(k/2.0-3.0)
    def sigk3p(k, x):
        p = k/2.0
        return sigk3coef(k)*(p-3.0)*x**(k/2.0-4.0)
    M = np.zeros((len(XC), K))
    B = np.zeros((len(XC), K))
    for i, k in enumerate(ks):
        B[:, i] = sigk(k, XC)
    for i, k in enumerate(ks):
        g = (lambda kk: (lambda t: sigk3(kk, t)))(k)
        gp = (lambda kk: (lambda t: sigk3p(kk, t)))(k)
        for j, x in enumerate(XC):
            M[j, i] = pint_S_g(g, x, gp)
    # BCs as linear constraints on coefficients a: pick basis of nullspace
    # sigma(0)=0 -> a_k=0 for k<=0; here all k>=1 so any combo has sigma(0)=0 automatically
    # sigma'(1)=0 : sum a_k (k/2) =0 ; sigma''(1)=0: sum a_k (k/2)(k/2-1)=0
    C = np.zeros((2, K))
    for i, k in enumerate(ks):
        p = k/2.0
        C[0, i] = p
        C[1, i] = p*(p-1)
    Nsp = null_space(C)     # K x (K-2)
    return M, B, Nsp


for K in [8, 12, 16, 20, 24, 28]:
    M, B, Nsp = build(K)
    Bn = B @ Nsp
    Mn = M @ Nsp
    G = (Bn*WQ[:, None]).T @ Mn
    H = (Bn*WQ[:, None]).T @ Bn
    ev = eig(G, -H, right=False)
    ev = np.real_if_close(ev)
    real = np.sort(np.array([v.real for v in ev if abs(v.imag) < 1e-6*max(1, abs(v.real))]))
    print("K=%2d  first 6 real eigs:" % K, np.round(real[:6], 5))