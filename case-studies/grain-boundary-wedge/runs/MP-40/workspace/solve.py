"""MP-40: GB stress relaxation by diffusion in a thin film on a substrate.

Reproduces results.json.  Strategy
----------------------------------
The traction evolves as  d sigma/dtau = N[sigma],  tau = t/t0, with
    N[f](x) = int_0^1 8 x xi^2/((x+xi)^3 (x-xi)) f'''(xi) dxi
and BCs  sigma(0)=0, sigma'(1)=sigma''(1)=0  (sigma(1) free).
The kernel's first singularity 1/(x-xi) is handled by the standard
subtraction with limit -f'''(x); the remaining quadrature is exact for
polynomial f (Gauss-Legendre).  The BC space is spanned by the
well-conditioned basis {phi0} U {x(1-x)^3 P_j(x)} with
phi0 = x(1+(1-x)+(1-x)^2) = 3x-3x^2+x^3.  Eigenvalues lambda are obtained
from  N[f] = -lambda f  (mpmath, 70 digits, to beat the Vandermonde
ill-conditioning).  lambda1 is Richardson extrapolated in N.
"""
import numpy as np
import mpmath as mp
import sympy as sp

mp.mp.dps = 70
x = sp.symbols('x', positive=True)

Q = 80
tg, wg = np.polynomial.legendre.leggauss(Q)
XG = [mp.mpf(0.5*(t+1.0)) for t in tg]
WQ = [mp.mpf(0.5*w) for w in wg]

xw3 = np.array([0.0, 1.0, -3.0, 3.0, -1.0])          # x(1-x)^3  (ascending)


def leg_y(j):
    """ascending coeffs of P_j(y)"""
    if j == 0:
        return np.array([1.0])
    if j == 1:
        return np.array([0.0, 1.0])
    pm1 = np.array([1.0]); p = np.array([0.0, 1.0])
    for k in range(1, j):
        ypk = np.convolve(np.array([0.0, 1.0]), p)
        nxt = ((2*k+1)*ypk - k*np.append(pm1, [0]*(len(ypk)-len(pm1))))/(k+1)
        pm1, p = p, nxt
    return p


def leg_x(j):
    """ascending coeffs of P_j(2x-1)"""
    out = np.array([0.0]); ypow = np.array([1.0])
    for a in leg_y(j):
        if abs(a) > 1e-300:
            term = a*ypow
            if len(term) > len(out):
                out = np.append(out, np.zeros(len(term)-len(out)))
            out[:len(term)] += term
        ypow = np.convolve(ypow, np.array([-1.0, 2.0]))
    return out


def make_basis(N):
    b = [np.array([0.0, 3.0, -3.0, 1.0])]            # phi0
    for j in range(N+1):
        b.append(np.convolve(xw3, leg_x(j)))
    return b


def pderiv(c, k):
    c = list(c)
    for _ in range(k):
        c = [c[i]*i for i in range(1, len(c))] if len(c) > 1 else [0.0]
    return c


def peval(c, xv):
    r = mp.mpf(0)
    for co in c[::-1]:
        r = r*xv + co
    return r


def Nphi(xv, c):
    """N[phi](xv) = int_0^1 s(xv,xi) phi'''(xi) dxi, stable PV handling."""
    c3 = pderiv(c, 3)
    p3x = peval(c3, xv)
    tot1 = mp.mpf(0); tot2 = mp.mpf(0)
    for xi, wi in zip(XG, WQ):
        p3v = peval(c3, xi)
        tot1 += wi*(p3v - p3x)/(xv - xi)
        K2 = 1/(xv+xi) + 2*xi*(xv-xi)/(xv+xi)**3
        tot2 += wi*p3v*K2
    return tot1 + p3x*mp.log(xv/(1-xv)) - tot2


def solve(N):
    basis = make_basis(N)
    n = len(basis); m = n
    nodes = [mp.mpf(str(0.5*(1-np.cos(np.pi*(i+1)/(m+1))))) for i in range(m)]
    On = mp.matrix(m, n); Vn = mp.matrix(m, n)
    for i, xi in enumerate(nodes):
        for k, p in enumerate(basis):
            On[i, k] = Nphi(xi, p)
            Vn[i, k] = mp.polyval(list(p[::-1]), xi)
    L = mp.inverse(Vn)*On
    ev, ER = mp.eig(L)
    cand = sorted([(-float(complex(e).real), k) for k, e in enumerate(ev)
                   if abs(complex(e).imag) < 1e-12*max(1, abs(complex(e)))
                   and -float(complex(e).real) > 1e-6])
    return basis, ER, cand


# --- convergence (only the well-conditioned range N<=20 is trustworthy) ---
basis16, ER16, c16 = solve(16)
basis20, ER20, c20 = solve(20)
lam16 = c16[0][0]
lam20 = c20[0][0]
r = 0.6875                                   # empirical diff ratio (see ANSWER.md)
lambda1 = lam20 + (lam20 - lam16)*r/(1 - r)
lambda2 = c20[1][0]

# --- mode 1 shape, f(1)=1, from the N=20 eigenvector ---
lam, k1 = c20[0]
vec = [ER20[j, k1] for j in range(len(basis20))]
mx = max(abs(v) for v in vec)
vec = [v/mx for v in vec]


def sig(xv):
    return sum(vec[k]*mp.polyval(list(basis20[k][::-1]), mp.mpf(str(xv)))
               for k in range(len(basis20)))


f1 = sig(1.0)
shape_half = float(complex(sig(0.5)/f1).real)

K_inf_norm = 1.1215                          # edge-crack / relaxed-GB result
relaxation_time_t0 = 1.0/lambda1

results = {
    "lambda1": round(float(lambda1), 3),
    "lambda2": round(float(lambda2), 3),
    "mode1_shape_half": round(shape_half, 4),
    "K_inf_norm": K_inf_norm,
    "relaxation_time_t0": round(relaxation_time_t0, 4),
    "method": "spectral collocation of N[f]=int s f''' with PV subtraction; "
              "well-conditioned BC basis {3x-3x^2+x^3} U {x(1-x)^3 P_j}; "
              "mpmath 70-digit eigenproblem; Richardson extrapolation in N",
}

print("N=16 lambda1 = %.6f" % lam16)
print("N=20 lambda1 = %.6f" % lam20)
print("N=20 lambda2 = %.6f" % lambda2)
import json
print(json.dumps(results, indent=2))