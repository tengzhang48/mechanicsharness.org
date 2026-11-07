import numpy as np
from numpy.polynomial.legendre import leggauss
from numpy.polynomial import polynomial as npoly
from scipy.linalg import eig

gl_x, gl_w = leggauss(600)
XI = 0.5*(gl_x+1.0)
WQ = 0.5*gl_w


def S_integral_poly(cg, x):
    # P int_0^1 S(x,xi) g(xi) dxi for polynomial g with coeffs cg (increasing)
    gx = npoly.polyval(x, cg)
    n = len(cg)
    qc = np.zeros(n-1)
    for j in range(1, n):
        for m in range(j):
            qc[m] += cg[j]*x**(j-1-m)
    intq = sum(qc[m]/(m+1) for m in range(len(qc))) if len(qc) else 0.0
    S1 = -intq + gx*np.log(abs(x/(1.0-x)))
    gxi = npoly.polyval(XI, cg)
    Sreg = -1.0/(x+XI) - 2.0*XI*(x-XI)/(x+XI)**3
    return S1 + np.sum(WQ*gxi*Sreg)


def Phi1(x):
    # Phi(x,1) = -ln|1-x| + ln(1+x) + 6x/(1+x) - 2x^2/(1+x)^2 ; int_0^1 S dxi = Phi(x,1)-4
    return -np.log(abs(1-x)) + np.log(1+x) + 6*x/(1+x) - 2*x**2/(1+x)**2


# check g=1 (cg=[1]); int S*1 = Phi(x,1)-4
print("quadrature check  int S dxi :")
for x in [0.1, 0.3, 0.5, 0.7, 0.9]:
    print("  x=%.2f  quad=%.10f  analytic=%.10f" % (x, S_integral_poly([1.0], x), Phi1(x)-4.0))


def deriv_poly(c, k):
    c = np.array(c, float)
    for _ in range(k):
        c = np.polynomial.polynomial.polyder(c)
    return c


def cheb_basis_poly(k):
    # phi_k = x (1-x)^3 * T_k(2x-1),  T_k via numpy
    Tk = np.polynomial.chebyshev.cheb2poly([0]*k + [1])
    # T_k(y), y=2x-1 -> poly in x
    cx = np.polynomial.polynomial.polyval  # not needed
    # substitute y=2x-1
    p = np.array([0.0])
    for i, c in enumerate(Tk):
        # c*(2x-1)^i
        term = np.array([c])
        for _ in range(i):
            term = np.convolve(term, np.array([-1.0, 2.0]))
        if len(term) > len(p):
            p = np.concatenate([p, np.zeros(len(term)-len(p))])
        p[:len(term)] += term
    xfac = np.array([0.0, 1.0])          # x
    one_m_x3 = np.array([1.0, -3.0, 3.0, -1.0])
    p = np.convolve(p, xfac)
    p = np.convolve(p, one_m_x3)
    return p


# GL nodes on (0,1) for projection
def opmatrix(K):
    # basis polys
    B = [cheb_basis_poly(k) for k in range(K)]
    # evaluate N[phi_k] at GL nodes
    Nmat = np.zeros((len(XI), K))
    for k in range(K):
        c3 = deriv_poly(B[k], 3)
        Nmat[:, k] = np.array([S_integral_poly(c3, x) for x in XI])
    Bmat = np.zeros((len(XI), K))
    for k in range(K):
        Bmat[:, k] = npoly.polyval(XI, B[k])
    G = (Bmat*WQ[:, None]).T @ Nmat
    H = (Bmat*WQ[:, None]).T @ Bmat
    return G, H, B, Bmat


for K in [8, 12, 16, 20]:
    G, H, B, Bmat = opmatrix(K)
    ev = eig(G, -H, right=False)
    ev = np.real_if_close(ev)
    pos = np.sort(np.array([v.real for v in ev if abs(v.imag) < 1e-4*max(1, abs(v.real)) and v.real > 0]))
    print("K=%d  first 6 positive rates:" % K, np.round(pos[:6], 5))