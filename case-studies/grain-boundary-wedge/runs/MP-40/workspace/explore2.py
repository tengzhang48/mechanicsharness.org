import numpy as np
from numpy.polynomial.legendre import leggauss
from numpy.polynomial import polynomial as npoly

gl_x, gl_w = leggauss(500)
XI = 0.5*(gl_x+1.0)
W  = 0.5*gl_w


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
    return S1 + np.sum(W*gxi*Sreg)


def deriv_poly(c, k):
    c = np.array(c, float)
    for _ in range(k):
        c = np.polynomial.polynomial.polyder(c)
    return c


def basis_poly(k):
    c = np.array([1.0, -3.0, 3.0, -1.0])
    return np.concatenate([np.zeros(k+1), c])


def N_of(c, xs):
    c3 = deriv_poly(c, 3)
    return np.array([S_integral_poly(c3, xx) for xx in xs])


for name, c in [("x(1-x)^3", basis_poly(0)), ("x^2(1-x)^3", basis_poly(1)),
                ("x^3(1-x)^3", basis_poly(2))]:
    xsf = np.array([0.02, 0.25, 0.5, 0.75, 0.98])
    val = N_of(c, xsf)
    fval = npoly.polyval(xsf, c)
    print(name)
    print("   x        ", np.round(xsf, 4))
    print("   f        ", np.round(fval, 6))
    print("   N[f]     ", np.round(val, 6))
    print("   N/f      ", np.round(val/fval, 6))