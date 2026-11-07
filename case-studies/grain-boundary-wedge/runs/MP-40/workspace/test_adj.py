import numpy as np
from numpy.polynomial.legendre import leggauss

XI, WQ = np.polynomial.legendre.leggauss(2000)
XI = 0.5*(XI+1.0); WQ = 0.5*WQ
gl2, w2 = leggauss(2000)
XO = 0.5*(gl2+1.0); WO = 0.5*w2


def Dfun(x):
    return -np.log(abs(1-x)) + np.log(1+x) + 6*x/(1+x) - 2*x**2/(1+x)**2 - 4.0


def pint(g, gp, xv):
    gx = g(xv); gxi = g(XI)
    diff = xv[:, None] - XI[None, :]
    close = np.abs(diff) < 1e-12
    num = gxi[None, :] - gx[:, None]
    ratio = np.where(close, -gp(xv)[:, None], num/np.where(close, 1.0, diff))
    integ = 8.0*xv[:, None]*XI[None, :]**2*ratio/(xv[:, None]+XI[None, :])**3
    integ = np.where(close, 0.0, integ)
    return integ @ WQ + gx*Dfun(xv)


def bump(x):
    return np.exp(-1.0/np.maximum(x*(1-x), 1e-30))


def make(f, h=1e-4):
    # return f and its 1st/3rd derivatives callables
    f3 = lambda t: (f(t+2*h)-2*f(t+h)+2*f(t-h)-f(t-2*h))/(2*h**3)
    g1 = lambda t: (f(t+h)-f(t-h))/(2*h)
    return f, g1, f3


f = bump
ff = lambda t: bump(np.asarray(t, float))
g, g1, g3 = make(ff)

Nf = pint(g3, lambda t: (g3(t+1e-4)-g3(t-1e-4))/2e-4, XO)

# another test function
def ff2(t):
    t = np.asarray(t, float)
    return np.exp(-1.0/np.maximum(t*(1-t), 1e-30))*(t-0.5)

g2_1 = lambda t: (ff2(t+1e-4)-ff2(t-1e-4))/2e-4
def g2_3(t):
    h = 1e-3
    return (ff2(t+2*h)-2*ff2(t+h)+2*ff2(t-h)-ff2(t-2*h))/(2*h**3)
Nf2 = pint(g2_3, lambda t: (g2_3(t+1e-4)-g2_3(t-1e-4))/2e-4, XO)

A = np.sum(WO*ff2(XO)*Nf)     # <psi, N phi>
B = np.sum(WO*ff(XO)*Nf2)     # <phi, N psi>
print("<psi,Nphi> =", A)
print("<phi,Npsi> =", B)
print("difference    =", A-B, " relative:", abs(A-B)/max(abs(A), abs(B)))

# also check the singular (Hilbert) part only:  P int g'''(xi)/(x-xi) dxi
def hilb(g, gp, xv):
    gx = g(xv); gxi = g(XI)
    diff = xv[:, None]-XI[None, :]
    close = np.abs(diff) < 1e-12
    num = gxi[None, :]-gx[:, None]
    ratio = np.where(close, -gp(xv)[:, None], num/np.where(close, 1.0, diff))
    ratio = np.where(close, -gp(xv)[:, None], ratio)
    return (ratio @ WQ) + gx*np.log(np.abs(xv/(1-xv)))

NfH = hilb(g3, lambda t: (g3(t+1e-4)-g3(t-1e-4))/2e-4, XO)
Nf2H = hilb(g2_3, lambda t: (g2_3(t+1e-4)-g2_3(t-1e-4))/2e-4, XO)
AH = np.sum(WO*ff2(XO)*NfH)
BH = np.sum(WO*ff(XO)*Nf2H)
print("Hilbert part: <psi,Nphi>=", AH, " <phi,Npsi>=", BH, " diff=", AH-BH)