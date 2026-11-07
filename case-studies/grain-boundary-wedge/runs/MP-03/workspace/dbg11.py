import numpy as np
from eig_g import tquad
from gbcore import Rfun
from numpy.polynomial.chebyshev import chebval

def solve_mode(Nv=48, Mt=800, lam_target=None):
    """Return (lams, a_coeffs of v) for pencil; eigenvector for eigenvalue nearest lam_target."""
    import scipy.linalg as sla
    def Tbasis(xi, mmax):
        return np.cos(np.outer(np.arccos(np.clip(2*xi-1,-1,1)), np.arange(mmax+1)))
    k = np.arange(1, Nv+2)
    xc = 0.5*(1-np.cos(np.pi*k/(Nv+2)))
    A = np.zeros((Nv+1,Nv+1)); F = np.zeros((Nv+1,Nv+1))
    for i,x in enumerate(xc):
        t, wt = tquad(0.0,1.0,Mt); xi = 1-t**2
        Bq = Tbasis(xi,Nv); Bx = Tbasis(np.array([x]),Nv)[0]
        Gq = xi[:,None]/t[:,None]*Bq
        Gx = (x/np.sqrt(1-x))*Bx
        den = x-1+t**2
        I1 = 2*(wt[:,None]*(t[:,None]*Gq - t[:,None]*Gx[None,:])/den[:,None]).sum(axis=0)
        I2 = Gx*np.log(x/(1-x))
        R = Rfun(x,xi)
        I3 = 2*((wt*R)[:,None]*(t[:,None]*Gq)).sum(axis=0)
        A[i] = I1+I2+I3
        w1 = 2*wt[:,None]*(1-t**2)[:,None]*Bq
        Ig = w1.sum(axis=0); Itg = (w1*xi[:,None]).sum(axis=0)
        tx, wtx = tquad(np.sqrt(1-x),1.0,Mt)
        xix = 1-tx**2; Bqx = Tbasis(xix,Nv)
        Ix = 2*(wtx[:,None]*(1-tx**2)[:,None]*Bqx*((x-xix)[:,None]**2/2)).sum(axis=0)
        F[i] = Itg*x + (-0.5*Ig)*x**2 + Ix
    lam, V = sla.eig(A, -F)
    sel = (np.abs(lam.imag)<1e-8*np.maximum(np.abs(lam.real),1e-3)) & (lam.real>1e-6)
    lams = np.sort(lam[sel].real)
    lam1 = lams[0]
    # eigenvector: column with eigenvalue closest to lam1
    j = np.argmin(np.abs(lam - lam1))
    a = np.real(V[:, j])
    # f-mode evaluator: f(x) = C1 x + C2 x^2 + int_0^x (x-xi)^2/2 g dxi
    def fmode(x):
        t, wt = tquad(0.0,1.0,Mt); xi = 1-t**2
        vq = chebval(2*xi-1, a)
        # int g = 2 int (1-t^2) v dt ; int xi g = 2 int (1-t^2) xi v dt
        Ig = 2*np.sum(wt*(1-t**2)*vq)
        Itg = 2*np.sum(wt*(1-t**2)*xi*vq)
        tx, wtx = tquad(np.sqrt(1-x),1.0,Mt); xix = 1-tx**2
        vx = chebval(2*xix-1, a)
        Ix = 2*np.sum(wtx*(1-tx**2)*vx*(x-xix)**2/2)
        return Itg*x - 0.5*Ig*x**2 + Ix
    return lams, a, fmode

lams, a, fmode = solve_mode(Nv=48)
f1 = fmode(1.0)
print("lams:", lams[:4])
print("f(1) =", f1)
for x in [0.1,0.25,0.5,0.75,0.9,0.99]:
    print(f"f({x})/f(1) = {fmode(x)/f1:+.6f}")
# residual check: A[g] + lam f on independent dense grid
import numpy as np
from gbcore import Rfun
from numpy.polynomial.chebyshev import chebval
Mt=800
def resid(x, lam):
    t, wt = tquad(0.0,1.0,Mt); xi = 1-t**2
    vq = chebval(2*xi-1, a)
    Gq = xi/t*vq
    vx = chebval(2*x-1, a)
    Gx = x/np.sqrt(1-x)*vx
    den = x-1+t**2
    I1 = 2*np.sum(wt*(t*Gq - t*Gx)/den)
    I2 = Gx*np.log(x/(1-x))
    I3 = 2*np.sum(wt*Rfun(x,xi)*t*Gq)
    return (I1+I2+I3) + lam*fmode(x)
for x in [0.02,0.1,0.3,0.5,0.7,0.9,0.98]:
    print(f"x={x}: |resid| = {abs(resid(x, lams[0])):.2e}   (|lam f| ~ {abs(lams[0]*fmode(x)):.3f})")
