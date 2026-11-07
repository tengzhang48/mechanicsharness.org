import numpy as np
from numpy.polynomial import Polynomial as Poly
from numpy.polynomial.legendre import leggauss

# Kernel (given): S(z,xi) = 1/(z-xi) - 1/(z+xi) - 2 xi (z-xi)/(z+xi)^3
# Rewrite: S = 1/(z-xi) + Rreg,  Rreg = (xi^2 - 4 z xi - z^2)/(z+xi)^3
def Rreg(Z, Xi):
    return (Xi**2 - 4*Z*Xi - Z**2)/(Z+Xi)**3

def Stilde(Z, Xi):
    return 1.0/(Z-Xi) - 1.0/(Z+Xi) - 2*Xi*(Z-Xi)/(Z+Xi)**3

# PV integral of polynomial P over [0,1]:  PV int_0^1 P(xi)/(Z-xi) dxi
# = -[ int_0^1 Q dxi + P(Z) ln((1-Z)/Z) ],  Q=(P(xi)-P(Z))/(xi-Z)
def pv_poly_over_Zminus(P, Z):
    c = P.coef                      # ascending
    n = len(c)-1
    integ = 0.0
    if n >= 1:
        # synthetic division by (xi - Z): Q coeffs q_0..q_{n-1}, q_{n-1}=c_n, q_{k-1}=c_k+Z q_k
        q = np.zeros(n)
        q[n-1] = c[n]
        for k in range(n-1, 0, -1):
            q[k-1] = c[k] + Z*q[k]
        integ = np.sum([q[k]/(k+1) for k in range(n)])
    val = -(integ + P(Z)*np.log((1.0-Z)/Z))
    return val

def apply_R(P, Z, xg, wg):
    """R[f](Z) = -int_0^1 S(Z,xi) f'''(xi) dxi,  with P=psi'''."""
    sing = pv_poly_over_Zminus(P, Z)
    reg = np.sum(wg * Rreg(Z, xg) * P(xg))
    return -(sing + reg)

def make_basis(N):
    """psi_0 = 1+(Z-1)^3 ; psi_k = Z^k (Z-1)^3, k=1..N"""
    base = Poly([-1.0, 3.0, -3.0, 1.0])   # (Z-1)^3
    psi = []
    psi.append(Poly([1.0,0,0,0]) + base) # 1+(Z-1)^3
    for k in range(1, N+1):
        psi.append(base * Poly([0.0]*k + [1.0]))
    return psi

def eigen(N, verbose=True):
    psi = make_basis(N)
    P3 = [p.deriv(3) for p in psi]
    xg, wg = leggauss(400)
    xg = 0.5*(xg+1.0); wg = 0.5*wg
    # collocation points (Chebyshev interior of (0,1)); square => #points = #basis
    M = len(psi)
    th = (np.arange(1,M+1))*np.pi/(M+1)
    Zc = 0.5*(1-np.cos(th))
    Rmat = np.zeros((len(Zc), len(psi)))
    Mmat = np.zeros((len(Zc), len(psi)))
    for i,Z in enumerate(Zc):
        for j in range(len(psi)):
            Rmat[i,j] = apply_R(P3[j], Z, xg, wg)
            Mmat[i,j] = psi[j](Z)
    from scipy.linalg import eig
    lam, vec = eig(Rmat, Mmat)
    idx = np.argsort(lam.real)
    lam = lam[idx]; vec = vec[:,idx]
    return lam, vec, Zc, psi

if __name__ == "__main__":
    import sys
    N = int(sys.argv[1]) if len(sys.argv)>1 else 20
    lam, vec, Zc, psi = eigen(N)
    order = np.argsort([l.real for l in lam])
    # pick real positive smallest
    cand = [(l.real, k) for k,l in enumerate(lam) if abs(l.imag)<1e-6*max(1,abs(l.real)) and l.real>1e-6]
    cand.sort()
    for (lr,k) in cand[:4]:
        c = vec[:,k]
        c0 = c[0]
        # mode shape
        Zs = np.array([0.0,0.25,0.5,0.75,1.0])
        vals = np.array([sum(c[j]*psi[j](Z) for j in range(len(psi))) for Z in Zs])
        fn = vals/ (sum(c[j]*psi[j](1.0) for j in range(len(psi))))
        print("lam=%.6f" % lr, "f(0,.25,.5,.75,1)=", np.round(fn,6), "imag=%.1e"%lam[k].imag)