import numpy as np
from numpy.polynomial import Polynomial as Poly
from numpy.polynomial.legendre import leggauss
from explore import apply_R
from scipy.linalg import eig

xg, wg = leggauss(1500); xg=0.5*(xg+1); wg=0.5*wg

def phi0_poly():
    Mcon = np.zeros((6,6)); rhs=np.array([0,1,0,0,0,0.0])
    for k in range(6):
        p = Poly([0]*k+[1.0])
        Mcon[0,k]=p(0); Mcon[1,k]=p(1); Mcon[2,k]=p.deriv(1)(1)
        Mcon[3,k]=p.deriv(2)(1); Mcon[4,k]=p.deriv(3)(1); Mcon[5,k]=p.deriv(3)(0)
    return Poly(np.linalg.solve(Mcon,rhs))
phi0=phi0_poly()

def make_basis(N, pmin=3):
    psi=[phi0]
    for k in range(N):
        psi.append(Poly([0]*(pmin+k)+[1.0]) * Poly([1.0,-4,6,-4,1.0]))
    return psi

def galerkin(N,pmin=3):
    psi=make_basis(N,pmin); M=len(psi)
    # R[psi_j] evaluated on the Gauss grid
    Rp=np.zeros((M,len(xg)))
    for j in range(M):
        p3=psi[j].deriv(3)
        for k,Z in enumerate(xg):
            Rp[j,k]=apply_R(p3,Z,xg,wg)
    # Hmm apply_R with Z on the grid: use analytic. But apply_R signature may vectorize? do loop
    Mm=np.zeros((M,M)); Km=np.zeros((M,M))
    Pv=np.array([[psi[i](Z) for Z in xg] for i in range(M)])
    for i in range(M):
        for j in range(M):
            Mm[i,j]=np.sum(wg*Pv[i]*Pv[j])
            Km[i,j]=np.sum(wg*Pv[i]*Rp[j])
    lam,vec=eig(Km,Mm)
    cand=sorted([(l.real,k) for k,l in enumerate(lam) if abs(l.imag)<1e-6*max(1,abs(l.real)) and l.real>1e-6])
    return [c[0] for c in cand[:4]], psi, vec, cand

if __name__ == "__main__":
    for pmin in [3,4]:
        for N in [6,10,14,18]:
            vals,psi,vec,cand=galerkin(N,pmin)
            print("pmin=%d N=%2d lams:"%(pmin,N), np.round(vals,5))