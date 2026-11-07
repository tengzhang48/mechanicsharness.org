import numpy as np
from numpy.polynomial.legendre import leggauss
from numpy.polynomial import Polynomial as Poly
from explore import apply_R
from gal import galerkin, make_basis

xg2, wg2 = leggauss(3000); xg2=0.5*(xg2+1); wg2=0.5*wg2

def check(N,pmin=3):
    vals,psi,vec,cand=galerkin(N,pmin)
    lr,k=cand[0]; c=vec[:,k].real
    M=len(psi)
    Zt=np.linspace(0.02,0.98,17)
    f=np.array([sum(c[j]*psi[j](Z) for j in range(M)) for Z in Zt])
    Rf=np.array([sum(c[j]*apply_R(psi[j].deriv(3),Z,xg2,wg2) for j in range(M)) for Z in Zt])
    print("N=%d pmin=%d lam=%.6f"%(N,pmin,lr))
    print(" f:", np.round(f,6))
    print(" R[f]/f:", np.round(Rf/f,4))

check(14,3)
check(14,4)