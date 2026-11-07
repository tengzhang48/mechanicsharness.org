import numpy as np
from numpy.polynomial import Polynomial as Poly
from numpy.polynomial.legendre import leggauss
from explore import apply_R
from scipy.linalg import eig

xg, wg = leggauss(1000); xg=0.5*(xg+1); wg=0.5*wg

def phi0_poly():
    Mcon = np.zeros((6,6)); rhs=np.array([0,1,0,0,0,0.0])
    for k in range(6):
        p = Poly([0]*k+[1.0])
        Mcon[0,k]=p(0); Mcon[1,k]=p(1); Mcon[2,k]=p.deriv(1)(1)
        Mcon[3,k]=p.deriv(2)(1); Mcon[4,k]=p.deriv(3)(1); Mcon[5,k]=p.deriv(3)(0)
    return Poly(np.linalg.solve(Mcon,rhs))

phi0 = phi0_poly()

def make_basis3(N, pmin):
    psi=[phi0]
    for k in range(N):
        psi.append(Poly([0]*(pmin+k)+[1.0]) * Poly([1.0,-4,6,-4,1.0]))
    return psi

def run(N,pmin):
    psi=make_basis3(N,pmin); M=len(psi)
    Zc,_=leggauss(M); Zc=0.5*(Zc+1)
    Rmat=np.zeros((M,M)); Mmat=np.zeros((M,M))
    for i,Z in enumerate(Zc):
        for j in range(M):
            Rmat[i,j]=apply_R(psi[j].deriv(3),Z,xg,wg)
            Mmat[i,j]=psi[j](Z)
    lam,vec=eig(Rmat,Mmat)
    cand=sorted([(l.real,k) for k,l in enumerate(lam) if abs(l.imag)<1e-7*max(1,abs(l.real)) and l.real>1e-6])
    return [c[0] for c in cand[:4]], psi, vec, cand

for pmin in [2,3,4,5]:
    for N in [8,12,16]:
        vals,_,_,_ = run(N,pmin)
        print("pmin=%d N=%2d  lams:"%(pmin,N), np.round(vals,5))