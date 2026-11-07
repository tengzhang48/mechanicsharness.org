import numpy as np
from numpy.polynomial import Polynomial as Poly
from numpy.polynomial.legendre import leggauss
from explore import apply_R
from scipy.linalg import eig

xg, wg = leggauss(800); xg=0.5*(xg+1); wg=0.5*wg

# Build basis:  phi0: degree5 poly with phi(0)=0, phi(1)=1, phi'(1)=phi''(1)=phi'''(1)=0, phi'''(0)=0
Mcon = np.zeros((6,6))
rhs = np.array([0,1,0,0,0,0.0])  # phi(0),phi(1),phi'(1),phi''(1),phi'''(1),phi'''(0)
for k in range(6):
    p = Poly([0]*k+[1.0])
    Mcon[0,k]=p(0); Mcon[1,k]=p(1); Mcon[2,k]=p.deriv(1)(1)
    Mcon[3,k]=p.deriv(2)(1); Mcon[4,k]=p.deriv(3)(1); Mcon[5,k]=p.deriv(3)(0)
c0 = np.linalg.solve(Mcon, rhs)
phi0 = Poly(c0)
print("phi0 coeffs", np.round(c0,6))

def make_basis2(N):
    psi=[phi0]
    for k in range(0,N):
        psi.append(Poly([0]*(k+2)+[1.0]) * Poly([1.0,-4,6,-4,1.0]))  # Z^{k+2}(Z-1)^4
    return psi

def build(N):
    psi = make_basis2(N)
    M=len(psi)
    Zc,_ = leggauss(M); Zc=0.5*(Zc+1)
    Rmat=np.zeros((M,M)); Mmat=np.zeros((M,M))
    for i,Z in enumerate(Zc):
        for j in range(M):
            Rmat[i,j]=apply_R(psi[j].deriv(3),Z,xg,wg)
            Mmat[i,j]=psi[j](Z)
    return Rmat,Mmat,Zc,psi

for N in [14]:
    Rmat,Mmat,Zc,psi=build(N)
    lam,vec=eig(Rmat,Mmat)
    cand=sorted([(l.real,k) for k,l in enumerate(lam) if abs(l.imag)<1e-8*max(1,abs(l.real)) and l.real>1e-6])
    print("N=%d"%N, "first real lambdas:", np.round([c[0] for c in cand[:5]],6))
    lr,k=cand[0]; c=vec[:,k]
    Zt=np.linspace(0.01,0.99,25)
    f=np.array([sum(c[j]*psi[j](Z) for j in range(len(psi))) for Z in Zt])
    Rf=np.array([sum(c[j]*apply_R(psi[j].deriv(3),Z,xg,wg) for j in range(len(psi))) for Z in Zt])
    print("lam=",lr)
    print("f:", np.round(f,5))
    print("R[f]/f:", np.round(Rf/f,4))