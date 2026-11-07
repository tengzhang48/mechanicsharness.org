import numpy as np
from numpy.polynomial import Polynomial as Poly
from numpy.polynomial.legendre import leggauss
from explore import apply_R, make_basis

xg, wg = leggauss(600); xg=0.5*(xg+1); wg=0.5*wg

def build(N, pts='gl'):
    psi = make_basis(N)
    P3 = [p.deriv(3) for p in psi]
    M = len(psi)
    if pts=='gl':
        zc,_ = leggauss(M); Zc = 0.5*(zc+1)
    else:
        Zc = np.linspace(0.02,0.98,M)
    Rmat = np.zeros((M,M)); Mmat=np.zeros((M,M))
    for i,Z in enumerate(Zc):
        for j in range(M):
            Rmat[i,j]=apply_R(P3[j],Z,xg,wg)
            Mmat[i,j]=psi[j](Z)
    return Rmat,Mmat,Zc,psi

N=16
Rmat,Mmat,Zc,psi = build(N)
from scipy.linalg import eig
lam, vec = eig(Rmat, Mmat)
cand = sorted([(l.real,k) for k,l in enumerate(lam) if abs(l.imag)<1e-8*max(1,abs(l.real)) and l.real>1e-6])
lr,k = cand[0]
c=vec[:,k]
res = Rmat@c - lam[k]*Mmat@c
print("lam=",lam[k], "||res||/||Rmat c||=", np.linalg.norm(res)/np.linalg.norm(Rmat@c))
# check on fine grid
Zt = np.linspace(0.02,0.98,25)
f = np.array([sum(c[j]*psi[j](Z) for j in range(N+1)) for Z in Zt])
Rf = np.array([sum(c[j]*apply_R(psi[j].deriv(3),Z,xg,wg) for j in range(N+1)) for Z in Zt])
ratio = Rf/f
print("R[f]/f on grid:", np.round(ratio,4))