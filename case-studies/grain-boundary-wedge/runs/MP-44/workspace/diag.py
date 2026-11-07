import numpy as np
from numpy.polynomial import Polynomial as Poly
from numpy.polynomial.legendre import leggauss
from explore import apply_R, make_basis, eigen

xg, wg = leggauss(400); xg=0.5*(xg+1); wg=0.5*wg
N=24
lam, vec, Zc, psi = eigen(N)
P3 = [p.deriv(3) for p in psi]
cand = [(l.real,k) for k,l in enumerate(lam) if abs(l.imag)<1e-6*abs(l.real) and l.real>1e-6]
cand.sort()
lr,k = cand[0]
c = vec[:,k]
f1 = sum(c[j]*psi[j](1.0) for j in range(len(psi)))
f = lambda Z: sum(c[j]*psi[j](Z) for j in range(len(psi)))/f1
Rf = lambda Z: sum(c[j]*apply_R(P3[j],Z,xg,wg) for j in range(len(psi)))/f1
print("lam1=",lr)
for Z in [0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95]:
    print("Z=%.2f f=%.6f  R[f]/f=%.6f"%(Z, f(Z), Rf(Z)/f(Z)))