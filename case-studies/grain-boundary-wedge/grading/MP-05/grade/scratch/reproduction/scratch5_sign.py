import numpy as np
from core import *

N = 32
lam, V, res, A, B, xc = solve_eig(N)
mu = -lam  # eigenvalues of A a = mu B a, i.e. K[rho] = mu f
# pick the smallest-|mu| real one
cand = [(abs(m.real), i) for i, m in enumerate(mu) if abs(m.imag) < 1e-8 and m.real != 0]
cand.sort()
i0 = cand[0][1]
a = np.real(V[:, i0])
m0 = np.real(mu[i0])
print("mu =", m0)

# rho_hat(x) on fine grid
xs = np.linspace(0.001, 0.999, 400)
Q = cheb_basis(xs, N)
rho = np.sqrt(1 - xs) * (Q @ a)
# f from f_matrix at these xs
fvals = f_matrix(xs, N) @ a
# brute K[rho]: PV + images with fine quadrature
t, w = np.polynomial.legendre.leggauss(3000)
xi = (t + 1) / 2; ww = w / 2
Qi = cheb_basis(xi, N)
rhoi = np.sqrt(1 - xi) * (Qi @ a)
for x in [0.2, 0.5, 0.8]:
    ker = 1.0/(x-xi) - 1.0/(x+xi) - 2*xi*(x-xi)/(x+xi)**3
    pv = np.log(x/(1-x))
    reg = np.sum(ww * rhoi * (ker - 1.0/(x-xi)))
    Krho = reg + rhoi@xi*0 + pv*rho_at(x, a, N)
    fv = f_at(x, a, N)
    print(f"x={x}: brute K[rho]={Krho:.8f}  A*basis*f?  mu*f={m0*fv:.8f}  ratio={Krho/fv if fv!=0 else 0:.6f}")

def rho_at(x, a, N):
    return np.sqrt(1-x) * (cheb_basis(np.array([x]), N) @ a)[0]

def f_at(x, a, N):
    return (f_matrix(np.array([x]), N) @ a)[0]
