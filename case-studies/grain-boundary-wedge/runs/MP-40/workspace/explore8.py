import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.linalg import eig

exec(open('<round-runs>/runs/MP-40/explore7.py').read().split('for K in')[0])

G, H = opmatrix(16)
print("G symmetry resid:", np.linalg.norm(G-G.T)/np.linalg.norm(G))
ev, V = eig(G, -H)
idx = np.argsort(np.abs(ev.imag) + np.where(ev.real < 0, 1e6, 0) + 0*ev.real)
# pick the eigenvector with smallest positive real eigenvalue
cand = [(v.real, k) for k, v in enumerate(ev) if abs(v.imag) < 1e-6*max(1, abs(v.real)) and abs(v.real) > 1e-8]
cand.sort()
print("candidate lambda:", [round(c[0], 4) for c in cand[:6]])
lam, k0 = cand[0]
a = V[:, k0].real
xs = np.linspace(1e-4, 0.2, 8)
vals = []
for x in xs:
    vals.append(sum(a[k]*np.polynomial.polynomial.polyval(x, cheb_basis_poly(k)) for k in range(16)))
vals = np.array(vals)
# normalize by value at x=0.2
print("x:", np.round(xs, 4))
print("sigma:", np.round(vals/vals[-1], 5))
slope = np.polyfit(np.log(xs), np.log(np.abs(vals)), 1)
print("fitted power p near 0:", slope[0])
print("lambda1 =", lam)