import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.linalg import null_space

exec(open('<round-runs>/runs/MP-40/explore13.py').read().split('print("Correct')[0])

K = 22
G, H, Q, NQ = build(K)
ev, V = np.linalg.eig(np.linalg.solve(H, G))
cands = [( -ev[k].real, k) for k in range(len(ev)) if abs(ev[k].imag) < 1e-6*max(1, abs(ev[k])) and -ev[k].real > 0]
cands.sort()
for lam1, k in cands[:4]:
    print("lam=%.6f" % lam1)
lam1, k1 = cands[0]
a1 = V[:, k1]
polys = [to01(chebT_poly(k)) for k in range(K+1)]
C = np.zeros((3, K+1))
for k, p in enumerate(polys):
    C[0, k] = np.polynomial.polynomial.polyval(0.0, p)
    C[1, k] = np.polynomial.polynomial.polyval(1.0, deriv(p, 1))
    C[2, k] = np.polynomial.polynomial.polyval(1.0, deriv(p, 2))
Z = null_space(C)
Bmat = np.zeros((len(XI), K+1))
for k, p in enumerate(polys):
    Bmat[:, k] = np.polynomial.polynomial.polyval(XI, p)
Bn = Bmat @ Z
R = np.linalg.lstsq(Bn, Q, rcond=None)[0]
coef = Z @ R @ a1


def modef(x):
    return np.polynomial.polynomial.polyval(x, coef)


f1 = modef(1.0)
print("mode f(1)=", f1)
print("x   f(x)/f(1)")
for x in [0.5, 0.75, 0.9]:
    print(" %0.2f %.8f" % (x, modef(x)/f1))
print("near x=0:")
for x in [1e-3, 2e-3, 5e-3, 1e-2, 2e-2, 5e-2]:
    print("  x=%.0e  f/f(1)=%.6e" % (x, modef(x)/f1))
xs = np.array([1e-3, 2e-3, 5e-3, 1e-2, 2e-2])
mv = np.abs(modef(xs))
print("power fit near 0:", np.polyfit(np.log(xs), np.log(mv), 1))
print("near x=1:  1-f(x)/f(1) vs (1-x)")
for x in [0.999, 0.995, 0.99, 0.98, 0.95]:
    print("  1-x=%.0e  1-f/f1=%.6e" % (1-x, 1-modef(x)/f1))
xe = np.array([1e-3, 2e-3, 5e-3, 1e-2, 2e-2])
print("power fit near 1:", np.polyfit(np.log(xe), np.log(np.abs(1-modef(1-xe)/f1)), 1))