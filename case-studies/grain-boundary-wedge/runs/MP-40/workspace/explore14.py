exec(open('<round-runs>/runs/MP-40/explore13.py').read().split('print("Correct')[0])

import numpy as np
from numpy.linalg import eig as geig

K = 22
G, H, Q, NQ = build(K)
ev, V = geig_eig = np.linalg.eig(np.linalg.solve(H, G))
# eigenvalue of solve(H,G) = eigenvalue of pencil G a = mu H a ; our lambda = -mu
mu = ev
lam = -mu.real
order = np.argsort(lam)
print("smallest magnitudes of |lam|, real parts:")
for k in order[:8]:
    print("   lam = %.6f  %+.6f i   imag/|lam|=%.2e" % (ev[k].real*-1, (-ev[k]).imag, abs(ev[k].imag)/max(1e-30, abs(ev[k]))))

# first real positive mode
cands = [( -ev[k].real, k) for k in range(len(ev)) if abs(ev[k].imag) < 1e-6*max(1.0, abs(ev[k])) and -ev[k].real > 0]
cands.sort()
lam1, k1 = cands[0]
lam2, k2 = cands[1]
print("\nlam1=%.6f lam2=%.6f ratio=%.4f" % (lam1, lam2, lam2/lam1))
a1 = V[:, k1]
mode = NQ @ a1          # sigma values at XI nodes (sinc NQ=N applied to Q; Q a1 = eigenfunction)
# actually eigenfunction = Q @ a1
sig = Q @ a1
xs = np.array([0.0, 0.1, 0.25, 0.5, 0.75, 0.9, 1.0])
vals = []
for x in xs:
    # interpolate from XI (Gauss nodes) via polynomial fit of degree K-3
    deg = min(K-2, 30)
    coef = np.polyfit(XI, sig, deg)
    vals.append(np.polyval(coef, x))
vals = np.array(vals)
print("x  :", xs)
print("sig:", np.round(vals, 6))
print("normalized f(1)=1:", np.round(vals/vals[-1], 6))
print("f(0.5) =", vals[3]/vals[-1])