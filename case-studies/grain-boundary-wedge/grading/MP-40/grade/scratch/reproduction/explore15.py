exec(open('<round-runs>/runs/MP-40/explore13.py').read().split('print("Correct')[0])

import numpy as np

for K in [18, 20, 22, 24, 26, 28, 30]:
    G, H, Q, NQ = build(K)
    M = np.linalg.solve(H, G)
    ev, V = np.linalg.eig(M)
    lam = -ev.real
    order = np.argsort(lam + np.where(lam < 0, 1e9, 0))
    best = [k for k in order if abs(ev[k].imag) < 1e-6*max(1, abs(ev[k])) and -ev[k].real > 0][:4]
    vals = [(-ev[k].real, V[:, k]) for k in best]
    vals.sort(key=lambda t: t[0])
    line = "K=%2d " % K
    for lam1, v1 in vals:
        line += " %.4f" % lam1
    print(line)

# mode shape for K=26 evaluated via basis polynomials
K = 26
G, H, Q, NQ = build(K)
ev, V = np.linalg.eig(np.linalg.solve(H, G))
cands = [( -ev[k].real, k) for k in range(len(ev)) if abs(ev[k].imag) < 1e-6*max(1, abs(ev[k])) and -ev[k].real > 0]
cands.sort()
k1 = cands[0][1]
a1 = V[:, k1]
# rebuild polynomials and Z,R to evaluate mode
polys = [to01(chebT_poly(k)) for k in range(K+1)]
C = np.zeros((3, K+1))
for k, p in enumerate(polys):
    C[0, k] = np.polynomial.polynomial.polyval(0.0, p)
    C[1, k] = np.polynomial.polynomial.polyval(1.0, deriv(p, 1))
    C[2, k] = np.polynomial.polynomial.polyval(1.0, deriv(p, 2))
Z = null_space(C)
# recover R from Q: Q = Bn@R
Bmat = np.zeros((len(XI), K+1))
for k, p in enumerate(polys):
    Bmat[:, k] = np.polynomial.polynomial.polyval(XI, p)
Bn = Bmat @ Z
R = np.linalg.lstsq(Bn, Q, rcond=None)[0]     # R s.t. Bn@R=Q
coef = Z @ R @ a1                              # coefficients in polys
def modef(x):
    return np.polynomial.polynomial.polyval(x, coef)
xs = np.array([0.5, 0.75, 0.9, 1.0])
print("mode shape (normalized f(1)=1):")
m = modef(xs)/modef(1.0)
for x, v in zip(xs, m):
    print("   f(%.2f)=%.6f" % (x, v))
xsm = np.array([0.002, 0.005, 0.01, 0.02, 0.05])
mv = np.abs(modef(xsm))
pfit = np.polyfit(np.log(xsm), np.log(mv), 1)
print("near-0 power fit p =", pfit[0], " coeff:", np.exp(pfit[1]))