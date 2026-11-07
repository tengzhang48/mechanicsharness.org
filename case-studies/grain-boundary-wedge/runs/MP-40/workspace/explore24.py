exec(open('<round-runs>/runs/MP-40/explore23.py').read().split('for D in')[0])

res1 = {}
for D in [22, 24, 26, 28, 30, 32, 34, 36, 38, 40]:
    lam, nodes, Z, On, Vn = run(D)
    res1[D] = lam[0]
    print("D=%2d  lam1=%.5f  lam2=%.5f" % (D, lam[0], lam[1]), flush=True)

print("\nlambda1 convergence:")
Ds = sorted(res1)
for a, b, c in zip(Ds, Ds[1:], Ds[2:]):
    d1 = res1[b]-res1[a]; d2 = res1[c]-res1[b]
    ait = res1[c] - d2*d2/(d2-d1)
    print("  D=%2d,%2d,%2d  d1=%.5f d2=%.5f  Aitken lam_inf=%.5f" % (a, b, c, d1, d2, ait))

# fit lam(D)=linf+C*D^-p on last three
import numpy as np
l1, l2, l3 = res1[Ds[-3]], res1[Ds[-2]], res1[Ds[-1]]
D1, D2, D3 = Ds[-3], Ds[-2], Ds[-1]


def f(p):
    u, v, w = D1**(-p), D2**(-p), D3**(-p)
    return (u-v)/(v-w)


from scipy.optimize import brentq
plo, phi = 0.3, 4.0
try:
    p = brentq(lambda p: f(p) - (l1-l2)/(l2-l3), plo, phi)
    u, v, w = D1**(-p), D2**(-p), D3**(-p)
    C = (l1-l2)/(u-v)
    linf = l3 - C*w
    print("\n  fit: p=%.4f  C=%.3f  lam_inf=%.5f" % (p, C, linf))
except Exception as e:
    print("fit failed", e)