import numpy as np
from scipy.linalg import eig
from solver2 import *

def run(N):
    x, V, Vinv, W, D1v, D2v, D3v, A, B = assemble(N)
    mu, Z = eig(A, B)
    out = []
    for k in range(len(mu)):
        m = mu[k]
        if abs(m.imag) > 1e-8 * max(1.0, abs(m.real)) or m.real >= 0 or abs(m) > 1e7:
            continue
        v = Z[:, k]
        if abs(v[0]) < 1e-12:  # normalize f value at x=1
            continue
        out.append((m.real, np.real(v / v[0])))
    out.sort(key=lambda t: t[0])
    return x, out, (V, Vinv, W, D1v, D2v, D3v)

for N in [80, 120, 160, 220, 280]:
    x, out, aux = run(N)
    print(f"N={N}: {len(out)} decaying real modes; largest mu (slowest decay):")
    for m, v in out[:4]:
        lam = -m
        # mode value at 0.5 (grid interp) and residual sanity: f'(1),f''(1),f(0)
        print(f"   lam={lam:.8f}   f(0.5)~{np.interp(0.5, x, v):+.6f}")
