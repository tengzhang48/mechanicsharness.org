exec(open('<round-runs>/runs/MP-40/explore23.py').read().split('def run(D)')[0])

import mpmath as mp
mp.mp.dps = 50
import numpy as np


def make_basis(M):
    basis = []
    w = np.array([1.0, -1.0])
    one_plus = np.polyadd(np.array([1.0]), w)
    one_plus = np.polyadd(one_plus, np.convolve(w, w))
    basis.append(np.convolve(np.array([0.0, 1.0]), one_plus))
    for j in range(3, M+1):
        pw = np.array([1.0])
        for _ in range(j):
            pw = np.convolve(pw, w)
        basis.append(np.convolve(np.array([0.0, 1.0]), pw))
    return basis


def run_hp(M):
    basis = make_basis(M)
    n = len(basis)
    m = n
    nodes = [mp.mpf(0.5)*(1-mp.cos(mp.pi*(i+1)/(m+1))) for i in range(m)]
    On = mp.matrix(m, n); Vn = mp.matrix(m, n)
    for i in range(m):
        xi = nodes[i]
        for k, p in enumerate(basis):
            acc = mp.mpf(0)
            for d, c in enumerate(p):
                if abs(c) > 1e-14:
                    acc += mp.mpf(c)*Nval(d, xi)
            On[i, k] = acc
            Vn[i, k] = mp.polyval(list(p[::-1]), xi)
    A = On; B = Vn
    C = mp.inverse(B)*A
    ev = mp.eig(C, left=False, right=False)
    lam = sorted([-complex(e).real for e in ev if abs(complex(e).imag) < 1e-6*max(1, abs(complex(e))) and -complex(e).real > 1e-6])
    return lam, nodes, basis, A, B


prev = {}
for M in [12, 16, 20, 24, 28, 32, 36, 40]:
    lam, nodes, basis, A, B = run_hp(M)
    prev[M] = lam[:3]
    print("M=%2d  " % M + "  ".join("%.5f" % v for v in lam[:4]), flush=True)

print()
for key in [0, 1, 2]:
    vv = [(M, prev[M][key]) for M in sorted(prev) if len(prev[M]) > key]
    print("quantity", key)
    for a, b in zip(vv, vv[1:]):
        print("   M=%2d val=%.5f  diff=%.5f" % (b[0], b[1], b[1]-a[1]))