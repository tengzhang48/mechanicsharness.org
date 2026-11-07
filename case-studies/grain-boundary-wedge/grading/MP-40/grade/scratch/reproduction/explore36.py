exec(open('<round-runs>/runs/MP-40/explore33.py').read().split('# cross-check')[0])

import numpy as np
import mpmath as mp
mp.mp.dps = 60

xw3 = np.array([0.0, 1.0, -3.0, 3.0, -1.0])


def leg_y(j):
    # ascending coeffs of P_j(y)
    if j == 0:
        return np.array([1.0])
    if j == 1:
        return np.array([0.0, 1.0])
    pm1 = np.array([1.0]); p = np.array([0.0, 1.0])
    for k in range(1, j):
        # (k+1)P_{k+1} = (2k+1) y P_k - k P_{k-1}
        ypk = np.convolve(np.array([0.0, 1.0]), p)
        nxt = ((2*k+1)*ypk - k*np.append(pm1, [0]*(len(ypk)-len(pm1))))
        nxt = nxt/(k+1)
        pm1, p = p, nxt
    return p


def leg_x(j):
    # P_j(2x-1) ascending coeffs in x
    cy = leg_y(j)
    out = np.array([0.0])
    vec = np.array([1.0])
    ypow = np.array([1.0])  # (2x-1)^0
    for m, a in enumerate(cy):
        if abs(a) > 1e-300:
            term = a*ypow
            if len(term) > len(out):
                out = np.append(out, np.zeros(len(term)-len(out)))
            out[:len(term)] += term
        ypow = np.convolve(ypow, np.array([-1.0, 2.0]))
    return out


def make_basis(N):
    b = [np.array([0.0, 3.0, -3.0, 1.0])]
    for j in range(N+1):
        b.append(np.convolve(xw3, leg_x(j)))
    return b


def study(N):
    basis = make_basis(N)
    n = len(basis); m = n
    nodes = [0.5*(1-np.cos(np.pi*(i+1)/(m+1))) for i in range(m)]
    On = mp.matrix(m, n); Vn = mp.matrix(m, n)
    for i, xi in enumerate(nodes):
        for k, p in enumerate(basis):
            On[i, k] = mp.mpf(Nphi_fast(xi, p))
            Vn[i, k] = mp.mpf(peval(p, xi))
    L = mp.inverse(Vn)*On
    ev = mp.eig(L, left=False, right=False)
    pos = sorted([-float(complex(e).real) for e in ev if abs(complex(e).imag) < 1e-9*max(1, abs(complex(e))) and -float(complex(e).real) > 1e-6])
    return basis, nodes, On, Vn, L, ev, pos


print("N   n   lam1        lam2        lam3", flush=True)
Ns = [8, 12, 16, 20, 24, 28, 32, 36, 40, 44]
res = {}
for N in Ns:
    basis, nodes, On, Vn, L, ev, pos = study(N)
    res[N] = (basis, nodes, On, Vn, L, ev, pos)
    print("%2d %3d  %s" % (N, len(basis), "  ".join("%.6f" % v for v in pos[:3])), flush=True)

lam1 = np.array([res[N][6][0] for N in Ns])
print("\nlam1:", list(zip(Ns, ["%.6f" % v for v in lam1])))
print("diffs:", ["%.5f" % v for v in np.diff(lam1)])
# Richardson fit lam1 ~ a + c/N^p over last 5
from scipy.optimize import curve_fit


def model(N, a, c, p):
    return a + c*np.array(N, float)**(-p)


try:
    popt, _ = curve_fit(model, Ns[-5:], lam1[-5:], p0=[12, 1, 2], maxfev=20000)
    print("Richardson a=%.6f c=%.4f p=%.4f" % tuple(popt))
except Exception as e:
    print("fit err", e)