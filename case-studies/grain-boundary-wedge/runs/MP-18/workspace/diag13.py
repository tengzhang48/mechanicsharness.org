import numpy as np, march
from march import Smat, D1, D2

def smoothness(v, x):
    # measure fraction of energy in the smooth (least-oscillatory) part:
    # compare with its projection onto a low-order polynomial (deg<=3) in x
    A = np.vstack([x**0, x**1, x**2, x**3]).T
    coef, *_ = np.linalg.lstsq(A, v, rcond=None)
    proj = A @ coef
    return np.linalg.norm(proj) / np.linalg.norm(v)

print("N    eig    -re     smoothness   shape(0.1,0.5,1)")
for N in [30, 40, 50, 60, 80, 100, 140, 200]:
    M = Smat(N); D1m = D1(N); D2m = D2(N)
    K = -(D2m @ M @ D1m)
    ev, V = np.linalg.eig(K)
    x = np.arange(1, N + 1) / (N + 1)
    # candidate slow modes: top 4 real
    idx = np.argsort(ev.real)[::-1][:4]
    best = None
    for i in idx:
        v = V[:, i].real
        sm = smoothness(v, x)
        if best is None or sm > best[2]:
            best = (i, ev.real[i], sm, v)
    i, re, sm, v = best
    vv = v / v[-1]
    print("%4d  %+.6f  sm=%.3f   %s" % (N, re, sm, np.round(np.interp([0.1,0.5,1.0], x, vv), 4)))