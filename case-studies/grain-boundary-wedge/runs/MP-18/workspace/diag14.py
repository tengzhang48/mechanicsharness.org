import numpy as np, march
from march import Smat, D1, D2

def smooth_v(N):
    M = Smat(N); D1m = D1(N); D2m = D2(N)
    K = -(D2m @ M @ D1m)
    ev, V = np.linalg.eig(K)
    x = np.arange(1, N + 1) / (N + 1)
    A = np.vstack([x**0, x**1]).T  # near-linear ramp; genuine mode is nearly linear
    idx = np.argsort(ev.real)[::-1][:5]
    best = None
    for i in idx:
        v = V[:, i].real
        sm = np.linalg.norm(A @ np.linalg.lstsq(A, v, rcond=None)[0]) / np.linalg.norm(v)
        if best is None or sm > best[2]:
            best = (i, ev.real[i], sm)
    return -best[1]

print("N     lambda1      lambda1*N/lnN   lambda1*N/ln(N)^2   lambda1*N^0.5")
for N in [30, 40, 60, 80, 120, 160, 240, 320, 480]:
    lam = smooth_v(N)
    print("%4d  %.7f   %.4f   %.4f   %.4f" %
          (N, lam, lam * N / np.log(N), lam * N / np.log(N)**2, lam * np.sqrt(N)))