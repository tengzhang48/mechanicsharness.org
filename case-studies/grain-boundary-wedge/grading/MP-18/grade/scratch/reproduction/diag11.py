import numpy as np, march
from march import Smat, D1, D2

for N in [40, 45, 65, 85]:
    M = Smat(N); D1m = D1(N); D2m = D2(N)
    K = -(D2m @ M @ D1m)
    ev = np.linalg.eigvals(K)
    re = np.sort(ev.real)[::-1]
    print("N=%d top reals:" % N, np.round(re[:6], 6))
    print("   #zero(<1e-8):", np.sum(np.abs(re) < 1e-8))