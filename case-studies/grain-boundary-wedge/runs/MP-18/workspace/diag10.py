import numpy as np, march
from march import Smat, D1, D2

print("N     lambda1 = -max re(K)    lambda1*N/lnN")
for N in [15, 25, 35, 45, 65, 85, 115, 155, 215]:
    M = Smat(N); D1m = D1(N); D2m = D2(N)
    K = -(D2m @ M @ D1m)
    ev = np.linalg.eigvals(K)
    lam1 = -np.max(ev.real)
    print("%4d   %.8f      %.4f" % (N, lam1, lam1 * N / np.log(N)))