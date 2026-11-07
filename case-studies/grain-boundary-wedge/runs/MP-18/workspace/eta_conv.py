import numpy as np
from eta_transform import transform_op


if __name__ == "__main__":
    print("eta-transform (p=1, uniform): genuine leading eigenvalue vs N")
    Ns = []
    l1 = []
    for N in [40, 80, 160, 320, 640, 960]:
        xc, A = transform_op(N, 1.0)
        r = np.sort(-np.linalg.eigvals(A).real)
        r = r[r > 1e-3]
        Ns.append(N)
        l1.append(r[0])
        print("   N=%4d  lam1=%.6f" % (N, r[0]))
    Ns = np.array(Ns, float)
    l1 = np.array(l1)
    A = np.vstack([np.ones(4), 1 / Ns[-4:], 1 / Ns[-4:] ** 2, 1 / Ns[-4:] ** 3]).T
    s, *_ = np.linalg.lstsq(A, l1[-4:], rcond=None)
    print("extrapolate (L+a/N+b/N^2+c/N^3) = %.6f" % s[0])