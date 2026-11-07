import numpy as np
from eta_transform import transform_op


if __name__ == "__main__":
    print("eta-transform p=1,1.2,1.5 : high N")
    for p in [1.0, 1.2, 1.5]:
        Ns, l1 = [], []
        for N in [320, 640, 960, 1280]:
            xc, A = transform_op(N, p)
            r = np.sort(-np.linalg.eigvals(A).real)
            r = r[r > 1e-3]
            Ns.append(N); l1.append(r[0])
        Ns = np.array(Ns, float); l1 = np.array(l1)
        A = np.vstack([np.ones(3), 1 / Ns[-3:], 1 / Ns[-3:] ** 2]).T
        s, *_ = np.linalg.lstsq(A, l1[-3:], rcond=None)
        print("p=%.1f  N=320..1280: %s  -> extrap %.5f" %
              (p, np.round(l1, 5), s[0]))