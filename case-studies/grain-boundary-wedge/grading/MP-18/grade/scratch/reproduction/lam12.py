import numpy as np
from d3_nonuniform import build_D3_nonuniform
from mfast import M_phys_fast


def dsolve(N):
    x = np.linspace(0, 1, N + 2)[1:-1]
    nodes = np.concatenate(([0.0], x, [1.0]))
    M = M_phys_fast(x, nodes)
    A = -(M @ build_D3_nonuniform(nodes))
    w, V = np.linalg.eig(A)
    return x, -w.real, V


if __name__ == "__main__":
    print("N       lam1(~60)    lam2(~1070)")
    l1, l2 = [], []
    for N in [240, 480, 960, 1440, 1920]:
        x, r, V = dsolve(N)
        r = r[r > 1]     # drop spurious zero
        rs = np.sort(r)
        a = rs[0]
        b = rs[np.argmin(np.abs(rs - 1070))]
        l1.append(a); l2.append(b)
        print("%5d   %.6f    %.5f" % (N, a, b))
    Ns = np.array([240, 480, 960, 1440, 1920], float)
    for name, v in [("lam1", np.array(l1)), ("lam2", np.array(l2))]:
        A = np.vstack([np.ones(3), 1 / Ns[-3:], 1 / Ns[-3:] ** 2]).T
        s, *_ = np.linalg.lstsq(A, v[-3:], rcond=None)
        print("%s extrap (last3, L+a/N+b/N^2) = %.5f" % (name, s[0]))