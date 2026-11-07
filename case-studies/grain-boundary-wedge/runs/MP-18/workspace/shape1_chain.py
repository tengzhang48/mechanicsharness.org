import numpy as np
from fd_nonuniform import build_D1D2
from mfast import M_phys_fast


def mode1_chain(N):
    x = np.linspace(0, 1, N + 2)[1:-1]
    nodes = np.concatenate(([0.0], x, [1.0]))
    D1, D2 = build_D1D2(nodes)
    M = M_phys_fast(x, nodes)
    A = -(M @ D1 @ D2)
    w, V = np.linalg.eig(A)
    r = -w.real
    r = np.where(r > 1, r, np.inf)
    i = np.argmin(r)
    v = V[:, i].real
    v = v / v[-1]
    return x, v


if __name__ == "__main__":
    print("chain operator, lambda1 mode shape (f(1)=1)")
    for N in [240, 480, 960]:
        x, v = mode1_chain(N)
        print("N=%4d  f(0.1)=%.5f f(0.25)=%.5f f(0.5)=%.6f f(0.75)=%.5f" %
              (N, np.interp(0.1, x, v), np.interp(0.25, x, v),
               np.interp(0.5, x, v), np.interp(0.75, x, v)))