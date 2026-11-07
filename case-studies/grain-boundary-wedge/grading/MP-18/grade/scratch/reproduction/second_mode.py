import numpy as np
from conv_phys import leading


def second_mode(N):
    x = np.linspace(0, 1, N + 2)[1:-1]
    nodes = np.concatenate(([0.0], x, [1.0]))
    lam, ev, M, D1, D2 = leading(nodes)
    A = -(M @ D1 @ D2)
    w, V = np.linalg.eig(A)
    idx = np.argsort(-w.real)
    # first two slow modes
    res = []
    for i in idx[:2]:
        v = V[:, i].real
        v = v / v[-1] if abs(v[-1]) > 1e-9 else v / np.max(np.abs(v))
        res.append((-w[i].real, np.interp(0.5, x, v), v))
    return x, res


if __name__ == "__main__":
    print("N      lam1        f1(0.5)     lam2        f2(0.5)")
    for N in [20, 40, 80, 160, 320, 480]:
        x, res = second_mode(N)
        print("%4d  %.7f  %.7f   %.7f  %.7f" %
              (N, res[0][0], res[0][1], res[1][0], res[1][1]))