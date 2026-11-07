import numpy as np
from conv_phys import leading


def modes_uniform(N, k=5):
    x = np.linspace(0, 1, N + 2)[1:-1]
    nodes = np.concatenate(([0.0], x, [1.0]))
    lam, ev, M, D1, D2 = leading(nodes)
    A = -(M @ D1 @ D2)
    w, V = np.linalg.eig(A)
    idx = np.argsort(-w.real)[:k]
    out = []
    for i in idx:
        v = V[:, i].real
        v = v / v[-1] if abs(v[-1]) > 1e-12 else v / np.max(np.abs(v))
        out.append((-w[i].real, v))
    return x, out


if __name__ == "__main__":
    print("Slowest modes on uniform grids; check convergence of shape at x=0.5")
    for N in [40, 80, 160, 320]:
        x, modes = modes_uniform(N, 3)
        line = "N=%3d  " % N
        for (lam, v) in modes:
            v05 = np.interp(0.5, x, v)
            line += "lam=%.5f f(0.5)=%.5f | " % (lam, v05)
        print(line)