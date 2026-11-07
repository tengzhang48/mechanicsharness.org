import numpy as np
from fd_nonuniform import build_D1D2
from mfast import M_phys_fast
from d3_nonuniform import build_D3_nonuniform


def rates_chain(N):
    x = np.linspace(0, 1, N + 2)[1:-1]
    nodes = np.concatenate(([0.0], x, [1.0]))
    D1, D2 = build_D1D2(nodes)
    M = M_phys_fast(x, nodes)
    A = -(M @ D1 @ D2)
    return np.sort(-np.linalg.eigvals(A).real)      # ascending: slowest first


def rates_D3(N):
    x = np.linspace(0, 1, N + 2)[1:-1]
    nodes = np.concatenate(([0.0], x, [1.0]))
    D3 = build_D3_nonuniform(nodes)
    M = M_phys_fast(x, nodes)
    A = -(M @ D3)
    return np.sort(-np.linalg.eigvals(A).real)[::-1]


if __name__ == "__main__":
    Ns = [40, 80, 160, 320, 640, 960]
    print("=== CHAIN operator: slowest rates ===")
    c1, c2 = [], []
    for N in Ns:
        r = rates_chain(N)
        c1.append(r[0]); c2.append(r[1])
        print("N=%4d  lam1=%.6f  lam2=%.6f" % (N, r[0], r[1]))
    print("=== D3 operator: mode near 60 and near 1070 ===")
    d60, d1070 = [], []
    for N in Ns:
        r = rates_D3(N)
        d60.append(r[np.argmin(np.abs(r - 60))])
        d1070.append(r[np.argmin(np.abs(r - 1070))])
        print("N=%4d  ~60=%.6f  ~1070=%.5f" % (N, d60[-1], d1070[-1]))