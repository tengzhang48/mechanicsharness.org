import numpy as np
from d3_nonuniform import build_D3_nonuniform
from mfast import M_phys_fast
from fd_nonuniform import build_D1D2


def modes(N, op):
    x = np.linspace(0, 1, N + 2)[1:-1]
    nodes = np.concatenate(([0.0], x, [1.0]))
    M = M_phys_fast(x, nodes)
    if op == 'd3':
        A = -(M @ build_D3_nonuniform(nodes))
    else:
        D1, D2 = build_D1D2(nodes)
        A = -(M @ D1 @ D2)
    w, V = np.linalg.eig(A)
    r = -w.real
    return x, r, V


if __name__ == "__main__":
    print("=" * 68)
    print("FINAL CONVERGENCE TABLE (uniform grid, physical kernel)")
    print("=" * 68)
    print("op=D3 operator  -M D3")
    print("   N       lambda1(~60)    f1(0.5)      f1(0.25)   f1(0.75)")
    for N in [120, 240, 480, 960]:
        x, r, V = modes(N, 'd3')
        i = np.argmin(np.where(r > 1, r, np.inf))       # slowest genuine
        v = V[:, i].real
        v = v / v[-1]
        print("  %5d   %10.5f    %9.6f   %8.5f   %8.5f" %
              (N, r[i], np.interp(0.5, x, v), np.interp(0.25, x, v),
               np.interp(0.75, x, v)))