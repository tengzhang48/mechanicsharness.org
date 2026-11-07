import numpy as np
from d3_nonuniform import build_D3_nonuniform
from mfast import M_phys_fast
from fd_nonuniform import build_D1D2


def spec(N, op='d3'):
    x = np.linspace(0, 1, N + 2)[1:-1]
    nodes = np.concatenate(([0.0], x, [1.0]))
    M = M_phys_fast(x, nodes)
    if op == 'd3':
        D3 = build_D3_nonuniform(nodes)
        A = -(M @ D3)
    else:
        D1, D2 = build_D1D2(nodes)
        A = -(M @ D1 @ D2)
    return np.sort(-np.linalg.eigvals(A).real)


if __name__ == "__main__":
    for op in ['d3', 'chain']:
        N = 160
        r = spec(N, op)
        # list all genuine rates < 100
        small = r[r < 100]
        print("%s N=%d  rates < 100: %s" % (op, N, np.round(small, 4)))
        print("   #rates<100 =", len(small))