import numpy as np
from d3_nonuniform import build_D3_nonuniform
from mfast import M_phys_fast
from fd_nonuniform import build_D1D2


def mod60(N, op='d3'):
    x = np.linspace(0, 1, N + 2)[1:-1]
    nodes = np.concatenate(([0.0], x, [1.0]))
    M = M_phys_fast(x, nodes)
    if op == 'd3':
        A = -(M @ build_D3_nonuniform(nodes))
    else:
        D1, D2 = build_D1D2(nodes)
        A = -(M @ D1 @ D2)
    r = np.sort(-np.linalg.eigvals(A).real)
    return r[np.argmin(np.abs(r - 60))]


if __name__ == "__main__":
    Ns = [120, 240, 480, 960]
    out = {}
    for op in ['chain', 'd3']:
        out[op] = np.array([mod60(N, op) for N in Ns])
    np.save("mod60_data.npy", (Ns, out['chain'], out['d3']))
    print("N        chain        d3")
    for N, c, d in zip(Ns, out['chain'], out['d3']):
        print("%5d  %.6f   %.6f" % (N, c, d))