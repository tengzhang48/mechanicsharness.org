import numpy as np
from fd_nonuniform import build_D1D2
from phys_oper import M_phys
from d3_nonuniform import build_D3_nonuniform
import grids


def solve_D3(nodes):
    xcols = np.asarray(nodes[1:-1], float)
    D3 = build_D3_nonuniform(nodes)
    M = M_phys(xcols, nodes)
    A = -(M @ D3)
    ev = np.linalg.eigvals(A)
    return -ev.real


if __name__ == "__main__":
    print("=== -M D3 with CORRECT nonuniform D3, uniform grid ===")
    for N in [20, 40, 80]:
        x = np.linspace(0, 1, N + 2)[1:-1]
        nodes = np.concatenate(([0.0], x, [1.0]))
        r = solve_D3(nodes)
        rs = np.sort(r)
        print("N=%3d  smallest |rates|: %s" % (N, np.round(rs[np.argsort(np.abs(rs))][:6], 5)))