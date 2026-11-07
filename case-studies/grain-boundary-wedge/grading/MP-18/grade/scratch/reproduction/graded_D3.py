import numpy as np
from fd_nonuniform import build_D1D2
from phys_oper import M_phys
from d3_nonuniform import build_D3_nonuniform
import grids


def rates_D3(nodes):
    """Return DECAY RATES (positive) sorted descending: slowest first."""
    xcols = np.asarray(nodes[1:-1], float)
    D3 = build_D3_nonuniform(nodes)
    M = M_phys(xcols, nodes)
    A = -(M @ D3)
    ev = np.linalg.eigvals(A).real
    return np.sort(-ev)[::-1]


if __name__ == "__main__":
    print("=== -M D3 on graded grids (x=eta^p): track the ~60 mode ===")
    for p in [1.0, 1.5, 2.0, 2.5]:
        print("-- p=%.1f --" % p)
        for N in [40, 80, 160, 240]:
            nodes = grids.graded_nodes(N, p)
            r = rates_D3(nodes)
            near = r[np.argmin(np.abs(r - 60))]
            print("   N=%3d  nearest-60 = %.6f   top4=%s" %
                  (N, near, np.round(np.sort(r)[:4], 3)))