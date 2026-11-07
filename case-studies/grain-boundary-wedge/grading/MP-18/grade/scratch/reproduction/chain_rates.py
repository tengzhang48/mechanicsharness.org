import numpy as np
from solve_D3n import solve_D3
from fd_nonuniform import build_D1D2
from phys_oper import M_phys
from chain import lam1_chain


def chain_rates(N):
    x = np.linspace(0, 1, N + 2)[1:-1]
    nodes = np.concatenate(([0.0], x, [1.0]))
    xcols = nodes[1:-1]
    D1, D2 = build_D1D2(nodes)
    M = M_phys(xcols, nodes)
    A = -(M @ D1 @ D2)
    ev = np.linalg.eigvals(A).real
    return np.sort(ev)[::-1]        # decay rates, largest first


if __name__ == "__main__":
    print("=== chain operator: slowest 4 decay rates vs N ===")
    for N in [40, 80, 160, 280]:
        r = chain_rates(N)
        print("N=%3d  top rates: %s" % (N, np.round(r[:4], 5)))