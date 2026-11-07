import numpy as np
from fd_nonuniform import build_D1D2
from phys_oper import M_phys
import grids


def lam1_chain(nodes):
    """lambda1 of A = -D2 M D1  (equivalent to -M D1 D2 via similarity)."""
    xcols = np.asarray(nodes[1:-1], float)
    D1, D2 = build_D1D2(nodes)
    M = M_phys(xcols, nodes)
    A = -(D2 @ M @ D1)
    ev = np.linalg.eigvals(A)
    return -np.max(ev.real), ev


if __name__ == "__main__":
    print("=== chain -D2 M D1, PHYSICAL kernel, uniform ===")
    for N in [20, 40, 80, 160]:
        x = np.linspace(0, 1, N + 2)[1:-1]
        nodes = np.concatenate(([0.0], x, [1.0]))
        lam, ev = lam1_chain(nodes)
        print("N=%3d  lam1=%.8f  top=%s" %
              (N, lam, np.round(np.sort(-ev.real)[::-1][:4], 5)))