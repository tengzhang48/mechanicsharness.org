import numpy as np
from fd_nonuniform import build_D1D2
from statement_oper import Smat_statement
import march, grids


def eig_leading(nodes, Mmat):
    xcols = np.asarray(nodes[1:-1], float)
    D1, D2 = build_D1D2(nodes)
    A = -(D2 @ Mmat @ D1)
    ev = np.linalg.eigvals(A)
    # slowest decay = max real part
    return ev


def run(kernel, N, p):
    nodes = grids.graded_nodes(N, p)
    xcols = nodes[1:-1]
    if kernel == 'stmt':
        M = Smat_statement(xcols, nodes)
    else:
        M = march.Smat(N) if abs(p - 1.0) < 1e-12 else None
    return nodes, M


if __name__ == "__main__":
    # uniform grid with statement kernel, watch the leading eigenvalue vs N
    print("=== statement kernel, UNIFORM grid ===")
    for N in [10, 20, 40, 80]:
        x = np.linspace(0, 1, N + 2)[1:-1]
        nodes = np.concatenate(([0.0], x, [1.0]))
        M = Smat_statement(x, nodes)
        ev = eig_leading(nodes, M)
        print("N=%3d  max re=%+.6e   top5=%s" %
              (N, np.max(ev.real), np.round(np.sort(ev.real)[::-1][:5], 4)))