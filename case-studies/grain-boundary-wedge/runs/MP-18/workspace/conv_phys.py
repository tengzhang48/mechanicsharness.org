import numpy as np
from fd_nonuniform import build_D1D2
from phys_oper import M_phys
import grids


def leading(nodes, Mfun=M_phys):
    xcols = np.asarray(nodes[1:-1], float)
    D1, D2 = build_D1D2(nodes)
    M = Mfun(xcols, nodes)
    A = -(M @ D1 @ D2)          # s-form  s_tau = -M D1 D2 s
    ev = np.linalg.eigvals(A)
    j = np.argmax(ev.real)
    return -ev.real[j], ev, M, D1, D2


if __name__ == "__main__":
    print("=== UNIFORM grid, physical kernel, s-form (-M D1 D2) ===")
    for N in [10, 20, 40, 80]:
        x = np.linspace(0, 1, N + 2)[1:-1]
        nodes = np.concatenate(([0.0], x, [1.0]))
        lam, ev, M, D1, D2 = leading(nodes)
        print("N=%3d  lambda1=%.8f  top=%s" %
              (N, lam, np.round(np.sort(-ev.real)[::-1][:4], 5)))