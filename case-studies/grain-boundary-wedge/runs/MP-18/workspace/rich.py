import numpy as np
from chain import lam1_chain
from fd_nonuniform import build_D1D2
from phys_oper import M_phys


def chain_top(N, k=3):
    x = np.linspace(0, 1, N + 2)[1:-1]
    nodes = np.concatenate(([0.0], x, [1.0]))
    xcols = nodes[1:-1]
    D1, D2 = build_D1D2(nodes)
    M = M_phys(xcols, nodes)
    A = -(M @ D1 @ D2)
    ev = np.linalg.eigvals(A).real
    return np.sort(ev)[::-1][:k]     # decay rates (positive), slowest first


if __name__ == "__main__":
    Ns = [20, 40, 80, 160, 320, 480]
    mode1 = []
    mode2 = []
    for N in Ns:
        r = chain_top(N, 3)
        mode1.append(-r[0])
        mode2.append(-r[1])
    mode1 = np.array(mode1)
    mode2 = np.array(mode2)
    print("N       mode1        mode2")
    for N, a, b in zip(Ns, mode1, mode2):
        print("%4d   %.6f    %.6f" % (N, a, b))

    # Richardson: mode = L + c/N^p.  Use 3 finest points to solve for L.
    def richardson(Ns, vals, p):
        A = np.vstack([np.ones(3), np.array(Ns) ** (-p)]).T
        # 3 points, 2 unknowns -> least squares
        sol, *_ = np.linalg.lstsq(A, vals, rcond=None)
        return sol[0]
    for p in [0.5, 1.0, 1.5, 2.0]:
        print("p=%.1f  mode1_extrap=%.6f   mode2_extrap=%.6f" %
              (p, richardson(Ns[-3:], mode1[-3:], p),
               richardson(Ns[-3:], mode2[-3:], p)))