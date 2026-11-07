import numpy as np
from conv_phys import leading
import grids


def track_uniform(Nmax_list):
    print("N     lambda1        lam*N/lnN   lam*sqrt(N)   lam*N")
    for N in Nmax_list:
        x = np.linspace(0, 1, N + 2)[1:-1]
        nodes = np.concatenate(([0.0], x, [1.0]))
        lam, ev, M, D1, D2 = leading(nodes)
        print("%4d  %.7f   %.4f   %.4f   %.4f" %
              (N, lam, lam * N / np.log(N), lam * np.sqrt(N), lam * N))


def track_graded(Nmax_list, p):
    print("p=%.1f  N     lambda1        lam*N/lnN   lam*sqrt(N)" % p)
    for N in Nmax_list:
        nodes = grids.graded_nodes(N, p)
        lam, ev, M, D1, D2 = leading(nodes)
        print("      %4d  %.7f   %.4f   %.4f" %
              (N, lam, lam * N / np.log(N), lam * np.sqrt(N)))


if __name__ == "__main__":
    import sys
    mode = sys.argv[1] if len(sys.argv) > 1 else "uniform"
    if mode == "uniform":
        track_uniform([20, 40, 80, 160, 240])
    else:
        p = float(sys.argv[2]) if len(sys.argv) > 2 else 3.0
        track_graded([20, 40, 80, 160, 240, 400], p)