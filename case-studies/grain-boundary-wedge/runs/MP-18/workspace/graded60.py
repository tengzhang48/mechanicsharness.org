import numpy as np
from d3_nonuniform import build_D3_nonuniform
from mfast import M_phys_fast
import grids


def mode60_graded(N, p, a=None):
    if a is None:
        nodes = grids.graded_nodes(N, p)
    else:
        t = np.linspace(0, 1, N + 2)
        nodes = a * t + (1 - a) * t ** p
        nodes[0] = 0.0
        nodes[-1] = 1.0
    xcols = nodes[1:-1]
    D3 = build_D3_nonuniform(nodes)
    M = M_phys_fast(xcols, nodes)
    A = -(M @ D3)
    r = np.sort(-np.linalg.eigvals(A).real)
    return r[np.argmin(np.abs(r - 60))]


if __name__ == "__main__":
    print("graded x=eta^p, track mode ~60")
    for p in [2.0, 3.0]:
        print("-- p=%.1f --" % p)
        for N in [40, 80, 160, 320, 640]:
            print("   N=%3d  lam~60=%.6f" % (N, mode60_graded(N, p)))
    print("blended x=a t+(1-a)t^p")
    for (a, p) in [(0.5, 2.0), (0.7, 2.0)]:
        print("-- a=%.1f p=%.1f --" % (a, p))
        for N in [40, 80, 160, 320, 640]:
            print("   N=%3d  lam~60=%.6f" % (N, mode60_graded(N, p, a=a)))