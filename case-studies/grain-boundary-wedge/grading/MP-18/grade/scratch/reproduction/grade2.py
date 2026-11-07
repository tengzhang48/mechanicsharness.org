import numpy as np
from fd_nonuniform import build_D1D2
from phys_oper import M_phys


def nodes_geometric(N, d0):
    """N+2 nodes on [0,1], first spacing d0, geometric growth until spacing
    reaches ~1/N, then uniform."""
    # geometric: d_k = d0 * r^k ; find r so that after n0 steps spacing ~ uniform
    # Simpler: build with a smooth stretch that has first spacing exactly d0.
    M = N + 2
    # use a sinh-like clustering: x_k = d0 * k + (1-d0*(M-1)) * (k/(M-1))^? 
    # We use a blended map: x(t) = a*t + (1-a)*t^p with t uniform, a small.
    return None


def nodes_map(N, a, p):
    """x(t) = a t + (1-a) t^p,  t = k/N uniform, k=0..N+1."""
    t = np.linspace(0, 1, N + 2)
    x = a * t + (1 - a) * t ** p
    x[0] = 0.0
    x[-1] = 1.0
    return x


def lam1(nodes):
    xcols = np.asarray(nodes[1:-1], float)
    D1, D2 = build_D1D2(nodes)
    M = M_phys(xcols, nodes)
    A = -(D2 @ M @ D1)
    ev = np.linalg.eigvals(A)
    return -np.max(ev.real), np.sort(-ev.real)[::-1]


if __name__ == "__main__":
    print("map x=a t + (1-a) t^p ; test convergence")
    for (a, p) in [(0.5, 2.0), (0.7, 3.0), (0.3, 2.0)]:
        print("-- a=%.1f p=%.1f --" % (a, p))
        for N in [40, 80, 160, 240]:
            nodes = nodes_map(N, a, p)
            lam, top = lam1(nodes)
            print("   N=%3d lam1=%.6f  top5=%s  min spacing=%.2e" %
                  (N, lam, np.round(top[:3], 4), np.diff(nodes).min()))