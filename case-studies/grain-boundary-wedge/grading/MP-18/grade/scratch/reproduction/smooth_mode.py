import numpy as np
from fd_nonuniform import build_D1D2
from mfast import M_phys_fast
from d3_nonuniform import build_D3_nonuniform


def build_grid(N, stretch):
    """Custom grid via a mapping with a resolution parameter `stretch`:
    x(t) = t + stretch * t*(1-t)*sin? no. Use x(t)=t^q with q>1 -> cluster at 0.
    But to avoid the 1/h^2 blowup, we use a *mild* clustering and higher N."""
    t = np.linspace(0, 1, N + 2)
    x = t ** stretch
    x[0], x[-1] = 0.0, 1.0
    return x


def smooth_mode(N, stretch=1.0, op='chain'):
    x = build_grid(N, stretch)
    nodes = x
    xcols = nodes[1:-1]
    M = M_phys_fast(xcols, nodes)
    if op == 'chain':
        D1, D2 = build_D1D2(nodes)
        A = -(M @ D1 @ D2)
    else:
        A = -(M @ build_D3_nonuniform(nodes))
    w, V = np.linalg.eig(A)
    r = -w.real
    # find the mode whose shape is closest to the smooth positive profile f(0.5)~0.75
    best = None
    for i in range(len(r)):
        v = V[:, i].real
        if abs(v[-1]) < 1e-9:
            continue
        vv = v / v[-1]
        f05 = np.interp(0.5, xcols, vv)
        if f05 > 0.4 and f05 < 1.2:
            score = abs(r[i])
            if best is None or abs(f05 - 0.75) < best[3]:
                best = (r[i], f05, vv, abs(f05 - 0.75))
    return (best[0], best[1]) if best else (np.nan, np.nan)


if __name__ == "__main__":
    print("smooth mode (f(0.5)~0.75) eigenvalue vs N and stretch, chain op")
    for stretch in [1.0, 1.3, 1.6]:
        print("-- stretch=%.1f --" % stretch)
        for N in [40, 80, 160, 320]:
            lam, f05 = smooth_mode(N, stretch, 'chain')
            print("   N=%3d  lam=%.6f  f(0.5)=%.5f" % (N, lam, f05))