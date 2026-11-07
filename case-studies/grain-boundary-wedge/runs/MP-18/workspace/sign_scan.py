import numpy as np
from fd_nonuniform import build_D1D2

# Generic nonlocal operator from a kernel with the exact antiderivative machinery,
# but allow sign variations of the three terms to identify the decaying one.
# S(x,y) = a/(x-y) + b/(x+y) + c * (-2 y (x-y)/(x+y)^3)
# Use numeric PV quadrature (robust) to build M for any (a,b,c).

from numpy.polynomial.legendre import leggauss


def M_generic(xcols, nodes, a, b, c, ngauss=500):
    N = len(xcols)
    M = np.zeros((N, N))
    gx, gw = leggauss(ngauss)

    def Sper(z, y):
        return (a / (z - y) + b / (z + y)
                + c * (-2.0 * y * (z - y) / (z + y) ** 3))

    for i in range(N):
        z = xcols[i]
        for j in range(N):
            lo0, mid, hi = nodes[j], nodes[j + 1], nodes[j + 2]
            tot = 0.0
            for lo, hi in ((lo0, mid), (mid, hi)):
                pts = [lo, hi]
                if lo < z < hi:
                    pts = [lo, z, hi]
                for k in range(len(pts) - 1):
                    l, hg = pts[k], pts[k + 1]
                    ys = 0.5 * (hg - l) * gx + 0.5 * (hg + l)
                    w = 0.5 * (hg - l) * gw
                    phi = (ys - lo0) / (mid - lo0) if hi == mid else (hi - ys) / (hi - mid)
                    vals = Sper(z, ys)
                    vals = np.where(np.abs(ys - z) < 1e-14, 0.0, vals)
                    tot += np.sum(w * vals * phi)
            M[i, j] = tot / np.pi
    return M


if __name__ == "__main__":
    import itertools
    N = 30
    x = np.linspace(0, 1, N + 2)[1:-1]
    nodes = np.concatenate(([0.0], x, [1.0]))
    D1, D2 = build_D1D2(nodes)
    for (a, b, c) in [(1, -1, 1), (-1, -1, 1), (1, 1, 1), (-1, 1, 1),
                      (1, -1, -1), (-1, 1, -1), (1, 1, -1), (-1, -1, -1)]:
        M = M_generic(x, nodes, a, b, c)
        A = -(M @ D1 @ D2)
        ev = np.linalg.eigvals(A)
        print("a=%+d b=%+d c=%+d  max re=%+.4e  min re=%+.4e" %
              (a, b, c, np.max(ev.real), np.min(ev.real)))