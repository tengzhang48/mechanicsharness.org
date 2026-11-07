import numpy as np
from numpy.polynomial.legendre import leggauss

# Exact antiderivatives of  S_phys(x,y) = 1/(y-x) + 1/(x+y) + 2y(x-y)/(x+y)^3
# (sympy-verified):
#   F1 = Int S_phys dy
#   F2 = Int F1 dy
#   Int (alpha + beta y) S_phys dy = alpha F1 + beta (y F1 - F2)


def F1(x, y):
    Lm = np.log(np.abs(y - x) + 1e-300)
    Lp = np.log(np.abs(y + x) + 1e-300)
    return -2.0 * x * (3.0 * y + 2.0 * x) / (y + x) ** 2 - (Lp - Lm)


def F2(x, y):
    Lm = np.log(np.abs(y - x) + 1e-300)
    Lp = np.log(np.abs(y + x) + 1e-300)
    return (-2.0 * x ** 2 + (y + x) * (y * Lm - x * Lm - (y + 7.0 * x) * Lp)) / (y + x)


def _int_lin(x, lo, hi, alpha, beta):
    def G(y):
        return alpha * F1(x, y) + beta * (y * F1(x, y) - F2(x, y))
    return G(hi) - G(lo)


def M_phys(xcols, nodes):
    N = len(xcols)
    M = np.zeros((N, N))
    for i in range(N):
        z = xcols[i]
        for j in range(N):
            a, b, c = nodes[j], nodes[j + 1], nodes[j + 2]
            val = _int_lin(z, a, b, -a / (b - a), 1.0 / (b - a))
            val += _int_lin(z, b, c, c / (c - b), -1.0 / (c - b))
            M[i, j] = val / np.pi
    return M


def M_phys_numeric(xcols, nodes, ngauss=600):
    """Reference by direct PV-aware Gauss quadrature of S_phys."""
    N = len(xcols)
    M = np.zeros((N, N))
    gx, gw = leggauss(ngauss)

    def Sper(z, yy):
        return 1.0 / (yy - z) + 1.0 / (z + yy) + 2.0 * yy * (z - yy) / (z + yy) ** 3

    for i in range(N):
        z = xcols[i]
        for j in range(N):
            lo0, mid, hi = nodes[j], nodes[j + 1], nodes[j + 2]
            tot = 0.0
            for lo, hg in ((lo0, mid), (mid, hi)):
                pts = [lo, hg]
                if lo < z < hg:
                    pts = [lo, z, hg]
                for k in range(len(pts) - 1):
                    l, hgt = pts[k], pts[k + 1]
                    ys = 0.5 * (hgt - l) * gx + 0.5 * (hgt + l)
                    w = 0.5 * (hgt - l) * gw
                    phi = (ys - lo0) / (mid - lo0) if hg == mid else (hi - ys) / (hi - mid)
                    vals = Sper(z, ys)
                    vals = np.where(np.abs(ys - z) < 1e-14, 0.0, vals)
                    tot += np.sum(w * vals * phi)
            M[i, j] = tot / np.pi
    return M


if __name__ == "__main__":
    N = 12
    x = np.linspace(0, 1, N + 2)[1:-1]
    nodes = np.concatenate(([0.0], x, [1.0]))
    Ma = M_phys(x, nodes)
    Mn = M_phys_numeric(x, nodes, ngauss=800)
    print("M_phys exact vs numeric: max diff", np.max(np.abs(Ma - Mn)),
          " rel", np.max(np.abs(Ma - Mn)) / np.max(np.abs(Ma)))
    # also compare to overall-negative of statement kernel
    from statement_oper import Smat_statement
    Ms = -Smat_statement(x, nodes)
    print("M_phys vs -Smat_statement: max diff", np.max(np.abs(Ma - Ms)))