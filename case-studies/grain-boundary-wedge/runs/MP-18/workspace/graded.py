import numpy as np

# ---- Reuse the validated exact kernel antiderivatives from march.py ----
from march import _prim, _P1, _yP1, _P2, _yP2, _P3, _yP3


def Smat_graded(xcols, nodes):
    """M_ij = (1/pi) PV Int_0^1 S(x_i,y) phi_j(y) dy on an arbitrary grid.
    xcols : interior collocation nodes (len N)
    nodes : full node array length N+2 including endpoints x=0 and x=1.
    Hat basis phi_j supported on [nodes[j], nodes[j+2]] for j=0..N-1.
    Uses exact kernel antiderivatives identical to march.Smat."""
    N = len(xcols)
    M = np.zeros((N, N))
    for i in range(N):
        z = xcols[i]
        for j in range(N):
            a, b, c = nodes[j], nodes[j + 1], nodes[j + 2]
            val = (_prim(b, z, -a / (b - a), 1.0 / (b - a))
                   - _prim(a, z, -a / (b - a), 1.0 / (b - a)))
            val += (_prim(c, z, c / (c - b), -1.0 / (c - b))
                    - _prim(b, z, c / (c - b), -1.0 / (c - b)))
            M[i, j] = val / np.pi
    return M


def Smat_direct_quad(xcols, nodes, ngauss=200):
    """Reference: direct Gauss quadrature with PV handling, calling the *actual*
    kernel S(x,y) from the problem statement (not the antiderivative)."""
    from numpy.polynomial.legendre import leggauss
    N = len(xcols)
    M = np.zeros((N, N))

    def S(z, y):
        return 1.0 / (y - z) + 1.0 / (y + z) - 2.0 * z * (y - z) / (y + z) ** 3

    gx, gw = leggauss(ngauss)
    for i in range(N):
        z = xcols[i]
        for j in range(N):
            a, b, c = nodes[j], nodes[j + 1], nodes[j + 2]
            tot = 0.0
            for lo, hi in ((a, b), (b, c)):
                ys = 0.5 * (hi - lo) * gx + 0.5 * (hi + lo)
                w = 0.5 * (hi - lo) * gw
                for yv, wv in zip(ys, w):
                    if abs(yv - z) < 1e-13:
                        continue
                    phi = (yv - lo) / (b - lo) if hi <= b else (hi - yv) / (hi - b)
                    tot += wv * S(z, yv) * phi
            M[i, j] = tot / np.pi
    return M


def _sens(a, b, c, y):
    """Sensitivity of exact quadrature: d/dy of the product antiderivative? No.
    Instead just evaluate the integrand at y for sanity."""
    return None


if __name__ == "__main__":
    # validation on a uniform grid vs march.Smat
    import march
    N = 12
    h = 1.0 / (N + 1)
    xcols = np.arange(1, N + 1) * h
    nodes = np.concatenate(([0.0], xcols, [1.0]))
    Ma = Smat_graded(xcols, nodes)
    Mb = march.Smat(N)
    print("graded==uniform march.Smat ?", np.allclose(Ma, Mb, atol=1e-12),
          " max diff", np.max(np.abs(Ma - Mb)))
    # validate against direct Gaussian quadrature of the raw kernel S(x,y)
    Mc = Smat_direct_quad(xcols, nodes, ngauss=300)
    print("antiderivative vs direct quad: max diff", np.max(np.abs(Ma - Mc)),
          " rel", np.max(np.abs(Ma - Mc)) / np.max(np.abs(Ma)))