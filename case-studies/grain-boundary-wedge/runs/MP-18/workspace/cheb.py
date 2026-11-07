import numpy as np
from numpy.polynomial.chebyshev import chebpts1
from numpy.polynomial.legendre import leggauss


def cheb_grid(N):
    """Chebyshev-Gauss-Lobatto points mapped to [0,1]."""
    k = np.arange(N + 1)
    x = 0.5 * (1 - np.cos(np.pi * k / N))     # increasing, x[0]=0, x[N]=1
    return x


def cheb_diff_matrix(x):
    """Spectral differentiation matrix D (d/dx) on the mapped Chebyshev grid x."""
    N = len(x) - 1
    xi = np.cos(np.pi * np.arange(N + 1) / N)   # [-1,1], decreasing; x=(1-xi)/2
    c = np.ones(N + 1)
    c[0] = 2.0
    c[N] = 2.0
    X = np.tile(xi, (N + 1, 1)).T
    dX = X - X.T
    Dxi = np.outer(c, 1.0 / c) / (dX + np.eye(N + 1))
    Dxi = Dxi - np.diag(Dxi.sum(axis=1))
    # d/dx = (dxi/dx) d/dxi = (-2) d/dxi
    D = -2.0 * Dxi
    return D


if __name__ == "__main__":
    for N in [10, 20, 40]:
        x = cheb_grid(N)
        D = cheb_diff_matrix(x)
        f = x ** 3
        err = np.max(np.abs(D @ f - 3 * x ** 2))
        print("N=%d cheb D err on x^3: %.3e   x[0..3]=%s" % (N, err, np.round(x[:4], 5)))