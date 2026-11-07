import numpy as np
from numpy.polynomial.legendre import leggauss


def build_D3_uniform(N, h=None):
    """Third-derivative operator on interior nodes x_i=i h, i=1..N, h=1/(N+1),
    with BCs s(0)=0, s'(1)=0, s''(1)=0, using ghost nodes at BOTH ends.
    Unknowns: s_1..s_N.
      left ghost s_{-1}: from central 3rd deriv stencil we need s_0 and s_{-1};
         s_0=0 and s_{-1} determined by ... we use the interior 5-point?
      Standard: interior 3rd deriv central: (s_{i+2}-2s_{i+1}+2s_{i-1}-s_{i-2})/(2h^3)
      left end i=1: use one-sided (needs s_{-1}); s_0=0 gives relation but s_{-1}
         unknown -> use 4-point one-sided on i=1,2,3 with s_0=0.
      right end i=N: BCs s'(1)=s''(1)=0 eliminate ghosts s_{N+1},s_{N+2}.
    """
    if h is None:
        h = 1.0 / (N + 1)
    D = np.zeros((N, N))
    for i in range(N):
        if i >= 2 and i <= N - 3:
            # central 3rd derivative: -s_{i-2}+2s_{i-1}-2s_{i+1}+s_{i+2} over 2h^3
            D[i, i - 2] = -1.0 / (2 * h ** 3)
            D[i, i - 1] = 2.0 / (2 * h ** 3)
            D[i, i + 1] = -2.0 / (2 * h ** 3)
            D[i, i + 2] = 1.0 / (2 * h ** 3)
        elif i == 0:
            # one-sided 3rd deriv at x_1 using s_0=0, s_1,s_2,s_3
            # f''' (~) = (-s_0 + 3 s_1 - 3 s_2 + s_3)/h^3 ; s_0=0
            D[i, i] = 3.0 / h ** 3
            D[i, i + 1] = -3.0 / h ** 3
            D[i, i + 2] = 1.0 / h ** 3
        elif i == 1:
            # second node one-sided: use s_0=0, s_1..s_4? fall back to central-ish
            # (-2 s_0 + ... ) use formula for f''' at x_2 with s_0=0:
            # ( -s_0 + 3s_1 -3s_2 + s_3 )? that's x_1. For x_2 use
            # ( s_0 - 2 s_1 + 0? ) ... just use central formula with s_0=0 ghost:
            D[i, i - 1] = -1.0 / (2 * h ** 3)    # s_0 term (value 0) dropped -> wait
            # central at i would need s_{-1}? no, i=1 central needs s_{-1},s_0,s_2,s_3
            # s_{-1} unknown -> fall back one-sided up to s_4
            D[i, :] = 0.0
            # f''' at x_2 from s_0..s_4 with s_0=0. Use undetermined coefficients.
            xs = np.array([0, h, 2 * h, 3 * h, 4 * h])
            k = 2
            coef = _fd_coeffs(xs, k, 3)
            for c, val in zip(coef, [0, 1, 2, 3, 4]):
                if val == 0:
                    continue
                D[i, val - 1] += c
        else:
            # near right end: BCs eliminate ghosts. For i=N-2,N-1 use one-sided.
            xs = np.array([(i - 2) * h, (i - 1) * h, i * h, (i + 1) * h, (i + 2) * h])
            # indices of nodes: i-2..i+2 -> interior indices i-3..i+1
            k = 2
            coef = _fd_coeffs(xs, k, 3)
            D[i, :] = 0.0
            for c, off in zip(coef, [-2, -1, 0, 1, 2]):
                j = i + off
                if j <= 0:
                    continue      # s_0=0
                if j > N:
                    # beyond x=1: eliminate via BCs
                    # s'(1)=s''(1)=0 -> s_{N+1}=s_{N-1}, s_{N+2}=s_{N-2}
                    if j == N + 1:
                        D[i, N - 1] += c
                        continue
                    if j == N + 2:
                        D[i, N - 2] += c
                        continue
                    raise RuntimeError
                D[i, j - 1] += c
    return D


def _fd_coeffs(xs, k, m):
    """Finite-difference coefficients for the m-th derivative at xs[k]."""
    import math
    n = len(xs)
    A = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            A[i, j] = (xs[j] - xs[k]) ** i
    b = np.zeros(n)
    b[m] = math.factorial(m)
    return np.linalg.solve(A, b)


if __name__ == "__main__":
    N = 10
    x = np.arange(1, N + 1) / (N + 1)
    D = build_D3_uniform(N)
    f = x ** 3
    print("D3 x^3 vs 6:", np.round(D @ f, 4))