"""Independent grading audit of the delivered solver; no submission files are changed."""
from pathlib import Path
import sys

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.linalg import eig
from scipy.special import eval_chebyu

sys.path.insert(0, str(Path(__file__).parent / "reproduction"))
import stmt_spec
import phys_oper
from d3_uniform import build_D3_uniform
from fd_nonuniform import build_D1D2
from mfast import M_phys_fast
from eta_transform import transform_op


def primitive_cos(n, a):
    n = np.asarray(n, dtype=float)
    return np.where(n == 0, a, np.sin(n * a) / np.where(n == 0, 1, n))


def traction_basis(x, n):
    """Integrate f'''=T_n(x)/sqrt(1-x^2), with f(0)=f'(1)=f''(1)=0."""
    x = np.atleast_1d(x)[:, None]
    n = np.atleast_1d(n)[None, :]
    a = np.arccos(x)

    def cos2_integral(a):
        return (0.5 * primitive_cos(n, a)
                + 0.25 * (primitive_cos(n - 2, a) + primitive_cos(n + 2, a)))

    cos1_integral = 0.5 * (primitive_cos(n - 1, a) + primitive_cos(n + 1, a))
    return 0.5 * (cos2_integral(np.pi / 2) - cos2_integral(a)
                  + 2 * x * cos1_integral - x ** 2 * primitive_cos(n, a))


def singular_collocation(N, quadrature=800):
    """Odd Chebyshev density; PV part exact, regular image term integrated in theta."""
    n = 2 * np.arange(N) + 1
    x = np.cos(np.arange(1, N + 1) * np.pi / (2 * N + 1))
    q, w = leggauss(quadrature)
    theta = (q + 1) * np.pi / 4
    w = w * np.pi / 4
    y = np.cos(theta)
    image = -2 * y[None, :] * (x[:, None] - y[None, :]) / (x[:, None] + y[None, :]) ** 3
    A = -np.pi * eval_chebyu(n[None, :] - 1, x[:, None])
    A += (image * w[None, :]) @ np.cos(theta[:, None] * n[None, :])
    B = traction_basis(x, n)
    rates, vectors = eig(-A, B)
    valid = np.flatnonzero((rates.real > 0) & (np.abs(rates.imag) < 1e-5))
    order = valid[np.argsort(rates[valid].real)]
    v = vectors[:, order[0]].real
    half = (traction_basis(0.5, n) @ v) / (traction_basis(1.0, n) @ v)
    return rates[order[:3]].real, half.item()


def main():
    print("Independent singular-integral collocation, same stated kernel and t0:")
    for N in [8, 16, 32, 64, 96]:
        rates, half = singular_collocation(N)
        print(f"N={N}: lambda={rates}, f(0.5)={half:.8f}")

    print("\nDelivered D3 manufactured-function checks:")
    for N in [10, 20, 40, 80, 160]:
        x = np.arange(1, N + 1) / (N + 1)
        D = build_D3_uniform(N)
        f = 1 - (1 - x) ** 3  # satisfies all three actual boundary conditions
        err = D @ f - 6
        singular = np.linalg.svd(D, compute_uv=False)
        _, _, Vh = np.linalg.svd(D)
        null = Vh[-1] / Vh[-1, -1]
        print(f"N={N}: BC-compatible cubic D3 error, max={max(abs(err)):.6g}, "
              f"last2={err[-2:]}; rank={np.linalg.matrix_rank(D)}, "
              f"sigma_min/sigma_max={singular[-1]/singular[0]:.3g}, "
              f"null f(0.5)={np.interp(0.5, x, null):.6f}")
        if N == 10:
            print("  D3 x^3:", D @ x ** 3)
            print("  D3 x^2:", D @ x ** 2)
            print("  last2 nonzero stencils:", [list(zip(np.flatnonzero(row), row[row != 0])) for row in D[-2:]])
        qnull = x * (2 - 1 / (N + 1) - x)
        print(f"  q=x*(2-h-x) null residual={max(abs(D @ qnull)):.3g}; "
              f"q'(1)={-1/(N+1):.6g}, q''(1)=-2")

    print("\nC6/C7 full spectra, original /pi convention:")
    for N in [20, 40, 80, 160]:
        x = np.arange(1, N + 1) / (N + 1)
        nodes = np.r_[0.0, x, 1.0]
        M = M_phys_fast(x, nodes)
        D1, D2 = build_D1D2(nodes)
        chain = -(D2 @ M @ D1)
        _, transformed = transform_op(N, 1.0)
        rchain = np.sort(-eig(chain, right=False).real)
        rtransformed = np.sort(-eig(transformed, right=False).real)
        mismatch = np.linalg.norm(transformed + M @ D1 @ D2) / np.linalg.norm(M @ D1 @ D2)
        print(f"N={N}: chain rates={rchain[:4]}; p=1 rates={rtransformed[:4]}; "
              f"p=1 same-grid operator relative mismatch={mismatch:.6g}")
    print("C5 stored original /pi rates multiplied by pi:",
          np.pi * np.array([-1054.505, -346.928, -59.927]))

    print("\nKernel and potential checks:")
    N = 12
    x = np.arange(1, N + 1) / (N + 1)
    nodes = np.r_[0.0, x, 1.0]
    M = stmt_spec.M_stmt(x, nodes)
    Mn = -phys_oper.M_phys_numeric(x, nodes, ngauss=800)
    print("C1 max relative matrix difference:", np.max(np.abs(M - Mn)) / np.max(np.abs(M)))
    for z, y in [(0.3, 0.7), (0.7, 0.3)]:
        S = 1 / (z - y) - 1 / (z + y) - 2 * y * (z - y) / (z + y) ** 3
        Kx = 1 / (z - y) - 1 / (z + y) - 2 * y * (z - y) / (z + y) ** 3
        print(f"z={z},y={y}: S={S}, dK/dx={Kx}")

    print("\nConstant slope sign test (PV analytic integral):")
    for z in [0.25, 0.5]:
        integral = -(stmt_spec.F1(z, 1.0) - stmt_spec.F1(z, 0.0))
        print(f"x={z}: integral S dy={integral:.8f}; "
              "u=y-1 has positive slope but negative opening, u=1-y has negative slope and positive opening")


if __name__ == "__main__":
    main()
