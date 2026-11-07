"""Independent FFT finite-crack mixed-boundary check (no factorization input).

Run: OPENBLAS_NUM_THREADS=1 python3 finite_check.py
The given w is the opening: V = integral w dx. Loading is unit resolvent
traction, not the Laplace transform 1/s of the step loading.
"""

import argparse
import time
import numpy as np
from scipy.sparse.linalg import LinearOperator, cg


def solve(n, box, s_values, tolerance=2e-10):
    C = 1 / (2 * (1 - .2))
    b = .89 * (1 - 2 * .2) / (1 - .2)
    beta = (.4 - .2) * (1 - .2) / (.89 * (1 - 2 * .2) * (1 - .4))
    eta = b * beta
    dx = box / n
    # Tips lie halfway between samples, avoiding sampled infinite crack flux.
    x = (np.arange(n) - n // 2 + .5) * dx
    indices = np.flatnonzero(np.abs(x) < 1)
    nc = len(indices)
    k = 2 * np.pi * np.abs(np.fft.rfftfreq(n, dx))

    def operator(symbol):
        def action(v):
            full = np.zeros((2, n))
            full[:, indices] = v.reshape(2, nc)
            modes = np.fft.rfft(full, axis=1)
            result = np.empty_like(modes)
            result[0] = symbol[0] * modes[0] + symbol[1] * modes[1]
            result[1] = symbol[1] * modes[0] + symbol[2] * modes[1]
            return np.fft.irfft(result, n=n, axis=1)[:, indices].ravel()
        return LinearOperator((2 * nc, 2 * nc), matvec=action, dtype=float)

    rhs = np.r_[np.ones(nc), np.zeros(nc)]
    # Same-grid elastic reference removes much of the sqrt-tip quadrature bias.
    a0 = C * (1 + eta) * k
    reference_symbol = (a0, np.zeros_like(k), np.ones_like(k))
    reference_inverse = (1 / np.maximum(a0, C / box), np.zeros_like(k), np.ones_like(k))
    reference, info = cg(operator(reference_symbol), rhs,
                         M=operator(reference_inverse), rtol=tolerance, atol=0)
    if info:
        raise RuntimeError(f"reference CG failed: {info}")
    v_reference = dx * reference[:nc].sum()
    v_exact = np.pi / (2 * C * (1 + eta))
    print(f"GRID box={box:g} N={n} dx={dx:.9g} crack_points={nc} "
          f"V0_grid={v_reference:.12g} V0_infinite={v_exact:.12g}", flush=True)
    previous = reference
    coefficients = []
    for s in s_values:
        start = time.monotonic()
        m = np.sqrt(k * k + s)
        r = s / (m + k)**2
        h = b * k / (m + k)
        lam = 2 * beta * C * s / b
        # q = sqrt(lam)*z. Equations are -S=1 and -p0/sqrt(lam)=0.
        a = C * k * (1 + eta * r) + lam * h * h / m
        cross = -np.sqrt(lam) * h / m
        d = 1 / m
        determinant = a * d - cross * cross
        safe = determinant.copy()
        safe[0] = C / box * d[0]
        inverse = (d / safe, -cross / safe, a / safe)
        inverse[2][0] = 1 / d[0]
        A = operator((a, cross, d))
        count = [0]

        def callback(_):
            count[0] += 1

        solution, info = cg(A, rhs, M=operator(inverse), x0=previous,
                            rtol=tolerance, atol=0, maxiter=1800, callback=callback)
        residual = np.linalg.norm(A @ solution - rhs) / np.linalg.norm(rhs)
        volume = dx * solution[:nc].sum()
        coefficient = 2 / np.sqrt(np.pi) * np.sqrt(s) * (volume - v_reference)
        coefficient_raw = 2 / np.sqrt(np.pi) * np.sqrt(s) * (volume - v_exact)
        print(f"RESULT s={s:g} V={volume:.12g} A_grid={coefficient:.9f} "
              f"A_raw={coefficient_raw:.9f} sqrt_s_dx={np.sqrt(s)*dx:.6g} "
              f"iterations={count[0]} relres={residual:.3g} info={info} "
              f"seconds={time.monotonic()-start:.2f}", flush=True)
        if info:
            raise RuntimeError(f"mixed CG failed: {info}")
        previous = solution
        coefficients.append(coefficient)
    return np.asarray(coefficients)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sizes", type=int, nargs="+", default=[8192, 16384, 32768, 65536])
    parser.add_argument("--box", type=float, default=32)
    parser.add_argument("--s", type=float, nargs="+", default=[100, 1000, 10000, 100000])
    args = parser.parse_args()
    print("Independent full mixed-boundary FFT check; expected A~2.4397; face-only A~1.354", flush=True)
    values = []
    for size in args.sizes:
        values.append(solve(size, args.box, args.s))
    if len(values) >= 2:
        ratio = args.sizes[-1] / args.sizes[-2]
        extrapolated = (ratio * values[-1] - values[-2]) / (ratio - 1)
        print("First-order dx extrapolation from last two grids (assumed, not guaranteed):", flush=True)
        for s, a in zip(args.s, extrapolated):
            print(f"EXTRAP s={s:g} A_dx0={a:.9f}", flush=True)
        if len(args.s) >= 3:
            t = np.asarray(args.s)**(-.25)
            # Empirical subleading powers, not used in the numerical solution.
            for degree in (1, 2):
                fit = np.polynomial.polynomial.polyfit(t, extrapolated, degree)
                error = np.max(np.abs(np.polynomial.polynomial.polyval(t, fit) - extrapolated))
                print(f"FIT empirical powers s^(-1/4) degree={degree} A_limit={fit[0]:.9f} "
                      f"max_fit_residual={error:.3g}", flush=True)
    print("Failure modes: finite periodic box; unresolved diffusion/tip layer when sqrt(s)*dx "
          "is large; midpoint sqrt-tip bias and singular flux; finite-s asymptotics; "
          "CG residual is not a discretization error bound. A_grid subtracts a same-grid "
          "undrained reference, not a factorization-derived correction.", flush=True)
