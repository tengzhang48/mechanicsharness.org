import numpy as np
import solve


def state(V, a, eps, npanel, ny, nlag=64):
    f, Kt, T0, W, S2, X, Y, K0 = solve.residuals(V, a, eps, npanel=npanel, ny=ny, nlag=nlag)
    return f, Kt, W, S2, X, Y, K0


def newton(V, a0, eps0, npanel, ny, tol=1e-7, itmax=8, verbose=True):
    a, eps = a0, eps0
    for it in range(itmax):
        f0, K0, *_ = state(V, a, eps, npanel, ny)
        r0 = np.array([f0 - 0.1, K0])
        if np.linalg.norm(r0) < tol:
            break
        da, de = 1e-3 * max(1, abs(a)), 1e-3 * max(1, abs(eps))
        f1, K1, *_ = state(V, a + da, eps, npanel, ny)
        f2, K2, *_ = state(V, a, eps + de, npanel, ny)
        J = np.array([[f1 - f0, f2 - f0], [K1 - K0, K2 - K0]]) / np.array([[da, de], [da, de]])
        try:
            step = np.linalg.solve(J, -r0)
        except np.linalg.LinAlgError:
            break
        a, eps = a + step[0], eps + step[1]
        if verbose:
            print('  it%d a=%.6f eps=%.6f res=%.2e,%.2e' % (it, a, eps, r0[0], r0[1]))
    return a, eps


if __name__ == "__main__":
    for V, g in [(0.001, (0.30, 0.45)), (0.1, (0.52, 0.78)), (10.0, (1.5, 1.0))]:
        print('V=', V)
        a, eps = newton(V, g[0], g[1], (300, 190, 300), 190)
        print('  root a=%.6f eps=%.6f' % (a, eps))