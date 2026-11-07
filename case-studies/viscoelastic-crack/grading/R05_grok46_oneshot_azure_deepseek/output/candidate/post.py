"""Post-processing: normalized fields, bulk dissipation integral, and checks."""
import numpy as np
import scipy.sparse.linalg as spla
from scipy.interpolate import RectBivariateSpline
import solve

LAM = solve.LAM
DELTA_BAR = solve.DELTA_BAR
DBAR = DELTA_BAR


def solve_state(V, a, eps, npanel, ny, L=14.0):
    A, b, X, Y, Dx1, Dx2, Dy1, Dy2, nyy, nW, K0 = solve.assemble(
        V, a, eps, L=L, npanel=npanel, ny=ny)
    sol = spla.spsolve(A.tocsc(), b)
    W = sol[:nW].reshape(len(X), nyy + 1)
    S2solved = sol[nW:]
    C = np.pi * (1.0 - LAM) / LAM
    Aavg, skip = solve.build_avg_matrix(X, K0)
    WX = Dx1 @ W
    WY = W @ Dy1.T
    S2 = C * (WY - Aavg @ WY) - C * (a / np.pi) * skip[:, None]
    S1 = C * (WX - Aavg @ WX)
    return dict(X=X, Y=Y, W=W, WX=WX, WY=WY, S1=S1, S2=S2,
                S2solved=S2solved, K0=K0, C=C, Dx1=Dx1, Dy1=Dy1)


def bulk_ratio(V, a, eps, npanel, ny, L=14.0):
    F = solve_state(V, a, eps, npanel, ny, L)
    integ = F['S1'] ** 2 + F['S2'] ** 2
    half = np.trapezoid(np.trapezoid(integ, F['X'], axis=0), F['Y'])
    full = 2.0 * half
    return full / (DBAR * (1.0 - LAM) * V * np.pi), F


def sample(V, a, eps, sample_X, sample_yh, npanel, ny, L=14.0):
    F = solve_state(V, a, eps, npanel, ny, L)
    X, Y = F['X'], F['Y']
    spl = {k: RectBivariateSpline(X, Y, F[k], kx=3, ky=3)
           for k in ['WX', 'WY', 'S1', 'S2']}
    rows = []
    for Xs in sample_X:
        for yh in sample_yh:
            Ys = np.pi * yh
            WX = float(spl['WX'](Xs, Ys)[0, 0])
            WY = float(spl['WY'](Xs, Ys)[0, 0])
            S1 = float(spl['S1'](Xs, Ys)[0, 0])
            S2 = float(spl['S2'](Xs, Ys)[0, 0])
            rows.append(dict(X=Xs, y_over_h=yh,
                             T1=np.pi * WX + S1, T2=np.pi * WY + S2,
                             E1=np.pi * WX, E2=np.pi * WY))
    return rows, F


if __name__ == "__main__":
    import sys
    V, a, eps = 0.1, 0.5158215, 0.7683462
    r, F = bulk_ratio(V, a, eps, (240, 150, 240), 150)
    print('V=%g a=%.6f eps=%.6f R=%.6f R-1=%.6f  integral=%.6f' %
          (V, a, eps, a * a / DBAR, a * a / DBAR - 1, r))