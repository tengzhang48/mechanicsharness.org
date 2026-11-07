import numpy as np
import solve

LAM = solve.LAM
DELTA_BAR = solve.DELTA_BAR


def fields(V, a, eps, npanel, ny, L=14.0):
    A, b, X, Y, Dx1, Dx2, Dy1, Dy2, nyy, nW, K0 = solve.assemble(V, a, eps, L=L, npanel=npanel, ny=ny)
    import scipy.sparse.linalg as spla
    sol = spla.spsolve(A.tocsc(), b)
    W = sol[:nW].reshape(len(X), nyy + 1)
    S2_solved = sol[nW:]
    C = np.pi * (1.0 - LAM) / LAM
    Aavg, skip = solve.build_avg_matrix(X, K0)
    WX = Dx1 @ W
    WY = W @ Dy1.T
    S2 = C * (WY - Aavg @ WY) - C * (a / np.pi) * skip[:, None]
    S1 = C * (WX - Aavg @ WX)
    T2 = np.pi * WY + S2
    T1 = np.pi * WX + S1
    E1 = np.pi * WX
    E2 = np.pi * WY
    return dict(X=X, Y=Y, W=W, WX=WX, WY=WY, S1=S1, S2=S2, T1=T1, T2=T2, E1=E1, E2=E2,
                S2_solved=S2_solved, K0=K0, C=C, Dy1=Dy1, Dx1=Dx1, Yg=Y)


def bulk_ratio(V, a, eps, npanel, ny, L=14.0, ret=False):
    F = fields(V, a, eps, npanel, ny, L)
    X, Y = F['X'], F['Y']
    integ = F['S1'] ** 2 + F['S2'] ** 2
    # integrate over upper half strip: integral dX dY, then double for full strip
    IX = np.trapezoid(integ, X, axis=0)          # (ny+1,)
    half = np.trapezoid(IX, Y)
    full = 2.0 * half
    ratio = full / (DELTA_BAR * (1.0 - LAM) * V * np.pi)
    if ret:
        return ratio, F
    return ratio


def sample_fields(V, a, eps, sample_X, sample_yh, npanel, ny, L=14.0):
    from scipy.interpolate import RectBivariateSpline
    F = fields(V, a, eps, npanel, ny, L)
    X, Y = F['X'], F['Y']
    out = []
    spl = {}
    for key in ['WX', 'WY', 'S1', 'S2']:
        spl[key] = RectBivariateSpline(X, Y, F[key], kx=3, ky=3)
    for Xs in sample_X:
        for yh in sample_yh:
            Ys = np.pi * yh
            WX = float(spl['WX'](Xs, Ys)[0, 0])
            WY = float(spl['WY'](Xs, Ys)[0, 0])
            S1 = float(spl['S1'](Xs, Ys)[0, 0])
            S2 = float(spl['S2'](Xs, Ys)[0, 0])
            E1 = np.pi * WX
            E2 = np.pi * WY
            T1 = E1 + S1
            T2 = E2 + S2
            out.append(dict(X=Xs, y_over_h=yh, T1=T1, T2=T2, E1=E1, E2=E2))
    return out


if __name__ == "__main__":
    V = 0.1
    a, eps = 0.5155, 0.7466
    R = 1.0 + a * a / DELTA_BAR
    print('expected R=%.6f  R-1=%.6f' % (R, R - 1))
    for np_, ny in [((160, 100, 160), 100), ((240, 150, 240), 150), ((320, 200, 320), 200)]:
        rr, F = bulk_ratio(V, a, eps, np_, ny, ret=True)
        chk = np.max(np.abs(F['S2'][:, 0] - F['S2_solved']))
        print('res=%s ratio=%.6f  R-1=%.6f  diff=%.2e  S2chk=%.2e' %
              (np_, rr, R - 1, rr - (R - 1), chk))