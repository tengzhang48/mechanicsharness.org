"""Numerical verification checks for a solved state.

Residuals are evaluated on interior bulk nodes only; the interface rows and the
truncation boundaries are checked separately as boundary conditions.
"""
import numpy as np
import solve


def checks(V, a, eps, npanel, ny, L=14.0, verbose=True):
    F = solve.state_fields(V, a, eps, npanel, ny, L)
    X, Y, W = F['X'], F['Y'], F['W']
    WX, WY, S1, S2 = F['WX'], F['WY'], F['S1'], F['S2']
    Dx1, Dy1 = F['Dx1'], F['Dy1']
    K0, C = F['K0'], F['C']
    lam = solve.LAM
    dbar = solve.DELTA_BAR
    nx, nyx = len(X), len(Y)
    out = {}

    # interior bulk nodes: exclude Y=0 interface (j=0) and Y=pi (j=ny), and the
    # two truncation ends (i=0, i=nx-1)
    ii = slice(1, nx - 1)
    jj = slice(1, nyx - 1)

    # 1. bulk PDE Laplace(W)=0 on interior
    Dx2 = solve.diff2(X)
    Dy2 = solve.diff2(Y)
    lap = (Dx2 @ W + W @ Dy2.T)[ii, jj]
    wscale = np.max(np.abs(W[ii, jj]))
    out['laplace_abs'] = float(np.max(np.abs(lap)))
    out['laplace_rel'] = float(np.max(np.abs(lap)) / (wscale * np.pi ** 2))

    # 2. interface conditions
    T2 = np.pi * WY[:, 0] + S2[:, 0]
    ahead = X > 0
    coh = (X < 0) & (X > -eps)
    behind = (X < 0) & (X < -eps - 1.0) & (X > -L + 2.0)
    i_eps = int(np.argmin(np.abs(X + eps)))
    out['W_ahead0_max'] = float(np.max(np.abs(W[ahead, 0])))
    out['T2_cohesive_err'] = float(np.max(np.abs(T2[coh] - 1.0)))
    out['T2_behind_absmax'] = float(np.max(np.abs(T2[behind]))) if behind.any() else 0.0
    out['critical_sep'] = float(2.0 * W[i_eps, 0])
    out['critical_sep_err'] = float(2.0 * W[i_eps, 0] - dbar)

    # 3. outer boundary
    out['top_W_err'] = float(np.max(np.abs(W[:, -1] - a)))
    out['left_W_err'] = float(np.max(np.abs(W[0, :] - a)))
    out['right_W_err'] = float(np.max(np.abs(W[-1, :] - a / np.pi * Y)))

    # 4. constitutive operator residuals on interior
    WXY = (Dx1 @ W @ Dy1.T)[ii, jj]
    res2 = ((Dx1 @ S2 - K0 * S2) - C * WXY)[ii, jj] if False else \
        ((Dx1 @ S2 - K0 * S2) - C * (Dx1 @ W @ Dy1.T))[ii, jj]
    WXX = (Dx1 @ WX)[ii, jj]
    res1 = ((Dx1 @ S1 - K0 * S1) - C * (Dx1 @ WX))[ii, jj]
    out['S2_ode_abs'] = float(np.max(np.abs(res2)))
    out['S1_ode_abs'] = float(np.max(np.abs(res1)))
    out['S2_ode_rel'] = float(np.max(np.abs(res2)) / max(np.max(np.abs(C * WXY)), 1e-30))
    out['S1_ode_rel'] = float(np.max(np.abs(res1)) / max(np.max(np.abs(C * WXX)), 1e-30))

    # trivial constitutive decomposition T_i - E_i = S_i
    out['T_minus_E_err'] = float(max(np.max(np.abs((np.pi * WX + S1) - np.pi * WX - S1)),
                                      np.max(np.abs((np.pi * WY + S2) - np.pi * WY - S2))))
    # S1,S2 match the constitutive integral operator
    out['S2_solved_err'] = float(np.max(np.abs(S2[:, 0] - F['S2solved'])))

    # 5. energy balance
    ratio = solve.bulk_ratio(V, a, eps, npanel, ny, L)
    R = a * a / dbar
    out['bulk_ratio'] = ratio
    out['R'] = R
    out['energy_balance_err'] = ratio - (R - 1.0)

    if verbose:
        for k, v in out.items():
            print('  %-18s %.6e' % (k, v))
    return out


if __name__ == "__main__":
    for (V, a, eps, L) in [(0.001, 0.449402, 0.577000, 14.0),
                           (0.1, 0.515225, 0.768600, 14.0),
                           (10.0, 0.972934, 6.340000, 50.0)]:
        print('==== V=%g a=%.6f eps=%.6f L=%g ====' % (V, a, eps, L))
        checks(V, a, eps, (320, 200, 320), 200, L=L)