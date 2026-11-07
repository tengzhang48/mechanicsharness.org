"""Crack-limit solve: 1 - A[w] = 0 on (0,1), w = v/sqrt(1-x), extract K (2 independent ways)."""
import numpy as np
from numpy.polynomial.legendre import leggauss
from numpy.polynomial.chebyshev import chebval
from gbcore import Rfun


def crack_solve(N=48, Mt=400, Mc=None, verbose=True):
    """v(x) = sum a_m T_m(2x-1); w = v/sqrt(1-x); least-squares collocation of A[w]=1."""
    if Mc is None:
        Mc = 3 * N
    # collocation points: Chebyshev-Lobatto interior of order Mc (avoid 0 and 1)
    k = np.arange(1, Mc)
    xc = 0.5 * (1 - np.cos(np.pi * k / Mc))  # in (0,1)
    # quadrature in t
    t, wt = leggauss(Mt)
    t = 0.5 * (t + 1.0)
    wt = 0.5 * wt
    xiq = 1 - t**2  # quadrature nodes in xi
    # basis at quadrature nodes and collocation points
    Tq = np.cos(np.outer(np.arccos(2 * xiq - 1), np.arange(N + 1)))
    Tc = np.cos(np.outer(np.arccos(2 * xc - 1), np.arange(N + 1)))
    vq = Tq  # coeffs applied later
    rows = np.zeros((Mc - 1 + 1, N + 1))  # +1 row for w(0) = -1/4
    rhs = np.ones(Mc - 1 + 1)
    for i, x in enumerate(xc):
        wx = None  # placeholder
        v_x = Tc[i]
        w_x = v_x / np.sqrt(1 - x)  # basis contribution: w(x) = sum a_m T_m(x)/sqrt(1-x)
        tstar = np.sqrt(1 - x)
        num1 = Tq - np.outer(t, w_x)  # (Mt, N+1): v(1-t^2) - t*w(x)
        den1 = x - 1 + t**2
        I1 = 2 * (wt[:, None] * num1 / den1[:, None]).sum(axis=0)
        I2 = w_x * np.log(x / (1 - x))
        R = Rfun(x, xiq)
        I3 = 2 * (wt * R)[:, None] * Tq
        I3 = I3.sum(axis=0)
        rows[i] = I1 + I2 + I3
    # BC row: w(0) = v(0) = -1/4
    rows[-1] = np.cos(np.outer(np.arccos(-1.0), np.arange(N + 1))).ravel()
    rhs[-1] = -0.25
    a, *_ = np.linalg.lstsq(rows, rhs, rcond=None)
    # evaluate
    v_of = lambda xx: chebval(2 * xx - 1, a)
    C = v_of(1.0)
    kappa_C = np.pi * np.sqrt(2) * abs(C)
    # residual at random points
    xs = np.linspace(0.01, 0.99, 25)
    res = []
    for x in xs:
        w_x = v_of(x) / np.sqrt(1 - x)
        Tq2 = np.cos(np.outer(np.arccos(2 * (1 - t**2) - 1), np.arange(N + 1)))
        vq2 = Tq2 @ a
        num1 = vq2 - t * w_x
        den1 = x - 1 + t**2
        I1 = 2 * np.sum(wt * num1 / den1)
        I2 = w_x * np.log(x / (1 - x))
        I3 = 2 * np.sum(wt * Rfun(x, 1 - t**2) * vq2)
        res.append(I1 + I2 + I3 - 1.0)
    # K from traction ahead of tip: s(x) = 1 - A[w](x) ~ kappa/sqrt(2(x-1)) for x>1
    xa = 1 + np.array([0.02, 0.05, 0.1, 0.2])
    sa = []
    for x in xa:
        # direct (no PV since x>1)
        wq = (Tq @ a) / t  # w(1-t^2) = v/t
        I1 = 2 * np.sum(wt * wq * t / (x - (1 - t**2)))  # 2*int t*w/(x-xi) dt... careful
        I3 = 2 * np.sum(wt * Rfun(x, 1 - t**2) * (Tq @ a))
        sa.append(1 - (I1 + I3))
    sa = np.array(sa)
    kappa_fit = np.mean(sa * np.sqrt(2 * (xa - 1)))
    if verbose:
        print(f"crack N={N}: |residual| max {np.max(np.abs(res)):.2e}, "
              f"kappa(C)={kappa_C:.6f}, kappa(fit ahead)={kappa_fit:.6f}, "
              f"handbook edge crack=1.1215, v(0)={v_of(0.0):.6f} (want -0.25)")
    return kappa_C, kappa_fit, v_of


if __name__ == "__main__":
    for N in [24, 32, 48, 64]:
        crack_solve(N=N)
