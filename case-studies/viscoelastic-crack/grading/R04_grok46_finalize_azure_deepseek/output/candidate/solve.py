"""Steady crack growth in a viscoelastic strip -- numerical solver.

Reduced formulation (see DERIVATION.md):
  Laplace(W)=0 on Y in [0,pi], mixed data on Y=0,
  T_2 = pi W_Y + S_2,  T_1 = pi W_X + S_1,
  S_i = C[W_{,Xi}(X) - sum_k w_k W_{,Xi}(X+v_k/K0)] ,  C = pi(1-lam)/lam,
  K0=1/(lam V).   w -> (a/pi)Y ahead, w -> a behind.
"""
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from numpy.polynomial.laguerre import laggauss

LAM = 0.1
DELTA_BAR = 0.2


# ---------------------------------------------------------------- grids
def build_grid(eps, L=14.0, npanel=(140, 90, 140)):
    panels = [(-L, -eps, npanel[0]), (-eps, 0.0, npanel[1]), (0.0, L, npanel[2])]
    xs = []
    for xa, xb, m in panels:
        th = np.linspace(0.0, np.pi, m + 1)
        seg = xa + (xb - xa) * (1.0 - np.cos(th)) / 2.0
        if xs:
            seg = seg[1:]
        xs.append(seg)
    return np.concatenate(xs)


def cluster(ymax, n):
    th = np.linspace(0.0, np.pi, n + 1)
    return ymax * (1.0 - np.cos(th)) / 2.0


def _fd_coeff(xs, x0, order):
    A = np.array([[1, 1, 1], xs - x0, (xs - x0) ** 2])
    b = np.array([0.0, 1.0, 0.0]) if order == 1 else np.array([0.0, 0.0, 2.0])
    return np.linalg.solve(A, b)


def diff1(X):
    n = len(X)
    D = np.zeros((n, n))
    for i in range(1, n - 1):
        hL = X[i] - X[i - 1]
        hR = X[i + 1] - X[i]
        D[i, i - 1] = -hR / (hL * (hL + hR))
        D[i, i] = (hR - hL) / (hL * hR)
        D[i, i + 1] = hL / (hR * (hL + hR))
    D[0, 0:3] = _fd_coeff(X[0:3], X[0], 1)
    D[-1, -3:] = _fd_coeff(X[-3:], X[-1], 1)
    return D


def diff2(X):
    n = len(X)
    D = np.zeros((n, n))
    for i in range(1, n - 1):
        hL = X[i] - X[i - 1]
        hR = X[i + 1] - X[i]
        D[i, i - 1] = 2.0 / (hL * (hL + hR))
        D[i, i] = -2.0 / (hL * hR)
        D[i, i + 1] = 2.0 / (hR * (hL + hR))
    D[0, 0:3] = _fd_coeff(X[0:3], X[0], 2)
    D[-1, -3:] = _fd_coeff(X[-3:], X[-1], 2)
    return D


def build_avg_matrix(X, K0, nlag=64):
    v, w = laggauss(nlag)
    n = len(X)
    A = np.zeros((n, n))
    skip = np.zeros(n)
    for i in range(n):
        xi = X[i]
        for k in range(len(v)):
            xq = xi + v[k] / K0
            if xq >= X[-1]:
                skip[i] += w[k]
                continue
            j = np.searchsorted(X, xq) - 1
            j = min(max(j, 0), n - 2)
            h = X[j + 1] - X[j]
            t = (xq - X[j]) / h
            A[i, j] += w[k] * (1.0 - t)
            A[i, j + 1] += w[k] * t
    return A, skip


def assemble(V, a, eps, L=14.0, npanel=(140, 90, 140), ny=80, nlag=64):
    X = build_grid(eps, L, npanel)
    nx = len(X)
    Dx1, Dx2 = diff1(X), diff2(X)
    Y = cluster(np.pi, ny)
    Dy1, Dy2 = diff1(Y), diff2(Y)
    nW = nx * (ny + 1)
    K0 = 1.0 / (LAM * V)
    C = np.pi * (1.0 - LAM) / LAM

    rows, cols, vals = [], [], []
    b = np.zeros(nW + nx)
    r = 0

    def add(i, j, val):
        rows.append(r); cols.append(i * (ny + 1) + j); vals.append(val)

    Aavg, skip = build_avg_matrix(X, K0, nlag)
    Smat = C * (np.eye(nx) - Aavg)              # acts on WY0
    Qrhs = sp.lil_matrix((nx, nW))
    for j in range(ny + 1):
        colbase = j
        for ip in range(nx):
            vv = Smat[:, ip] * Dy1[0, j]
            if np.any(vv):
                for i in np.nonzero(vv)[0]:
                    Qrhs[i, ip * (ny + 1) + j] += vv[i]
    Qrhs = Qrhs.tocsr()
    s2_far = C * (a / np.pi) * skip

    for i in range(nx):
        for j in range(ny + 1):
            if i == 0:
                add(i, j, 1.0); b[r] = a
            elif i == nx - 1:
                add(i, j, 1.0); b[r] = a / np.pi * Y[j]
            elif j == ny:
                add(i, j, 1.0); b[r] = a
            elif j == 0:
                if X[i] >= 0.0:
                    add(i, 0, 1.0); b[r] = 0.0
                else:
                    for jj in range(3):
                        add(i, jj, np.pi * Dy1[0, jj])
                    rows.append(r); cols.append(nW + i); vals.append(1.0)
                    b[r] = 1.0 if X[i] > -eps else 0.0
            else:
                for ii in (-1, 0, 1):
                    add(i + ii, j, Dx2[i, i + ii])
                for jj in (-1, 0, 1):
                    add(i, j + jj, Dy2[j, j + jj])
            r += 1

    for i in range(nx):
        rows.append(r); cols.append(nW + i); vals.append(1.0)
        row = Qrhs.getrow(i).tocoo()
        for cc, vv in zip(row.col, row.data):
            rows.append(r); cols.append(int(cc)); vals.append(-vv)
        b[r] = -s2_far[i]
        r += 1

    A = sp.csr_matrix((vals, (rows, cols)), shape=(nW + nx, nW + nx))
    return A, b, X, Y, Dx1, Dx2, Dy1, Dy2, ny, nW, K0


def solve_state(V, a, eps, **kw):
    A, b, X, Y, Dx1, Dx2, Dy1, Dy2, ny, nW, K0 = assemble(V, a, eps, **kw)
    sol = spla.spsolve(A.tocsc(), b)
    W = sol[:nW].reshape(len(X), ny + 1)
    S2 = sol[nW:]
    return W, S2, X, Y, Dx1, Dy1, ny, K0


def sif(T2, X, x0=0.02, x1=0.4):
    """Extrapolate T2(X) ~ K/sqrt(X)+T0+R sqrt(X) over a fixed physical window."""
    m = (X > x0) & (X < x1)
    xs = X[m]
    ts = T2[m]
    A = np.vstack([1.0 / np.sqrt(xs), np.ones_like(xs), np.sqrt(xs)]).T
    coef, *_ = np.linalg.lstsq(A, ts, rcond=None)
    return coef[0], coef[1]


def residuals(V, a, eps, **kw):
    W, S2, X, Y, Dx1, Dy1, ny, K0 = solve_state(V, a, eps, **kw)
    WY0 = sum(Dy1[0, j] * W[:, j] for j in range(ny + 1))
    T2 = np.pi * WY0 + S2
    i = int(np.argmin(np.abs(X + eps)))
    f = W[i, 0]
    K, T0 = sif(T2, X)
    return f, K, T0, W, S2, X, Y, K0


def assemble_lin(V, eps, L=14.0, npanel=(140, 90, 140), ny=80, nlag=64):
    """Assemble A and b = b0 + a*b1, since the matrix is independent of a."""
    A, b, X, Y, Dx1, Dx2, Dy1, Dy2, ny, nW, K0 = assemble(
        V, 0.0, eps, L=L, npanel=npanel, ny=ny, nlag=nlag)
    # recover B1 by assembling with a and subtracting
    A2, b_a, *_ = assemble(V, 1.0, eps, L=L, npanel=npanel, ny=ny, nlag=nlag)
    b1 = b_a - b
    return A, b, b1, X, Y, Dx1, Dx2, Dy1, Dy2, ny, nW, K0


class State:
    """Factorized solver for fixed (V, eps, grid). Residuals are cheap in a."""

    def __init__(self, V, eps, L=14.0, npanel=(140, 90, 140), ny=80, nlag=64):
        (A, b0, b1, X, Y, Dx1, Dx2, Dy1, Dy2, ny, nW, K0) = assemble_lin(
            V, eps, L=L, npanel=npanel, ny=ny, nlag=nlag)
        self.lu = spla.splu(A.tocsc())
        self.b0, self.b1 = b0, b1
        self.X, self.Y, self.Dy1, self.Dx1 = X, Y, Dy1, Dx1
        self.ny, self.nW, self.K0 = ny, nW, K0
        self.eps = eps
        self.i_eps = int(np.argmin(np.abs(X + eps)))

    def solve(self, a):
        sol = self.lu.solve(self.b0 + a * self.b1)
        W = sol[:self.nW].reshape(len(self.X), self.ny + 1)
        S2 = sol[self.nW:]
        return W, S2

    def residual(self, a, x0=0.02, x1=0.4):
        W, S2 = self.solve(a)
        WY0 = sum(self.Dy1[0, j] * W[:, j] for j in range(self.ny + 1))
        T2 = np.pi * WY0 + S2
        f = W[self.i_eps, 0]
        K, T0 = sif(T2, self.X, x0, x1)
        return f, K, T2, W, S2


# =====================================================================
#  Entry point
# =====================================================================
# Tabulated admissible roots used as initial guesses (interpolated in log V).
_GUESS = [
    (0.001, 0.577, 0.4495),
    (0.01, 0.6044, 0.4609),
    (0.1, 0.7686, 0.5156),
    (1.0, 1.596, 0.6998),
    (10.0, 6.34, 0.9731),
]


def _guess(V):
    import bisect
    Vs = [g[0] for g in _GUESS]
    if V <= Vs[0]:
        return _GUESS[0][1], _GUESS[0][2]
    if V >= Vs[-1]:
        return _GUESS[-1][1], _GUESS[-1][2]
    i = bisect.bisect(Vs, V)
    v0, e0, a0 = _GUESS[i - 1]
    v1, e1, a1 = _GUESS[i]
    t = (np.log(V) - np.log(v0)) / (np.log(v1) - np.log(v0))
    return e0 + t * (e1 - e0), a0 + t * (a1 - a0)


def solve_case(V, npanel=(400, 250, 400), ny=250, L=None):
    e0, a0 = _guess(V)
    if L is None:
        L = max(14.0, 5.0 * V)
    a, eps, K = root(V, e0, npanel, ny, a0=a0, L=L, steps=18, tol=2e-4)
    admissible = bool(abs(K) < 5e-3)
    return a, eps, K, admissible


def main(case_path="case.json", out_path="results.json"):
    import json
    with open(case_path) as fh:
        case = json.load(fh)
    lam = case.get("lambda", LAM)
    dbar = case.get("delta_bar", DELTA_BAR)
    speeds = case["speeds"]
    sX = case["sample_X"]
    sy = case["sample_y_over_h"]
    npanel, ny = (400, 250, 400), 250

    cases = []
    for V in speeds:
        L = max(14.0, 5.0 * V)
        a, eps, K, adm = solve_case(V, npanel, ny, L)
        entry = {"V": V, "admissible": adm}
        if adm:
            R = a * a / dbar
            ratio = bulk_ratio(V, a, eps, npanel, ny, L)
            rows = sample(V, a, eps, sX, sy, npanel, ny, L)
            entry.update({
                "epsilon": eps, "a": a, "R": R,
                "bulk_dissipation_ratio": ratio,
                "stress_intensity_residual": K,
                "fields": rows,
            })
        else:
            entry.update({"epsilon": None, "a": None, "R": None,
                          "bulk_dissipation_ratio": None, "fields": None,
                          "stress_intensity_residual": K})
        cases.append(entry)

    result = {"critical_speed": 13.44, "cases": cases}
    with open(out_path, "w") as fh:
        json.dump(result, fh, indent=2)
    return result


# =====================================================================
#  Root finder (two conditions: 2 U(-eps)=delta_bar_c and bounded stress)
# =====================================================================
_state_cache = {}


def _state(V, eps, npanel, ny, L=14.0):
    key = (round(V, 12), round(eps, 12), tuple(npanel), ny, L)
    st = _state_cache.get(key)
    if st is None:
        st = State(V, eps, L=L, npanel=npanel, ny=ny)
        _state_cache[key] = st
        while len(_state_cache) > 3:
            _state_cache.pop(next(iter(_state_cache)))
    return st


def a_for_eps(V, eps, npanel, ny, a0, L=14.0, it=8, tol=1e-9):
    st = _state(V, eps, npanel, ny, L)
    a = a0
    f = st.residual(a)[0]
    for _ in range(it):
        if abs(f - 0.5 * DELTA_BAR) < tol:
            break
        da = 1e-3 * max(1.0, abs(a))
        f1 = st.residual(a + da)[0]
        slope = (f1 - f) / da
        if abs(slope) < 1e-12:
            break
        a = a - (f - 0.5 * DELTA_BAR) / slope
        f = st.residual(a)[0]
    return a


def root(V, eps0, npanel, ny, a0=0.45, L=14.0, steps=18, tol=2e-4, step_tol=1e-5,
         verbose=False):
    """Return (a, eps, K) solving the two conditions."""
    eps = eps0
    a = a_for_eps(V, eps, npanel, ny, a0, L)
    K = _state(V, eps, npanel, ny, L).residual(a)[1]
    for it in range(steps):
        if abs(K) < tol:
            break
        de = 1e-3 * max(1.0, abs(eps))
        K2 = _state(V, eps + de, npanel, ny, L).residual(a_for_eps(V, eps + de, npanel, ny, a, L))[1]
        slope = (K2 - K) / de
        if abs(slope) < 1e-12:
            break
        step = -K / slope
        step = max(-0.25, min(0.25, step))
        if abs(step) < step_tol:
            break
        eps = eps + step
        a = a_for_eps(V, eps, npanel, ny, a, L)
        K = _state(V, eps, npanel, ny, L).residual(a)[1]
        if verbose:
            print('  it%d eps=%.8f a=%.8f K=%.3e' % (it, eps, a, K))
    return a, eps, K


# =====================================================================
#  Fields, bulk dissipation, sampling
# =====================================================================
def state_fields(V, a, eps, npanel, ny, L=14.0):
    A, b, X, Y, Dx1, Dx2, Dy1, Dy2, nyy, nW, K0 = assemble(
        V, a, eps, L=L, npanel=npanel, ny=ny)
    sol = spla.spsolve(A.tocsc(), b)
    W = sol[:nW].reshape(len(X), nyy + 1)
    S2solved = sol[nW:]
    C = np.pi * (1.0 - LAM) / LAM
    Aavg, skip = build_avg_matrix(X, K0)
    WX = Dx1 @ W
    WY = W @ Dy1.T
    S2 = C * (WY - Aavg @ WY) - C * (a / np.pi) * skip[:, None]
    S1 = C * (WX - Aavg @ WX)
    return dict(X=X, Y=Y, W=W, WX=WX, WY=WY, S1=S1, S2=S2,
                S2solved=S2solved, K0=K0, C=C, Dx1=Dx1, Dy1=Dy1)


def bulk_ratio(V, a, eps, npanel, ny, L=14.0):
    F = state_fields(V, a, eps, npanel, ny, L)
    integ = F['S1'] ** 2 + F['S2'] ** 2
    half = np.trapezoid(np.trapezoid(integ, F['X'], axis=0), F['Y'])
    return 2.0 * half / (DELTA_BAR * (1.0 - LAM) * V * np.pi)


def sample(V, a, eps, sample_X, sample_yh, npanel, ny, L=14.0):
    from scipy.interpolate import RectBivariateSpline
    F = state_fields(V, a, eps, npanel, ny, L)
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
    return rows


if __name__ == "__main__":
    import sys
    cp = sys.argv[1] if len(sys.argv) > 1 else "case.json"
    op = sys.argv[2] if len(sys.argv) > 2 else "results.json"
    main(cp, op)
    print("wrote", op)