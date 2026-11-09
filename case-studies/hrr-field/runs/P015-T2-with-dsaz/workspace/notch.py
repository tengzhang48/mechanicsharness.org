"""Power-law (Norton-Hoff / HRR-type) near-tip field for a sharp symmetric V-notch.

Material occupies |theta| < alpha (alpha = 2*pi/3 here).  In-plane polar coordinates,
the notch faces are traction free.  Plane strain, pressure-insensitive,
monotonic proportional loading, small strain.

Stresses scale as sigma ~ r^(s-2) = r^(-p), p = 2 - s.

The angular governing ODE is the plane-strain eigenvalue ODE of Hutchinson (1968),
eq. (31) (same as the one implemented in the knowledge unit's hrr.py), but here the
eigenvalue s is *solved for* from the notch boundary conditions f(alpha)=f'(alpha)=0
instead of being fixed at (2n+1)/(n+1).
"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, fsolve

# --- quadratic form of the effective stress (plane strain), reused from hrr.py ---
def quadratic_coeff(n, s):
    """se_tilde^2 = z^T M z, z = (f, f', f'').  Plane strain, Hutchinson eq. (31)."""
    M = np.zeros((3, 3))
    M[0, 0] = 0.75 * s * s * (2 - s) ** 2
    M[1, 1] = 3 * (1 - s) ** 2
    M[2, 2] = 0.75
    M[0, 2] = M[2, 0] = 0.75 * s * (2 - s)
    return M


def make_rhs(n, s):
    """f''''(f,f',f'',f''') from Hutchinson's plane-strain ODE (eq. 31)."""
    k = (n - 1.0) / 2.0
    a = n * (s - 2.0)
    b = a + 1.0
    M = quadratic_coeff(n, s)

    def rhs(x0, x1, x2, x3):
        z = np.array([x0, x1, x2])
        v = np.array([x1, x2, x3])
        Q = z @ M @ z
        dQ = 2.0 * (M @ z)
        dQv = dQ @ v
        vTMv = v @ M @ v
        rho = Q ** k
        g = k * Q ** (k - 1.0) * dQ
        g0, g1, g2 = g
        R0 = (k * Q ** (k - 1.0) * (2.0 * vTMv)
              + k * (k - 1.0) * Q ** (k - 2.0) * dQv ** 2
              + g0 * x2 + g1 * x3)
        P = g @ v
        C = s * (2.0 - s) * x0 + x2
        Cp = s * (2.0 - s) * x1 + x3
        coef = rho + g2 * C
        rest = (R0 * C + 2.0 * P * Cp + rho * s * (2.0 - s) * x2
                - a * (a + 2.0) * rho * C
                + 4.0 * (s - 1.0) * b * (P * x1 + rho * x2))
        return -rest / coef

    return rhs


def _deriv(rhs, y):
    return np.array([y[1], y[2], y[3], rhs(y[0], y[1], y[2], y[3])])


def integrate_rk4(n, s, p, alpha, nsteps=4000):
    """Fast fixed-step RK4. Returns (f, f') at alpha."""
    rhs = make_rhs(n, s)
    h = alpha / nsteps
    y = np.array([1.0, 0.0, p, 0.0])
    for _ in range(nsteps):
        k1 = _deriv(rhs, y)
        k2 = _deriv(rhs, y + 0.5 * h * k1)
        k3 = _deriv(rhs, y + 0.5 * h * k2)
        k4 = _deriv(rhs, y + h * k3)
        y = y + (h / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
        if not np.all(np.isfinite(y)):
            return np.nan, np.nan
    return y[0], y[1]


def rhs_array(n, s, x0, x1, x2, x3):
    """Vectorized f''''; x0..x3 are arrays (any shape)."""
    k = (n - 1.0) / 2.0
    a = n * (s - 2.0)
    b = a + 1.0
    M = quadratic_coeff(n, s)
    Q = (M[0, 0] * x0 * x0 + M[1, 1] * x1 * x1 + M[2, 2] * x2 * x2
         + 2.0 * M[0, 2] * x0 * x2)
    dQ0 = 2.0 * (M[0, 0] * x0 + M[0, 2] * x2)
    dQ1 = 2.0 * M[1, 1] * x1
    dQ2 = 2.0 * (M[2, 2] * x2 + M[0, 2] * x0)
    dQv = dQ0 * x1 + dQ1 * x2 + dQ2 * x3
    vTMv = (M[0, 0] * x1 * x1 + M[1, 1] * x2 * x2 + M[2, 2] * x3 * x3
            + 2.0 * M[0, 2] * x1 * x3)
    rho = Q ** k
    g0 = k * Q ** (k - 1.0) * dQ0
    g1 = k * Q ** (k - 1.0) * dQ1
    g2 = k * Q ** (k - 1.0) * dQ2
    R0 = (k * Q ** (k - 1.0) * (2.0 * vTMv)
          + k * (k - 1.0) * Q ** (k - 2.0) * dQv ** 2
          + g0 * x2 + g1 * x3)
    P = g0 * x1 + g1 * x2 + g2 * x3
    C = s * (2.0 - s) * x0 + x2
    Cp = s * (2.0 - s) * x1 + x3
    coef = rho + g2 * C
    rest = (R0 * C + 2.0 * P * Cp + rho * s * (2.0 - s) * x2
            - a * (a + 2.0) * rho * C
            + 4.0 * (s - 1.0) * b * (P * x1 + rho * x2))
    return -rest / coef


def integrate_many(n, s, ps, alpha, nsteps=3000):
    """Vectorized RK4 over an array of initial p = f''(0). Returns (f, f') at alpha."""
    ps = np.asarray(ps, dtype=float)
    h = alpha / nsteps
    x0 = np.ones_like(ps)
    x1 = np.zeros_like(ps)
    x2 = ps.copy()
    x3 = np.zeros_like(ps)

    def deriv(x0, x1, x2, x3):
        d = rhs_array(n, s, x0, x1, x2, x3)
        bad = ~np.isfinite(d)
        d = np.where(bad, 0.0, d)
        return x1, x2, x3, d

    for _ in range(nsteps):
        a1 = deriv(x0, x1, x2, x3)
        a2 = deriv(x0 + 0.5 * h * a1[0], x1 + 0.5 * h * a1[1],
                   x2 + 0.5 * h * a1[2], x3 + 0.5 * h * a1[3])
        a3 = deriv(x0 + 0.5 * h * a2[0], x1 + 0.5 * h * a2[1],
                   x2 + 0.5 * h * a2[2], x3 + 0.5 * h * a2[3])
        a4 = deriv(x0 + h * a3[0], x1 + h * a3[1], x2 + h * a3[2], x3 + h * a3[3])
        x0 = x0 + (h / 6.0) * (a1[0] + 2 * a2[0] + 2 * a3[0] + a4[0])
        x1 = x1 + (h / 6.0) * (a1[1] + 2 * a2[1] + 2 * a3[1] + a4[1])
        x2 = x2 + (h / 6.0) * (a1[2] + 2 * a2[2] + 2 * a3[2] + a4[2])
        x3 = x3 + (h / 6.0) * (a1[3] + 2 * a2[3] + 2 * a3[3] + a4[3])
    return x0, x1


def integrate(n, s, p, alpha):
    rhs = make_rhs(n, s)
    f = lambda t, y: [y[1], y[2], y[3], rhs(y[0], y[1], y[2], y[3])]
    sol = solve_ivp(f, [0.0, alpha], [1.0, 0.0, p, 0.0],
                    rtol=1e-12, atol=1e-14, dense_output=True, max_step=1e-3)
    return sol


def residual(n, s, p, alpha):
    sol = integrate(n, s, p, alpha)
    return np.array([sol.y[0, -1], sol.y[1, -1]])  # f(alpha), f'(alpha)


def solve_notch(n, alpha=2.0 * np.pi / 3.0, s0=1.5, p0=-1.0):
    """Solve (s, p) such that f(alpha)=f'(alpha)=0, f(0)=1, f'(0)=f'''(0)=0."""
    def F(x):
        s, p = x
        return residual(n, s, p, alpha)
    # fsolve from a good initial guess
    sol = fsolve(F, [s0, p0], full_output=True, xtol=1e-13)
    xs, info, ier, msg = sol
    return xs[0], xs[1], info


def scan_bg(n, alpha, sgrid, pmin, pmax, Np, outfile):
    import json
    rows = []
    for s in sgrid:
        # find p with f(alpha)=0 via sign changes on a grid
        ps = np.linspace(pmin, pmax, Np)
        vals = []
        for p in ps:
            try:
                vals.append(integrate(n, s, p, alpha).y[0, -1])
            except Exception:
                vals.append(np.nan)
        vals = np.asarray(vals)
        for i in range(Np - 1):
            if np.isfinite(vals[i]) and np.isfinite(vals[i+1]) and vals[i]*vals[i+1] < 0:
                try:
                    r = brentq(lambda p: integrate(n, s, p, alpha).y[0, -1],
                               ps[i], ps[i+1], xtol=1e-12)
                    sol = integrate(n, s, r, alpha)
                    rows.append({"s": float(s), "p": float(r),
                                 "fp": float(sol.y[1, -1])})
                except Exception:
                    pass
    with open(outfile, "w") as fh:
        json.dump(rows, fh)
    return rows
    """For fixed s, return sign changes of f(alpha) over p to locate candidate roots."""
    ps = np.linspace(pmin, pmax, N)
    vals = []
    for p in ps:
        try:
            vals.append(integrate(n, s, p, alpha).y[0, -1])
        except Exception:
            vals.append(np.nan)
    vals = np.asarray(vals)
    roots = []
    for i in range(N - 1):
        if np.isfinite(vals[i]) and np.isfinite(vals[i + 1]) and vals[i] * vals[i + 1] < 0:
            try:
                r = brentq(lambda p: integrate(n, s, p, alpha).y[0, -1],
                           ps[i], ps[i + 1], xtol=1e-13)
                roots.append(r)
            except Exception:
                pass
    return roots


if __name__ == "__main__":
    for n in [4.0]:
        s, p, info = solve_notch(n)
        print(f"n={n} notch alpha=2pi/3: s={s:.12f}  p_sing=2-s={2-s:.12f}  p_fpp0={p:.10f}")