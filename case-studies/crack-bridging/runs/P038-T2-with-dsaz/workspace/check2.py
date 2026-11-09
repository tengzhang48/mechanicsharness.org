"""Fast independent checks using the closed-form particular solution.

The linear softening ODE d'''' = -(1-d) (normalized) has particular solution 1
and homogeneous e^{+-z}, cos z, sin z.  We confirm the analytic solve.py result
with a fully independent collocation construction.
"""
import numpy as np
import solve as S

dc, b, sb, Dtil, Lc = S.dc, S.b, S.sb, S.Dtil, S.Lc

def solve_norm(Zl, Mhat):
    """Return d(z) on [0,Zl]; independent construction (complex exponentials)."""
    # d = 1 + A e^{z} + B e^{-z} + C cos z + D sin z, d'(...) etc.
    r = np.array([1.0, -1.0, 1j, -1j])
    def basis(z, k):
        return r**k * np.exp(r * z)
    rows = []
    rhs = []
    # z=0: d=0, d'=0
    rows.append(list(basis(0, 0))); rhs.append(-1.0)
    rows.append(list(basis(0, 1))); rhs.append(0.0)
    # z=Zl: d=1, d'=free?  No -- use d''(Zl)=Mhat, d'''(Zl)=0 AND d(Zl)=1?
    # We have 4 unknowns; use: d(0)=0, d'(0)=0, d''(Zl)=Mhat, d'''(Zl)=0.
    rows.append(list(basis(Zl, 2))); rhs.append(Mhat)
    rows.append(list(basis(Zl, 3))); rhs.append(0.0)
    Amat = np.array(rows, dtype=complex)
    coef = np.linalg.solve(Amat, np.array(rhs, dtype=complex))

    def d(z):
        z = np.atleast_1d(np.asarray(z, float))
        return 1.0 + sum(c * np.exp(rr * z) for c, rr in zip(coef, r))

    return d, coef

if __name__ == "__main__":
    M = np.sqrt(2 * S.b * S.Dtil * (S.Gam0 + sb * dc / 2))
    Mhat = M * Lc**2 / (Dtil * dc)
    # find Zl such that d(Zl)=1
    from scipy.optimize import brentq
    def f(Zl):
        d, _ = solve_norm(Zl, Mhat)
        return float(np.real(d(Zl)[0])) - 1.0
    Zl = brentq(f, 0.5, 4.0, xtol=1e-13)
    d, _ = solve_norm(Zl, Mhat)
    print("M (N m)      =", M)
    print("L (mm)       =", Zl * Lc * 1e3)
    print("delta(L) um  =", float(np.real(d(Zl)[0])) * dc * 1e6)
    print("delta(L/2) um=", float(np.real(d(Zl/2)[0])) * dc * 1e6)
    # compare to solve.py
    delta, _ = S.solve_delta(Zl * Lc, M)
    print("solve.py delta(L/2) um =", delta(Zl * Lc / 2)[0] * 1e6)