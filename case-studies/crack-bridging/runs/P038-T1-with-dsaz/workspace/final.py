"""Efficient definitive solve + independent check.

Condition (monotonic branch, first connector failure):
  connector adjacent to original slit, |x|=a, reaches delta = delta_c
  while the crack front, |x|=c=a+L, is at K_tip = K0 = sqrt(E' Gamma0).
Uniform closure sigma0 over a<|x|<c.

Normalized: E'=sigma0=delta_c=1; lengths by ell=E'delta_c/sigma0.
  alpha = a/ell,  K0n = sqrt(Gamma0/(sigma0 delta_c)).
"""

import numpy as np

E = 130e9
nu = 0.28
Ep = E / (1.0 - nu**2)
sigma0 = 160e6
delta_c = 9.00e-6
Gamma0 = 22.0
a_half = 2.40e-3
ell = Ep * delta_c / sigma0
an = a_half / ell
K0n = np.sqrt(Gamma0 / (sigma0 * delta_c))


def solve_unit(c, L, N):
    """Finite crack |x|<c, zone c-L<|x|<c both tips, E'=1.
    Returns (b_closure, K_closure, b_remote, K_remote) where b_* are the
    dislocation coefficients for the two unit loads (closure=1 remote=0, and
    closure=0 remote=1)."""
    theta_a = np.arccos((c - L) / c)
    th = np.pi * np.arange(1, N + 1) / (N + 1)
    n = np.arange(1, N + 1)
    sin_n_th = np.sin(np.outer(th, n))
    basis_t = sin_n_th / np.sin(th)[:, None]
    A = -(1.0 / 4.0) * basis_t
    zone = (th <= theta_a) | (th >= np.pi - theta_a)
    rhs_c = np.where(zone, 1.0, 0.0)
    rhs_r = -np.ones(N)
    b_c = np.linalg.solve(A, rhs_c)
    b_r = np.linalg.solve(A, rhs_r)
    Kc = 0.25 * np.sqrt(np.pi * c) * b_c.sum()
    Kr = 0.25 * np.sqrt(np.pi * c) * b_r.sum()
    return b_c, Kc, b_r, Kr


def state(Ln, N):
    c = an + Ln
    b_c, Kc, b_r, Kr = solve_unit(c, Ln, N)
    Sig = (K0n - Kc) / Kr
    b = b_c + Sig * b_r
    th_a = np.arccos((c - Ln) / c)                 # x = a
    n = np.arange(1, N + 1)
    dt = c * np.sum(b * np.sin(th_a * n) / n)
    return Sig, dt


def bisect(N, iters=20):
    lo, hi = 0.35, 0.75
    for _ in range(iters):
        mid = 0.5 * (lo + hi)
        _, dt = state(mid, N)
        if dt < 1.0:
            lo = mid
        else:
            hi = mid
    Ln = 0.5 * (lo + hi)
    Sig, dt = state(Ln, N)
    return Ln, Sig, dt


if __name__ == '__main__':
    print(f"E'={Ep:.6e}  ell={ell*1e3:.4f} mm  alpha={an:.5f}  K0n={K0n:.5f}")
    for N in [300, 600, 1000]:
        Ln, Sig, dt = bisect(N)
        print(f"N={N:5d}: L={Ln*ell*1e3:.6f} mm  Sigma={Sig*160:.5f} MPa  "
              f"delta(root)/dc={dt:.6f}")