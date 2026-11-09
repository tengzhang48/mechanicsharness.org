"""Finite flaw under opening tension with connectors (large-scale bridging).

Model
-----
Plane strain, center slit a<|x| (a = unbridged half-length = 2.40 mm).
Connectors (bridging zone) occupy a < |x| < c = a + L at each end, carrying
a constant closure traction t = sigma0 = 160 MPa (rectilinear law with
delta_c = 9 um).  The intrinsic crack front consumes Gamma0 = 22 J/m^2, i.e.
the front-state tip condition is K_tip = sqrt(E' Gamma0).

Monotonic branch: at first connector failure (outermost connector, |x| = c)
delta(c) = delta_c while all connectors remain < delta_c, so closure = sigma0
everywhere in the zone and K_tip = sqrt(E' Gamma0).  For each L we solve the
crack problem for the remote stress Sigma(L) that satisfies K_tip = K0; the
failing state is the L with delta(c) = delta_c.

Two independent solvers:
  A) spectral dislocation solver for the finite center crack a..c
     (reused/extended from the knowledge unit).
  B) semi-infinite crack in a K-field (front-scale model): valid when L/a,
     a/L are not small (deep LSB), independent construction.
"""
import numpy as np

# ---- material ----
E = 130e9
nu = 0.28
Ep = E / (1.0 - nu**2)          # plane-strain modulus
sigma0 = 160e6
delta_c = 9.00e-6
Gamma0 = 22.0
a_half = 2.40e-3
K0 = np.sqrt(Ep * Gamma0)       # = sqrt(E' Gamma0): tip condition
ell = Ep * delta_c / sigma0     # bridging material length

print(f"E' = {Ep:.6e} Pa")
print(f"K0 = sqrt(E'Gamma0) = {K0:.6e} Pa.sqrt(m)")
print(f"ell = E'delta_c/sigma0 = {ell*1e3:.4f} mm")
print(f"alpha = a/ell = {a_half/ell:.5f}")


# =====================================================================
# Solver A: spectral dislocation solver, finite center crack a<|x|<c,
# uniform closure sigma0 on the two zones a<|x|<c, remote stress Sigma.
# K_tip at the physical front x=c.
# =====================================================================
def theta_grid(m):
    return np.pi * np.arange(1, m + 1) / (m + 1)


def solve_finite(c, L, sigma0_=1.0, remote=0.0, free_remote=False,
                 N=400, E_=1.0):
    """Crack |x|<c, bridging zone c-L<|x|<c at both tips, closure sigma0_,
    remote stress `remote`.  free_remote: Sigma unknown, enforce K_tip=0.
    Returns dict(b, remote, K_tip, delta_at)."""
    theta_a = np.arccos((c - L) / c)
    free = bool(free_remote)
    n_unk = N + (1 if free else 0)
    m_coll = n_unk
    th = theta_grid(m_coll)
    n = np.arange(1, N + 1)
    sin_n_th = np.sin(np.outer(th, n))
    sin_th = np.sin(th)
    with np.errstate(divide='ignore', invalid='ignore'):
        basis_t = sin_n_th / sin_th[:, None]
    basis_d = c * sin_n_th / n
    zone = (th <= theta_a) | (th >= np.pi - theta_a)
    A = np.zeros((m_coll, n_unk))
    A[:, :N] = -(E_ / 4.0) * basis_t
    rhs = np.zeros(m_coll)
    rhs[zone] = sigma0_
    if free:
        A[:, N] = 1.0
        Krow = np.empty(n_unk)
        Krow[:N] = (E_ / 4.0) * np.sqrt(np.pi * c)
        Krow[N] = 0.0
        imid = int(np.argmin(np.abs(th - np.pi / 2.0)))
        A[imid, :] = Krow
        rhs[imid] = 0.0
    else:
        rhs = rhs - remote
    u = np.linalg.solve(A, rhs)
    b = u[:N]
    sig = u[N] if free else remote
    K_tip = (E_ / 4.0) * np.sqrt(np.pi * c) * b.sum()

    def delta_at(theta):
        theta = np.asarray(theta, dtype=float)
        return c * np.sum(b[None, :] * np.sin(np.outer(theta, n)) / n, axis=1)
    return dict(b=b, remote=sig, K_tip=K_tip, delta_at=delta_at, c=c, L=L)


# =====================================================================
# Solver B: semi-infinite crack in a K-field.  Tip at x=0; closure sigma0
# on 0<x<X (zone extends back to the slit at x=-a).  Applied K = Kfront
# (the local field at the front), which relates to the slit via the
# finite-crack geometry.  Here we compute the *exterior* relation: remote
# K_app applied to the front -> opening and K_tip.
# Governing (zone.py convention): delta(x) = (8 K0/E')sqrt(x/2pi)
#   - (4/(pi E')) int_0^X ln|...| p(xi) dxi ;  K_tip = K0 - sqrt(2/pi) *
#   int_0^X p/sqrt(xi) dxi.  With p = sigma0 and K_tip = K0 (front-set):
#   closure-induced K drop = sqrt(2/pi)*2 sigma0 sqrt(X)  (Dugdale).
# =====================================================================
# ---------------------------------------------------------------------
# Quick exploration: finite crack, fixed L, solve Sigma so K_tip = K0;
# inspect delta(x) profile and delta at outer connector.
# ---------------------------------------------------------------------
def sigma_for_ktip(Lfrac, c, Ktgt, N=500):
    """Given c and zone L=Lfrac*c (Lfrac=(c-a)/c), find remote Sigma s.t.
    K_tip = Ktgt (normalized).  Normalize: E'=1, sigma0=1, delta0=1.
    Here we work in *normalized* units separately below."""
    pass

if __name__ == '__main__':
    import sys
    # normalized demo (E'=sigma0=delta0=1): choose a and L, solve for Sigma
    # so that K_tip = K0n, look at delta profile.
    E_ = 1.0
    sig1 = 1.0
    an = a_half / ell
    # normalized tip condition: K_tip/(sigma0 sqrt(ell)) = sqrt(Gamma0/(sigma0 delta_c))
    K0n = np.sqrt(Gamma0 / (sigma0 * delta_c))
    print(f"normalized: alpha={an:.5f} K0n={K0n:.5f}  (sigma0*dc={sigma0*delta_c:.2f} J/m2)")
    for Ln in [0.4, 0.5, 0.6, 0.7, 0.8, 1.2]:
        c = an + Ln
        # linear: K_tip = Ktip(closure only, remote=0) + Sigma*Ktip(unit)
        s0 = solve_finite(c, Ln, sigma0_=sig1, remote=0.0, N=300, E_=E_)
        s1 = solve_finite(c, Ln, sigma0_=0.0, remote=1.0, N=300, E_=E_)
        # want Ktgt: Ktip = Kt_a + Sigma*Kt_b
        Sig = (K0n - s0['K_tip']) / s1['K_tip']
        sol = solve_finite(c, Ln, sigma0_=sig1, remote=Sig, N=300, E_=E_)
        th_a = np.arccos((c - Ln) / c)
        dt = sol['delta_at'](np.array([th_a]))[0]
        # delta at inner connectors
        th_in = np.arccos((c - Ln + 0.02*Ln) / c)
        dt_in = sol['delta_at'](np.array([th_in]))[0]
        print(f"Ln={Ln:.3f} c={c:.4f} Sigma={Sig:.4f} delta_outer={dt:.4f} delta_inner={dt_in:.4f}")
