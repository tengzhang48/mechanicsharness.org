"""
Per-mode (Fourier in x) plane-strain compliance of an isotropic elastic layer,
derived in closed form, plus an independent Chebyshev-collocation BVP solver
used only to validate the closed form.

Conventions
-----------
Layer occupies 0 < y < t. Fields ~ exp(i*xi*x).
Boundary conditions of the mode problem:
    at y = 0 : sigma_xy = q, sigma_yy = p   (mode amplitudes)
    at y = t : traction-free.
Returns complex 2x2 matrix M with
    [U(0), V(0)]^T = M @ [q, p]^T
where U, V are displacement amplitudes (x, y components) at y = 0.

Closed-form solution (lambda = |xi|, sigma = sign(xi), kappa = 3-4nu):
    U(y) = (a + b y) e^{lam (y-t)}    + (a2 + b2 y) e^{-lam y}
    V(y) = i sigma [(-a + kap b/lam - b y) e^{lam (y-t)}
                   + (a2 + kap b2/lam + b2 y) e^{-lam y}]
(obtained from the Navier equations with double roots r = +-lam; the secular
terms are included). Validated against a spectral collocation BVP solver.
"""
import numpy as np


def layer_mode_solution(xi, t, mu, nu, q, p):
    """Closed-form per-mode solution; returns (U, V) as callables of y."""
    lam = abs(xi)
    if lam == 0.0:
        raise ValueError("xi must be nonzero")
    sig = 1.0 if xi > 0 else -1.0
    kap = 3.0 - 4.0 * nu
    Lam = 2.0 * mu * nu / (1.0 - 2.0 * nu)     # plane-strain Lame lambda
    L2 = 2.0 * mu + Lam
    e0 = np.exp(-lam * t)
    f = sig * L2

    A4 = np.zeros((4, 4), dtype=complex)
    # unknown order: a, b, a2, b2
    # row 0: tau(0) = mu (U' + i xi V) at y=0
    A4[0, 0] = 2.0 * lam * e0
    A4[0, 1] = (1.0 - kap) * e0
    A4[0, 2] = -2.0 * lam
    A4[0, 3] = (1.0 - kap)
    # row 1: sig(0) = (2 mu + Lam) V' + i xi Lam U at y=0
    A4[1, 0] = 1j * (f * (-lam) + xi * Lam) * e0
    A4[1, 1] = 1j * f * (kap - 1.0) * e0
    A4[1, 2] = 1j * (f * (-lam) + xi * Lam)
    A4[1, 3] = 1j * f * (1.0 - kap)
    # row 2: tau(t)
    A4[2, 0] = 2.0 * lam
    A4[2, 1] = (1.0 - kap) + 2.0 * lam * t
    A4[2, 2] = -2.0 * lam * e0
    A4[2, 3] = ((1.0 - kap) - 2.0 * lam * t) * e0
    # row 3: sig(t)
    A4[3, 0] = 1j * (f * (-lam) + xi * Lam)
    A4[3, 1] = 1j * (f * (kap - 1.0 - lam * t) + xi * Lam * t)
    A4[3, 2] = 1j * (f * (-lam) + xi * Lam) * e0
    A4[3, 3] = 1j * (f * (1.0 - kap - lam * t) + xi * Lam * t) * e0
    A4[0] *= mu
    A4[2] *= mu
    a, b, a2, b2 = np.linalg.solve(A4, np.array([q, p, 0, 0], complex))

    def U(y):
        return (a + b * y) * np.exp(lam * (y - t)) + (a2 + b2 * y) * np.exp(-lam * y)

    def V(y):
        return 1j * sig * ((-a + kap * b / lam - b * y) * np.exp(lam * (y - t))
                           + (a2 + kap * b2 / lam + b2 * y) * np.exp(-lam * y))

    return U, V


def layer_compliance(xi, t, mu, nu):
    """[U(0), V(0)]^T = M @ [q, p]^T (mode exp(i*xi*x); tractions (q,p) at
    y=0, traction-free at y=t). See module docstring."""
    M = np.zeros((2, 2), dtype=complex)
    for col, (q, p) in enumerate([(1.0, 0.0), (0.0, 1.0)]):
        U, V = layer_mode_solution(xi, t, mu, nu, q, p)
        M[0, col] = U(0.0)
        M[1, col] = V(0.0)
    return M


# ----------------------------------------------------------------------------
# Independent validation solver: Chebyshev collocation of the per-mode BVP
# ----------------------------------------------------------------------------

def cheb_D(N):
    """Chebyshev differentiation matrix (Trefethen), nodes x_j = cos(pi j/N)."""
    x = np.cos(np.pi * np.arange(N + 1) / N)
    c = np.ones(N + 1)
    c[0] *= 2.0
    c[-1] *= 2.0
    c *= (-1.0) ** np.arange(N + 1)
    X = np.tile(x.reshape(-1, 1), (1, N + 1))   # X[i, j] = x_i
    dX = X - X.T                                 # dX[i, j] = x_i - x_j
    D = (c[:, None] / c[None, :]) / (dX + np.eye(N + 1))
    D -= np.diag(np.sum(D, axis=1))
    return x, D


def layer_compliance_bvp(xi, t, mu, nu, N=56):
    """Same per-mode problem solved by Chebyshev collocation (validation only).

    Unknowns: U, V at the N+1 Chebyshev nodes of [0, t]
    (node 0 <-> y = t, node N <-> y = 0).
    Equations: Navier at interior nodes; traction BCs at both ends.
    """
    Lam = 2.0 * mu * nu / (1.0 - 2.0 * nu)
    L2 = 2.0 * mu + Lam
    x, D = cheb_D(N)
    Dy = D * (2.0 / t)
    D2 = Dy @ Dy
    n = N + 1

    def e(j):
        v = np.zeros(n)
        v[j] = 1.0
        return v

    A = np.zeros((2 * n, 2 * n), dtype=complex)
    for j in range(1, N):                       # interior nodes
        r = 2 * j
        # mu(D2-xi^2)U - xi^2 (Lam+mu) U + i xi (Lam+mu) DV = 0
        A[r, :n] = mu * D2[j] - (2.0 * mu + Lam) * xi**2 * e(j)
        A[r, n:] = 1j * xi * (Lam + mu) * Dy[j]
        # mu(D2-xi^2)V + (Lam+mu)(i xi DU + D2 V) = 0
        A[r + 1, :n] = 1j * xi * (Lam + mu) * Dy[j]
        A[r + 1, n:] = (2.0 * mu + Lam) * D2[j] - mu * xi**2 * e(j)
    # BC rows: node N is y=0 (values q, p), node 0 is y=t (values 0, 0)
    for j in (N, 0):
        A[2 * j, :n] = mu * Dy[j]               # tau = mu(DU + i xi V)
        A[2 * j, n:] = 1j * xi * mu * e(j)
        A[2 * j + 1, :n] = 1j * xi * Lam * e(j)  # sig = L2 DV + i xi Lam U
        A[2 * j + 1, n:] = L2 * Dy[j]

    M = np.zeros((2, 2), dtype=complex)
    for col, (q, p) in enumerate([(1.0, 0.0), (0.0, 1.0)]):
        rhs = np.zeros(2 * n, dtype=complex)
        rhs[2 * N] = q
        rhs[2 * N + 1] = p
        sol = np.linalg.solve(A, rhs)
        M[0, col] = sol[N]        # U at y=0
        M[1, col] = sol[n + N]    # V at y=0
    return M
