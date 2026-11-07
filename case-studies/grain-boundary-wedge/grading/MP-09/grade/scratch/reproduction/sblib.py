"""Spectral solver in the tip-resolving coordinate s = sqrt(1-x)  (x in [0,1]).

Why: the eigenfunctions of A f = PV int_0^1 S(x,xi) f'''(xi) dxi,
BCs f(0)=0, f'(1)=f''(1)=0, behave like (1-x)^(7/2) near the tip
(i.e. f''' ~ (1-x)^(1/2)), giving only O(1/N) convergence in x-Chebyshev.
In s = sqrt(1-x):  f = sum a_k s^k, the BCs become
    f(s=1) = 0            (sigma(0)=0)
    f_s(0) = f_ss(0) = f_sss(0) = 0  (tip BCs; also kills the ker(A) mode s^4)
so the basis is phi_k(s) = s^k - 1, k = 5,6,7,...  -- polynomial in s,
exactly capturing (1-x)^(k/2) terms including the 7/2 singular one.

Third derivative in s-coordinates (d/dx = -(1/2s) d/ds):
    (s^k)''' = c_k s^(k-6),  c_k = -k(k-2)(k-4)/8   (k>=5; note c_4=0).
Operator:  (A phi_k)(s) = 2 c_k R_{k-5}(s),
    R_m(s) = PV int_0^1 S(1-s^2, 1-r^2) r^m dr          (d xi = -2r dr)
Split: S = 1/(x-xi) + S_reg,  S_reg = -1/(x+xi) - 2 xi (x-xi)/(x+xi)^3.
    R_m = M_m + G_m
    M_m = PV int_0^1 r^m/(r^2 - s^2) dr  (exact recurrence, forward-stable)
    G_m = int_0^1 r^m S_reg(1-s^2, 1-r^2) dr  (smooth in r; Gauss-Legendre)
"""
import numpy as np

def M_moments(s, mmax):
    """M_m = PV int_0^1 r^m/(r^2-s^2) dr, m=0..mmax. Forward recurrence."""
    M = np.zeros(max(mmax, 1)+1)
    if s < 1e-12:
        M[0] = -1.0            # limit of (1/2s) ln((1-s)/(1+s))
    else:
        M[0] = np.log((1.0-s)/(1.0+s))/(2.0*s)
    M[1] = 0.5*np.log((1.0-s*s)/(s*s))
    for m in range(2, mmax+1):
        M[m] = 1.0/(m-1) + s*s*M[m-2]
    return M[:mmax+1]

def S_reg(x, xi):
    """Regular part of the kernel (Cauchy term removed)."""
    return -1.0/(x+xi) - 2.0*xi*(x-xi)/(x+xi)**3

# Gauss-Legendre rule on [0,1] (enough points to resolve the near-mouth
# peak of width ~(a-1) at r=1 for all collocation s)
_NG = 4200
_rg, _wg = np.polynomial.legendre.leggauss(_NG)
_rg = 0.5*(_rg+1.0); _wg = 0.5*_wg

def G_moments(s, mmax):
    """G_m = int_0^1 r^m S_reg(1-s^2, 1-r^2) dr, m=0..mmax (Gauss)."""
    x = 1.0 - s*s
    xi = 1.0 - _rg*_rg
    h = S_reg(x, xi)                        # shape (NG,)
    # moments: G_m = sum_i w_i r_i^m h_i ; build powers iteratively
    G = np.zeros(mmax+1)
    rm = np.ones(_NG)
    for m in range(mmax+1):
        G[m] = np.dot(_wg, rm*h)
        rm = rm*_rg
    return G

def R_moment(s, m):
    M = M_moments(s, m)
    G = G_moments(s, m)
    return M + G

def R_all(s, mmax):
    return M_moments(s, mmax) + G_moments(s, mmax)

def c_coeff(k):
    return -k*(k-2)*(k-4)/8.0

def cheb_points_s(N):
    """N Chebyshev roots points on (0,1) for s, ascending."""
    i = np.arange(1, N+1)
    u = np.cos((2*i-1)*np.pi/(2*N))
    return 0.5*(1.0+u[::-1])

def build_s_matrices(N):
    """Collocation at N Chebyshev s-points; basis phi_k = s^k-1, k=5..N+4.
    Returns s-grid, Theta (A phi_k)(s_i), Fm phi_k(s_i)."""
    ss = cheb_points_s(N)
    ks = np.arange(5, N+5)
    Th = np.zeros((N, N))
    Fm = np.zeros((N, N))
    for i, s in enumerate(ss):
        R = R_all(s, N-1)          # need R_{k-5}, k up to N+4 -> m up to N-1
        Th[i, :] = 2.0*c_coeff(ks)*R[ks-5]
        Fm[i, :] = s**ks - 1.0
    # row scaling: at small s all phi_k ~ -1 -> rows degenerate; rescale rows
    w = np.maximum(np.max(np.abs(Fm), axis=1), 1e-300)
    Th = Th/w[:, None]
    Fm = Fm/w[:, None]
    return ss, Th, Fm

def f_from_coeffs(a, x):
    """f(x) = sum a_j (s^{k_j} - 1), s = sqrt(1-x), k_j = j+5."""
    s = np.sqrt(1.0-x)
    ks = np.arange(5, 5+len(a))
    return np.sum(a*(s**ks - 1.0), axis=-1) if np.ndim(s) else \
        float(np.sum(a*(s**ks - 1.0)))

def eig_spectral_s(N, tol=1e-12):
    """Generalized eigenvalues/vectors of the s-spectral pencil (Th, Fm).

    The raw pencil is numerically singular (cond(Fm) ~ 1e19 at N>=32: the
    columns s^k-1 become indistinguishable in floating point).  Fix:
    truncated-SVD projection -- SVD Fm = U S V^T, keep r singular values
    above tol*sigma_max, and QZ the reduced pencil
        (Ur^T Th Vr,  Ur^T Fm Vr).
    Eigenvectors a = Vr ah are mapped back to the phi_k = s^k-1 basis and
    normalized to f(1) = 1.  Returns (lams, vecs): positive decay rates in
    ascending order, columns of vecs the corresponding coefficient vectors.
    Stability: lam1 changes < 5e-7 for N = 60..110, tol = 1e-11..1e-13.
    """
    from scipy.linalg import eig as geig
    ss, Th, Fm = build_s_matrices(N)
    U, sv, Vt = np.linalg.svd(Fm)
    r = int(np.sum(sv > tol*sv[0]))
    Ur = U[:, :r]
    Vr = Vt[:r, :].conj().T
    Th_r = Ur.T @ Th @ Vr
    Fm_r = Ur.T @ Fm @ Vr
    w, Vh = geig(Th_r, Fm_r, homogeneous_eigvals=True)
    alpha, beta = w
    with np.errstate(divide='ignore', invalid='ignore'):
        lam = alpha/beta
    ok = (np.isfinite(lam)
          & (np.abs(beta) > 1e-9*np.max(np.abs(alpha)))
          & (np.abs(lam.imag) < 1e-7*np.maximum(1.0, np.abs(lam.real)))
          & (lam.real < 0))
    lams = -lam[ok].real
    order = np.argsort(lams)
    lams = lams[order]
    vecs = (Vr @ Vh[:, ok])[:, order]
    for j in range(vecs.shape[1]):
        vecs[:, j] /= np.sum(vecs[:, j])   # f(1) = -sum a  ->  normalize
    return lams, -vecs                     # now f(1) = +1

def mode1_spectral_s(N, tol=1e-12):
    """Slowest mode: returns (lam1, lam2, f_half) with f normalized f(1)=1,
    plus the coefficient vector a1 (phi_k basis) for profile evaluation."""
    lams, vecs = eig_spectral_s(N, tol)
    a1 = vecs[:, 0].real
    f_half = f_from_coeffs(a1, 0.5)
    return lams[0], lams[1], f_half, a1

def offgrid_residual(a1, lam, xs):
    """Residual of the continuous eigen-relation (A f)(x) + lam f(x) = 0 at
    arbitrary off-grid points x (strongest validation of the spectral solve).
    (A f)(x) = sum a_k 2 c_k R_{k-5}(s), s = sqrt(1-x)."""
    ks = np.arange(5, 5+len(a1))
    out = []
    for x in xs:
        s = np.sqrt(1.0-x)
        R = R_all(s, int(ks[-1]-5))
        Af = np.sum(a1*2.0*c_coeff(ks)*R[ks-5])
        out.append(Af + lam*f_from_coeffs(a1, x))
    return np.array(out)

# ================================================================ CRACK / K
def crack_moment_solve(M):
    """Steady-state crack problem  PV int_0^1 S(x,xi) vhat(xi) dxi = 1,
    with vhat(xi) = w(xi)/sqrt(1-xi),  w(xi) = sum_{m=0}^{M-1} c_m (1-xi)^m.
    Collocation at M Chebyshev points in s = sqrt(1-x): since d xi = 2r dr
    and 1-xi = r^2,  (T vhat)(s) = 2 PV int_0^1 S~ r^{2m} dr = 2 R_{2m}(s).
    Returns c (well-conditioned: c = O(1); w(1) = c_0)."""
    ss = cheb_points_s(M)
    A = np.zeros((M, M))
    mm = np.arange(M)
    for i, s in enumerate(ss):
        R = R_all(s, 2*M-2)
        A[i, :] = 2.0*R[2*mm]
    c = np.linalg.solve(A, np.ones(M))
    return c

def M_ext_moments(x, mmax):
    """M_m^ext = int_0^1 r^m/(r^2+q^2) dr with q^2 = x-1 > 0 (outside tip)."""
    q2 = x - 1.0
    M = np.zeros(mmax+1)
    q = np.sqrt(q2)
    M[0] = np.arctan(1.0/q)/q
    if mmax >= 1:
        M[1] = 0.5*np.log((1.0+q2)/q2)
    for m in range(2, mmax+1):
        M[m] = 1.0/(m-1) - q2*M[m-2]
    return M

def T_outside(c, x):
    """sigma_int(x) = int_0^1 S(x,xi) vhat dxi for x > 1 (regular kernel);
    evaluated with exact external moments in the (1-xi)^m basis."""
    Mc = len(c)
    Mex = M_ext_moments(x, 2*Mc-2)
    xi = 1.0 - _rg*_rg
    h = S_reg(x, xi)
    G = np.zeros(2*Mc-1)
    rm = np.ones(_NG)
    for m in range(2*Mc-1):
        G[m] = np.dot(_wg, rm*h)
        rm = rm*_rg
    U = 2.0*(Mex + G)
    return np.sum(c*U[np.arange(Mc)*2])

def crack_K_and_COD(M):
    """Returns dict: c0 (= w(1)), K_tip = -pi sqrt(2) c0, COD mouth
    E*COD(0)/(sigma0 h) = -8 pi sum c_m/(2m+1), and K from the outside
    stress  sigma_bar(x) = 1 - int S vhat,  K = sigma_bar sqrt(2 eps),
    Richardson-extrapolated in sqrt(eps)."""
    c = crack_moment_solve(M)
    c0 = c[0]
    K_tip = -np.pi*np.sqrt(2.0)*c0
    COD = -8.0*np.pi*np.sum(c/(2.0*np.arange(M)+1.0))
    eps = np.array([4e-2, 1e-2, 2.5e-3, 6.25e-4, 1.5625e-4])
    Ke = np.array([np.sqrt(2.0*e)*(1.0 - T_outside(c, 1.0+e)) for e in eps])
    Afit = np.vstack([np.ones_like(eps), np.sqrt(eps), eps]).T
    K_out = np.linalg.lstsq(Afit, Ke, rcond=None)[0][0]
    return {"c0": c0, "K_tip": K_tip, "COD_mouth": COD, "K_out": K_out,
            "c": c}
