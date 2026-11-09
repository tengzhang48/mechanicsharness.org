"""Semi-infinite crack-bridging zone solver (small-scale bridging).

Governing equations for a semi-infinite crack with bridging zone 0 < x < L
(tip at x = 0, zone trailing edge at x = L), closure law p = sigma0 + s*delta,
remote stress intensity K0 applied at infinity:

    delta(x) = (8*K0/E')*sqrt(x/(2*pi))
               - (4/(pi*E')) * int_0^L ln|(sqrt(x)+sqrt(xi))/(sqrt(x)-sqrt(xi))| p(xi) dxi

    K_tip = K0 - sqrt(2/pi) * int_0^L p(xi)/sqrt(xi) dxi

The kernel normalization (4/pi) and the K-kernel sqrt(2/pi) were calibrated
against the independent spectral solver in crack_bridging.py (residuals ~4e-5)
and satisfy the exact Dugdale checks K_tip = -sigma0*sqrt(8L/pi) and
delta_t = 8*sigma0*L/(pi*E') for uniform closure with K_tip = 0.

Unknown: q(u) with delta(L*u) = sqrt(L*u) * q(u), q smooth; collocation at
Chebyshev-sin^2 nodes u_j = sin^2(j*pi/(2N)), j = 1..N.  The logarithmic
singularity at xi = x is removed by subtraction using
I(x) = int_0^L ln|(sqrt(x)+sqrt(xi))/(sqrt(x)-sqrt(xi))| dxi.
"""

import numpy as np
from numpy.polynomial.legendre import leggauss

_GAUSS_CACHE = {}


def _leggauss(m):
    if m not in _GAUSS_CACHE:
        _GAUSS_CACHE[m] = leggauss(m)
    return _GAUSS_CACHE[m]


def I_integral(x, L, M=200):
    """I(x) = int_0^L ln|(sqrt(x)+sqrt(xi))/(sqrt(x)-sqrt(xi))| dxi."""
    if x <= 1e-14:
        return 0.0
    t, w = _leggauss(M)
    total = 0.0
    for x0, x1 in ((0.0, x), (x, L)):
        if x1 - x0 < 1e-300:
            continue
        xi = 0.5 * (x1 - x0) * t + 0.5 * (x1 + x0)
        wi = 0.5 * (x1 - x0) * w
        ker = np.log(np.abs((np.sqrt(x) + np.sqrt(xi)) / (np.sqrt(x) - np.sqrt(xi))))
        ker[~np.isfinite(ker)] = 0.0
        total += np.sum(wi * ker)
    return total


def _q_nodes(N):
    """Nodes u_j = sin^2(j pi/(2N)), j = 1..N on (0,1] (u_N = 1)."""
    j = np.arange(1, N + 1)
    return np.sin(0.5 * j * np.pi / N)**2


def _trap_weights(N):
    """Weights for int_0^1 g du ~ (pi/(2N)) sum g(u_j) sin(j pi/N)."""
    j = np.arange(1, N + 1)
    return (np.pi / (2.0 * N)) * np.sin(j * np.pi / N)


def _lagrange_matrix(nodes, eval_pts):
    """L[i,j] = l_j(eval_pts[i]) for the Lagrange basis through `nodes`."""
    n = len(nodes)
    Lm = np.ones((len(eval_pts), n))
    for j in range(n):
        for k in range(n):
            if k != j:
                Lm[:, j] *= (eval_pts - nodes[k]) / (nodes[j] - nodes[k])
    return Lm


def ssb_zone_solve(L, s=0.0, sigma0=1.0, K0=0.0, free_K0=False,
                   N=48, E=1.0):
    """Solve the zone problem; returns (q, K0, K_tip, info).

    free_K0=False : prescribed K0; N collocation equations for q_1..q_N.
    free_K0=True  : K0 extra unknown, enforced with K_tip = 0 (R-curve state).
    """
    u = _q_nodes(N)
    x = L * u
    Ix = np.array([I_integral(xi, L) for xi in x])

    # quadrature for the subtracted integral: Chebyshev-trapezoid nodes on (0,1)
    Mq = 4 * N
    j = np.arange(1, Mq + 1)
    uq = np.sin(0.5 * j * np.pi / Mq)**2
    wq = L * (np.pi / (2.0 * Mq)) * np.sin(j * np.pi / Mq)
    xq = L * uq
    Lmat = _lagrange_matrix(u, uq)
    delta_mat = np.sqrt(xq)[:, None] * Lmat     # delta(xi_m) = sum_j delta_mat[m,j] q_j

    with np.errstate(divide='ignore', invalid='ignore'):
        kern = np.log(np.abs((np.sqrt(x)[:, None] + np.sqrt(xq)[None, :])
                             / (np.sqrt(x)[:, None] - np.sqrt(xq)[None, :])))
    kern[~np.isfinite(kern)] = 0.0   # coincident nodes: subtracted integrand -> 0
    W = wq[None, :] * kern                       # (N, Mq)
    sumW = W.sum(axis=1)

    pre = 4.0 / (np.pi * E)
    free = bool(free_K0)
    n_unk = N + (1 if free else 0)
    A = np.zeros((n_unk, n_unk))
    rhs = np.zeros(n_unk)

    # row i:  sqrt(x_i) q_i + pre*s*[ I_i*sqrt(x_i) q_i + (W @ delta)_i - sumW_i*sqrt(x_i) q_i ]
    #          - (8 K0 /(E sqrt(2pi))) sqrt(x_i) = - pre*sigma0*I_i
    D = np.sqrt(x)
    A[:N, :N] = np.diag(D) \
        + pre * s * (np.diag(Ix * D) + W @ delta_mat - np.diag(sumW * D))
    rhs[:N] = -pre * sigma0 * Ix
    if free:
        A[:N, N] = -8.0 * D / (E * np.sqrt(2.0 * np.pi))
        # K_tip = 0:  K0 = sqrt(2/pi) * (2 sigma0 sqrt(L) + s * L * sum(wtrap q))
        wtrap = _trap_weights(N)
        A[N, :N] = -np.sqrt(2.0 / np.pi) * s * L * wtrap
        A[N, N] = 1.0
        rhs[N] = np.sqrt(2.0 / np.pi) * 2.0 * sigma0 * np.sqrt(L)
    else:
        rhs = rhs + 8.0 * K0 * D / (E * np.sqrt(2.0 * np.pi))

    u_vec = np.linalg.solve(A, rhs)
    q = u_vec[:N]
    K0_val = u_vec[N] if free else K0
    wtrap = _trap_weights(N)
    K_tip = K0_val - np.sqrt(2.0 / np.pi) * (2.0 * sigma0 * np.sqrt(L)
                                             + s * L * np.sum(wtrap * q))
    delta_t = np.sqrt(L) * q[-1]
    return q, K0_val, K_tip, dict(u=u, delta_t=delta_t)


def ssb_coefficients(beta, L=0.02, N=48, E=1.0):
    """Coefficients c1..c4 of eqs (23)-(26) for the linear spring law with
    beta = sL/E' prescribed (hardening beta>0, softening beta<0).

    c1, c3 from the K_tip = 0 state (K_applied = c1 sigma0 sqrt(L),
    delta_t = c3 sigma0 L/E'); c2, c4 from the pure-K0 subproblem
    (sigma0 = 0, K0 = 1)."""
    s = beta / L
    # K_tip = 0 state
    q0, K0v, Kt, info0 = ssb_zone_solve(L, s=s, sigma0=1.0, free_K0=True, N=N, E=E)
    c1 = K0v / np.sqrt(L)
    c3 = info0['delta_t'] * E / L
    # pure K0 subproblem: closure = s*delta only, K0 = 1
    qB, _, KB, infoB = ssb_zone_solve(L, s=s, sigma0=0.0, K0=1.0, N=N, E=E)
    c2 = 1.0 / KB
    c4 = infoB['delta_t'] * E / (KB * np.sqrt(L))
    return c1, c2, c3, c4


def ssb_rcurve_c1(lam, N=48, E=1.0):
    """Softening triangle law sigma = sigma0*(1 - delta/delta0) (s = -sigma0/delta0,
    so with sigma0=delta0=E'=1, s = -1 and beta = -lam):
    R-curve point with K_tip = 0.  Returns (c1, 2*G/(sigma0 delta0), delta_t/delta0)
    with c1 = K_applied/(sigma0 sqrt(L)), G = c1^2 sigma0^2 L/E'."""
    q, K0, Kt, info = ssb_zone_solve(lam, s=-1.0, sigma0=1.0, free_K0=True, N=N, E=E)
    c1 = K0 / np.sqrt(lam)
    return c1, 2.0 * c1**2 * lam, info['delta_t']
