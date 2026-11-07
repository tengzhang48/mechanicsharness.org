#!/usr/bin/env python3
"""Independent semi-infinite tip-layer collocation for K(0+).

Order of limits: at each Laplace parameter z>0, extract the drained-tip
intensity from r->0, then z->infinity (t->0+). After scaling X=x/sqrt(c/z),
the Laplace-domain inner problem is z-independent.

Crack on X<0: S=0, p=0.
Ligament on X>0: w=0, Q=p_y=0.
Far behind the tip, the outer undrained field has
  w ~ (2(1-nu_u)/G) * K_outer * sqrt(-2X_phys / pi) / something...
We work in units G=1, diffusion length=1, and match
  W ~ sqrt(-X) as X -> -infinity
in MP-37's normalized opening. Then
  K_true / K_outer = n11 / (1+eta)
with eta = (nu_u-nu)/(1-nu_u), 1+eta = (1-nu)/(1-nu_u).

Surface symbol (derived below from Biot plane strain, zero shear),
with k=|xi|, m=sqrt(k^2+1), eta=(nu_u-nu)/(1-nu_u),
b = alpha(1-2nu)/(1-nu), beta = (nu_u-nu)*(1-nu)/(alpha(1-2nu)*(1-nu_u)),
so b*beta = eta:

  S_hat = - (k/(2(1-nu))) * [1 + eta (m-k)/(m+k)] W_hat - b k/(m+k) p_hat
  Q_hat = - m p_hat + (beta k (m-k)/(1-nu)) W_hat

On the crack p=0 so
  S_hat = - (k/(2(1-nu))) D(k) W_hat
with D = 1+eta - 2 eta k^2 /(m(m+k))
which is the drained-stiffness factor: D(infty)=1, D(0)=1+eta (undrained).

The mixed problem is solved by expanding:
  W(X) = H(-X) * sqrt(max(-X,0)) * sum a_j L_j^{(-1/2)} mapped, or simpler:
Chebyshev on a truncated interval.

Method: represent W on (-R,0) as sqrt(-X) * polynomial in a mapped variable,
and p on (0,R) as sqrt(X) * polynomial (pressure is bounded at tip? actually
p ~ sqrt(r) so p=0 at tip; wait MP-37 says p ~ B_p sqrt(R) cos(theta/2),
so on the ligament p ~ B_p sqrt(X). Flux Q on crack ~ X^{-1/2} possibly.

Unknowns:
  W supported on X<0
  p supported on X>0  (p=0 on crack; Q=0 on ligament implies we use
                       p_ligament as unknown OR Q_crack as unknown)

We collocate S=0 on crack and Q=0 on ligament via Fourier quadrature.

Implementation: truncated Fourier integrals with
  W(X) = sqrt(-X) * sum_{j=0}^{N-1} a_j U_j(2(-X)/R - 1)  for -R<X<0, 0 else
  p(X) = sqrt(X)  * sum_{j=0}^{N-1} b_j U_j(2X/R - 1)     for 0<X<R, 0 else
plus a far-field undrained correction: the outer field has constant S=-1
in physical units, but here we enforce homogeneous S=0,p=0 on the crack
and match W ~ C sqrt(-X) at large |X|. For a semi-infinite crack under
no far loading in the inner variables, the driving is the outer square-root
opening, imposed as W + W_part, or equivalently we set the leading
Chebyshev coefficient by matching at the left endpoint.

Simpler driving: impose W(-R) corresponding to unit outer opening
  W_outer = sqrt(-X) * (normalized so coefficient of sqrt is 1),
and solve for the correction that cancels residual traction from the
difference between drained and undrained kernels.

Undrained outer: S = - (k/(2(1-nu_u))) W = - (k (1+eta)/(2(1-nu))) W
Drained inner operator: S = - (k D(k)/(2(1-nu))) W - b k/(m+k) p

If W = sqrt(-X) everywhere on the crack (Fourier known), undrained would
give S=0 only for the infinite-k (local) relation... actually a pure
sqrt(-X) opening of a semi-infinite crack in elasticity produces
uniform S=0 ahead? No: Westergaard, K-field, S ~ 1/sqrt(X) ahead, S=0 on crack.

For elasticity, W = C sqrt(-X) on X<0, W=0 on X>0 is EXACTLY the K-field,
with S=0 on the crack and S = (const) / sqrt(X) on the ligament.

In Fourier: W_hat ~ C Gamma(3/2) (0- i xi)^{-3/2} (one-sided).
S_hat = - mu k W_hat, which is not zero on the crack unless we include
the ligament stress. That's consistent: S is supported on X>0.

For the poroelastic inner problem, the SAME W=C sqrt(-X) does NOT produce
S=0 on the crack because the kernel is not purely |k|. The difference
kernel (D(k)-1) produces a residual traction on the crack which must be
cancelled by a correction opening and a ligament pressure.

We solve:
  W = C sqrt(-X) + W_c   (W_c = o(sqrt(|X|)) as X-> -infty)
  p = p_c on ligament
  S[W,p] = 0 on crack
  Q[W,p] = 0 on ligament

Normalize C so that the outer (undrained) K is 1 in units of sqrt(pi L)
after matching. Matching:
  outer W_phys ~ 4 (1-nu_u) K_outer / (G sqrt(2 pi)) * sqrt(r)
  inner drained K from W ~ 4 (1-nu) K_true / (G sqrt(2 pi)) * sqrt(r) near tip
  If W ~ n11 * sqrt(-X) far back and W ~ n11_tip * sqrt(-X) at tip...

MP-37: far-field W ~ sqrt(-X), inner W ~ n11 sqrt(-X), K_true/sqrt(pi L)=n11/(1+eta)

If far opening coefficient is 1, undrained K_outer corresponds to
K_outer / sqrt(pi L) related by the undrained COD-K relation:
  W ~ [4(1-nu_u)/sqrt(2 pi)] K_outer/G  * sqrt  in physical, scaled.

I'll compute the ratio K_true/K_outer = (1-nu_u)/(1-nu) * (W_tip_coeff / W_far_coeff)
= 1/(1+eta) * (W_near / W_far)

If W_far coeff = 1 and W_near coeff = n11, we get n11/(1+eta).

This script reports both the near-tip COD coefficient and an independent
estimate from ligament stress: K ~ lim sqrt(2 pi X) S, using drained
near-tip elasticity (S ~ K/sqrt(2 pi X)).
"""
import json
import numpy as np
from numpy.polynomial.chebyshev import chebvander

NU, NU_U, ALPHA = 0.2, 0.4, 0.89
ETA = (NU_U - NU) / (1 - NU_U)  # modulus-ratio coupling, independent of alpha
A = ALPHA * (1 - 2 * NU)
B = A / (1 - NU)
BETA = ETA / B  # since b*beta = eta


def kernel_factors(xi):
    k = np.abs(xi)
    m = np.sqrt(k * k + 1.0)
    t = np.divide(k, m + k, out=np.zeros_like(k), where=k > 0)
    D = 1 + ETA - 2 * ETA * k * t / np.maximum(m, 1e-30)
    D = np.where(k == 0, 1 + ETA, D)
    # S_hat = -k/(2(1-nu)) * D * W_hat - B * t * p_hat   (since k/(m+k)=t)
    # Q_hat = -m p_hat + (BETA k (m-k)/(1-nu)) W_hat
    mk = m - k
    return k, m, t, D, mk


def sqrt_opening_ft(xi, eps=1e-12):
    """Fourier transform of H(-X) sqrt(-X), convention int f e^{i xi X} dX.

    H(-x) (-x)^{1/2}  has FT = Gamma(3/2) e^{i 3 pi/4 sgn} / |xi|^{3/2}
    more carefully: int_{-inf}^0 sqrt(-x) e^{i xi x} dx
    Let y=-x: int_0^inf y^{1/2} e^{-i xi y} dy = Gamma(3/2) / (i xi)^{3/2}
    with Re(i xi)>0 i.e. branch of (i xi) with positive real part.
    """
    s = 1j * xi
    # (i xi)^{3/2} with Re > 0
    # i xi = e^{i pi/2} xi for xi>0; i xi = e^{-i pi/2} |xi| for xi<0
    # branch of z^{3/2} with positive real part: use principal sqrt of positive-real branch.
    z = 1j * xi + 0j
    # force Re(z^{1/2})>0 then z^{3/2}=z * z^{1/2}
    # For xi>0, i xi is on +i axis, slightly to the right if we take +eps.
    z = 1j * xi + eps  # shift so Re>0? 1j*xi + eps has Re=eps
    # Wait: the Laplace integral int_0 inf y^{1/2} e^{-p y} dy = Gamma(3/2)/p^{3/2}
    # with p = i xi, Re p >0. For xi real, approach from Re p>0:
    # p = eps + i xi, p^{3/2} with principal value (arg in (-pi/2,pi/2)).
    p = eps + 1j * xi
    return np.math.gamma(1.5) / (p ** 1.5)


def mapped_nodes(n, r):
    """Chebyshev nodes on (0,1) mapped to physical segments."""
    # nodes in (0,1): (1+cos(pi(2j+1)/(2n)))/2 from 1 toward 0? 
    theta = np.pi * (np.arange(n) + 0.5) / n
    u = np.cos(theta)  # (1,-1)
    s = 0.5 * (1 + u)  # (1,0)
    # cluster at both 0 and R: use s in (0,1) as fraction of R, but
    # also cluster at tip: take physical |X| = R * s^2 ? 
    # Use X = R * (1+u)/2 so clustered at endpoints of each panel.
    x_neg = -R * 0.5 * (1 + u)  # from -R to 0, clustered at -R and 0
    x_pos = R * 0.5 * (1 + u)
    return x_neg, x_pos, u


def basis_w(X, coef, R):
    """W = sqrt(-X) * sum c_j T_j(2(-X)/R - 1) on (-R,0)."""
    X = np.asarray(X)
    out = np.zeros_like(X, dtype=float)
    mask = (X < 0) & (X > -R)
    if not np.any(mask):
        return out
    xi = np.clip(2 * (-X[mask]) / R - 1, -1, 1)
    V = chebvander(xi, len(coef) - 1)
    out[mask] = np.sqrt(-X[mask]) * (V @ coef)
    return out


def basis_p(X, coef, R):
    X = np.asarray(X)
    out = np.zeros_like(X, dtype=float)
    mask = (X > 0) & (X < R)
    if not np.any(mask):
        return out
    xi = np.clip(2 * (X[mask]) / R - 1, -1, 1)
    V = chebvander(xi, len(coef) - 1)
    out[mask] = np.sqrt(X[mask]) * (V @ coef)
    return out


def sample_grid(R, n_spatial=4096):
    # uniform grid for FFT, covering [-Lfft, Lfft)
    L = 8 * R
    n = n_spatial
    x = (np.arange(n) - n // 2) * (L / n)
    dx = L / n
    xi = 2 * np.pi * np.fft.fftfreq(n, d=dx)
    return x, xi, dx


def apply_symbol(W, p, xi):
    k, m, t, D, mk = kernel_factors(xi)
    W_hat = np.fft.fft(np.fft.ifftshift(W))
    p_hat = np.fft.fft(np.fft.ifftshift(p))
    # Note: numpy fft is int f e^{-2pi i ...}; we used e^{i xi x} with xi=2pi freq
    # ifftshift/fft convention: if x is ifftshifted (0 at index 0), fft matches
    # int f e^{-i xi x} roughly. Sign of xi only swaps left/right; kernels depend on k=|xi|.
    S_hat = -(k / (2 * (1 - NU))) * D * W_hat - B * t * p_hat
    Q_hat = -m * p_hat + (BETA * k * mk / (1 - NU)) * W_hat
    S = np.fft.fftshift(np.fft.ifft(S_hat).real)
    Q = np.fft.fftshift(np.fft.ifft(Q_hat).real)
    return S, Q


def solve(n_coef=16, R=40.0, n_spatial=8192):
    x, xi, dx = sample_grid(R, n_spatial)
    # collocation points: Chebyshev on crack and ligament, excluding very near
    # the truncation edges.
    theta = np.pi * (np.arange(n_coef) + 0.5) / n_coef
    u = np.cos(theta)
    x_crack = -R * 0.5 * (1 + u)
    x_lig = R * 0.5 * (1 + u)
    # drop points too close to +/-R where truncation is dirty
    keep = np.abs(np.abs(x_crack) - 0)  # keep all; weight later

    # unit bases
    I = np.eye(n_coef)
    S_from_a = np.zeros((n_coef, n_coef))
    Q_from_a = np.zeros((n_coef, n_coef))
    S_from_b = np.zeros((n_coef, n_coef))
    Q_from_b = np.zeros((n_coef, n_coef))

    def interp_at(xq, field, xgrid):
        return np.interp(xq, xgrid, field)

    for j in range(n_coef):
        aj = I[j]
        W = basis_w(x, aj, R)
        p0 = np.zeros_like(x)
        S, Q = apply_symbol(W, p0, xi)
        S_from_a[:, j] = interp_at(x_crack, S, x)
        Q_from_a[:, j] = interp_at(x_lig, Q, x)

        bj = I[j]
        W0 = np.zeros_like(x)
        p = basis_p(x, bj, R)
        S, Q = apply_symbol(W0, p, xi)
        S_from_b[:, j] = interp_at(x_crack, S, x)
        Q_from_b[:, j] = interp_at(x_lig, Q, x)

    # Driving: far-field opening sqrt(-X) truncated, which we already include
    # as the j=0 Chebyshev T0=1 term if we set a0=1 and solve for the rest?
    # Better: RHS from NOT including the far field in unknowns; put
    # W_drive = sqrt(-X) for -R<X<0 (smoothly cut), p=0, compute residual
    # S on crack, Q on ligament, cancel with corrections of the SAME basis
    # including the constant-in-mapped-variable term. Then the total a0
    # is 1 + da0.

    a_drive = np.zeros(n_coef)
    a_drive[0] = 1.0  # T0=1 => W ~ sqrt(-X)
    Wd = basis_w(x, a_drive, R)
    Sd, Qd = apply_symbol(Wd, np.zeros_like(x), xi)
    rhs_S = -interp_at(x_crack, Sd, x)
    rhs_Q = -interp_at(x_lig, Qd, x)

    # We want total W = Wd + Wc, but Wc in same space as Wd would just
    # rescale. Split: freeze far coefficient? 
    # Actually the homogeneous problem admits the outer undrained field
    # only approximately because of truncation. We should NOT include
    # a0 as free if we want W~sqrt at -R. Constrain a_total[0]=1, solve
    # remaining a[1:] and all b.

    # Total operators on (a,b):
    # [S_a S_b][a;b] = 0, but a[0] fixed = 1
    A11 = S_from_a[:, 1:]
    A12 = S_from_b
    A21 = Q_from_a[:, 1:]
    A22 = Q_from_b
    M = np.block([[A11, A12], [A21, A22]])
    rhs = np.concatenate([rhs_S, rhs_Q])
    # residual of drive already is -S(Wd), so
    # S_a[:,1:] a[1:] + S_b b = - S_a[:,0]*1
    # which is rhs as defined if rhs uses Wd = basis e0.
    sol, residuals, rank, sv = np.linalg.lstsq(M, rhs, rcond=None)
    a = np.zeros(n_coef)
    a[0] = 1.0
    a[1:] = sol[: n_coef - 1]
    b = sol[n_coef - 1 :]

    Wtot = basis_w(x, a, R)
    ptot = basis_p(x, b, R)
    S, Q = apply_symbol(Wtot, ptot, xi)

    # Near-tip COD coefficient: W / sqrt(-X) as X->0-
    xtip = x[(x < 0) & (x > -0.05 * min(1.0, R))]
    if len(xtip) > 5:
        wtip = basis_w(xtip, a, R)
        coef_near = np.mean(wtip / np.sqrt(-xtip))
    else:
        coef_near = a.sum()  # T_j(1)=1, X=0- mapped u= -1 wait
    # At tip X=0-, mapped argument 2*0/R-1=-1, T_j(-1)=(-1)^j
    coef_tip = np.sum(a * ((-1.0) ** np.arange(n_coef)))
    # At far X=-R, mapped arg=1, T_j(1)=1, coef_far = sum a_j
    coef_far = np.sum(a)

    # Ligament stress intensity from S * sqrt(2 pi X) as X->0+
    xlig = x[(x > 0) & (x < 0.05)]
    Slim = np.interp(xlig, x, S)
    k_from_s = Slim * np.sqrt(2 * np.pi * xlig)

    # Residual measures
    crack = (x < -1e-3) & (x > -0.8 * R)
    lig = (x > 1e-3) & (x < 0.8 * R)
    resS = float(np.max(np.abs(S[crack])))
    resQ = float(np.max(np.abs(Q[lig])))

    K_ratio_cod = coef_tip / (1 + ETA) / coef_far
    # If far coeff is 1, K_true/sqrt(pi L)_outer-unit = coef_tip/(1+ETA)
    K_ratio_if_far1 = coef_tip / (1 + ETA)

    return {
        "n_coef": n_coef,
        "R": R,
        "n_spatial": n_spatial,
        "eta": ETA,
        "coef_tip": float(coef_tip),
        "coef_far": float(coef_far),
        "K_ratio_cod": float(K_ratio_cod),
        "K_ratio_if_far1": float(K_ratio_if_far1),
        "k_from_s_median": float(np.median(k_from_s)) if len(k_from_s) else None,
        "resS": resS,
        "resQ": resQ,
        "lstsq_residual": float(np.linalg.norm(M @ sol - rhs)),
        "sv_min": float(sv[-1]),
        "sv_max": float(sv[0]),
        "a": a.tolist(),
        "b": b.tolist(),
    }


if __name__ == "__main__":
    results = []
    for n_coef, R, ns in [
        (6, 16.0, 2048),
        (8, 24.0, 4096),
        (12, 32.0, 4096),
    ]:
        row = solve(n_coef=n_coef, R=R, n_spatial=ns)
        results.append(row)
        print(json.dumps({k: v for k, v in row.items() if k not in ("a", "b")}))
    with open("tip_collocation_out.json", "w") as f:
        json.dump(results, f, indent=2)
