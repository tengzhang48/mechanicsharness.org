"""Reproduction of the main quantitative results of

    Bao G., Suo Z. (1992). Remarks on Crack-Bridging Concepts.
    Applied Mechanics Reviews 45, 355-366.

Method
------
Small-scale and large-scale bridging problems for a center crack |x| < c in an
infinite elastic plate (plane stress/strain lumped into E') are solved with a
spectral integral-equation method:

  * crack opening  delta(c cos(phi)) = c * sum_n b_n sin(n phi) / n
  * (climb) dislocation density  sin(phi) * b(c cos(phi)) = sum_n b_n cos(n phi)
    (Appendix eq. A1 of the paper: b = -d delta/dx)
  * traction induced on the crack line:
        t(theta) = -(E'/4) * sum_n b_n sin(n theta) / sin(theta)
  * crack-tip stress intensity:
        K_tip = (E'/4) sqrt(pi c) * sum_n b_n

The kernel identities are exact (finite-part/Hadamard integrals of Chebyshev
polynomials); the unknown coefficients b_1..b_N are found by collocation at
theta_k = k*pi/(M+1).  Collocation at M = N points gives a square system for a
prescribed remote stress; the R-curve / strength states (K_tip = 0, free remote
stress) use N collocation points plus the extra equation K_tip = 0.

Verification of the kernel: uniform crack-face pressure p over the whole crack
gives delta(x) = (4p/E') sqrt(c^2-x^2) and K = p sqrt(pi c).

Units are normalized throughout: E' = sigma0 = delta0 = c = 1 unless noted.
"""

import numpy as np


# ---------------------------------------------------------------------------
# core solver
# ---------------------------------------------------------------------------

def _theta_grid(m):
    return np.pi * np.arange(1, m + 1) / (m + 1)


def solve_bridged(L, s=0.0, sigma0=1.0, remote=0.0, free_remote=False,
                  N=200, E=1.0, c=1.0):
    """Solve the symmetrically bridged center crack.

    Geometry: crack |x| < c; bridging zone c-L < |x| < c at both tips;
    traction-free faces for |x| < c-L; remote stress `remote` (if not free).

    Closure (bridging) law on the zone: p(delta) = sigma0 + s*delta(x),
    where delta is the *total* crack opening.  t_el(delta) = -remote - p.

    free_remote=False : `remote` is prescribed; unknowns b_1..b_N;
                        N collocation equations.
    free_remote=True  : remote stress is an extra unknown enforced together
                        with K_tip = 0 (R-curve / strength states);
                        N collocation equations + 1 K_tip equation.

    Returns dict(b, remote, K_tip, delta_at(theta) function).
    """
    theta_a = np.arccos((c - L) / c)          # zone = theta < theta_a
    n_mode = n_unk = N
    free = bool(free_remote)
    if free:
        n_unk += 1
    m_coll = n_unk                           # square system

    th = _theta_grid(m_coll)
    n = np.arange(1, N + 1)
    sin_n_th = np.sin(np.outer(th, n))       # (m_coll, N)
    sin_th = np.sin(th)
    with np.errstate(divide='ignore', invalid='ignore'):
        basis_t = sin_n_th / sin_th[:, None]         # sin(n th)/sin(th)
    basis_d = c * sin_n_th / n                      # delta(theta)
    # symmetric bridging zones at BOTH tips: theta < theta_a or theta > pi - theta_a
    zone = (th <= theta_a) | (th >= np.pi - theta_a)

    A = np.zeros((m_coll, n_unk))
    # traction operator  -(E/4) * sum b_n sin(n th)/sin(th)
    A[:, :N] = -(E / 4.0) * basis_t
    # Bridging pulls the faces together: sigma_yy^p = +p_b - remote on the zone,
    # with p_b = sigma0 + s*delta.  Move s*delta and remote to the left.
    A[zone, :N] -= s * basis_d[zone, :]
    rhs = np.zeros(m_coll)
    rhs[zone] = sigma0
    if free:
        A[:, N] = 1.0                        # + remote on the left side
        Krow = np.empty(n_unk)
        Krow[:N] = (E / 4.0) * np.sqrt(np.pi * c)
        Krow[N] = 0.0
        # replace the collocation row nearest mid-crack (least information)
        imid = int(np.argmin(np.abs(th - np.pi / 2.0)))
        A[imid, :] = Krow
        rhs[imid] = 0.0                      # K_tip = 0
    else:
        rhs = rhs - remote                   # prescribed remote moves to rhs

    u = np.linalg.solve(A, rhs)
    b = u[:N]
    sig = u[N] if free else remote

    K_tip = (E / 4.0) * np.sqrt(np.pi * c) * b.sum()

    def delta_at(theta):
        theta = np.asarray(theta, dtype=float)
        return c * np.sum(b[None, :] * np.sin(np.outer(theta, n)) / n, axis=1)

    return dict(b=b, remote=sig, K_tip=K_tip, delta_at=delta_at,
                E=E, c=c, L=L, s=s, sigma0=sigma0)


def trailing_edge_theta(L, c=1.0):
    return np.arccos((c - L) / c)


# ---------------------------------------------------------------------------
# Table 2: hardening bridging, coefficients c1, c2, c3, c4
# ---------------------------------------------------------------------------

def hardening_coefficients(beta, L=0.01, N=200):
    """Coefficients c1..c4 of eqs (23)-(26) for the linear spring law
    sigma = sigma0 + s*delta with beta = sL/E' prescribed (delegates to the
    semi-infinite zone solver; L and N kept for signature compatibility)."""
    from zone import ssb_coefficients
    return ssb_coefficients(beta, L=0.02, N=128)


TABLE2_8BETA_PI = np.array([0.25, 0.5, 1.0, 2.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0])
TABLE2_C2 = np.array([1.244, 1.478, 1.91, 2.68, 4.44, 6.54, 8.15, 9.51, 10.70, 11.76])
TABLE2_C1 = np.array([1.670, 1.737, 1.836, 1.984, 2.184, 2.307, 2.356, 2.386, 2.404, 2.414])


def reproduce_table2(L=0.01, N=200):
    """Reproduce Table 2 (hardening bridging): c1 from the solver, c2 from the
    paper's eq. (26), c2_tab = tabulated values."""
    beta = TABLE2_8BETA_PI * np.pi / 8.0
    out = []
    for btab, c2t, c1t in zip(beta, TABLE2_C2, TABLE2_C1):
        c1, c2, c3, c4 = hardening_coefficients(btab, L=L, N=N)
        out.append(dict(beta=btab, c1=c1, c2_from_eq26=c2,
                        c2_paper=c2t, c1_paper=c1t, c3=c3, c4=c4))
    return out


# ---------------------------------------------------------------------------
# Table 3: softening bridging (linear softening triangle, s = -sigma0/delta0)
# ---------------------------------------------------------------------------

TABLE3_4LAMBDA_PI = np.array([0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.466])
TABLE3_2G = np.array([0.192, 0.366, 0.521, 0.657, 0.771, 0.864, 0.933,
                      0.978, 0.998, 1.000])
TABLE3_C1 = np.array([1.564, 1.526, 1.487, 1.446, 1.401, 1.354, 1.303,
                      1.248, 1.189, 1.169])


def softening_rcurve_point(lam, N=200):
    """R-curve point of the softening triangle law sigma = sigma0*(1-delta/delta0)
    at lam = L*sigma0/(E'*delta0)  (= -beta since s = -sigma0/delta0).

    State: crack growing with K_tip = 0; returns
        c1 = K_applied/(sigma0*sqrt(L)),  two_G = 2*G/(sigma0*delta0) = 2*c1^2*lam,
        delta_t/delta0.
    (Delegates to the semi-infinite zone solver; N kept for compatibility.)"""
    from zone import ssb_rcurve_c1
    return ssb_rcurve_c1(lam, N=48)


def reproduce_table3(N=200):
    out = []
    for lam4pi, g2t, c1t in zip(TABLE3_4LAMBDA_PI, TABLE3_2G, TABLE3_C1):
        lam = lam4pi * np.pi / 4.0
        c1, g2, dt = softening_rcurve_point(lam, N=N)
        out.append(dict(lam=lam, c1=c1, two_G=g2, delta_t=dt,
                        c1_paper=c1t, two_G_paper=g2t))
    return out


def softening_plateau(N=200, tol=1e-4):
    """Steady state of the softening triangle: lam* with delta_t = delta0
    (eq. 28: lam* = 0.366) and plateau energy (eq. 27: G = sigma0*delta0/2)."""
    # delta_t increases monotonically with lam; bracket the root
    lo, hi = 0.30, 0.50
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        _, _, dt = softening_rcurve_point(mid, N=N)
        if dt < 1.0:
            lo = mid
        else:
            hi = mid
        if hi - lo < tol:
            break
    lam = 0.5 * (lo + hi)
    c1, g2, dt = softening_rcurve_point(lam, N=N)
    return lam, c1, 0.5 * g2          # lam*, c1(lam*), G_plateau/(sigma0 delta0)


# ---------------------------------------------------------------------------
# Rectilinear law closed forms (eqs 10, 11, 18, 19, 20) and Fig. 8
# ---------------------------------------------------------------------------

def c1_rectilinear(L=0.01, N=200):
    """eq. (19): K_applied = (8/pi)^{1/2} sigma0 sqrt(L) + K_tip, i.e.
    c1 = (8/pi)^{1/2} = 1.59577 (semi-infinite zone solver, exact)."""
    from zone import ssb_zone_solve
    q, K0, Kt, info = ssb_zone_solve(L if L else 0.02, s=0.0, sigma0=1.0,
                                     K0=0.0, N=48)
    return -Kt / np.sqrt(L if L else 0.02)


def L0_rectilinear(Gamma0_over_sigma0delta0=0.0):
    """eq. (20), units of E'*delta0/sigma0:
    L0 = (pi/8) [(1 + g)^{1/2} + g^{1/2}]^{-2},  g = Gamma0/(sigma0 delta0)."""
    g = Gamma0_over_sigma0delta0
    return (np.pi / 8.0) * (np.sqrt(1.0 + g) + np.sqrt(g))**(-2)


def steady_state_rectilinear(Gamma0_over_sigma0delta0, N=48):
    """Independent numerical route to eq. (20) -- does NOT use the closed-form
    Dugdale coefficients that derive eq. (20).

    For a trial zone size l (units of E'*delta0/sigma0) the semi-infinite zone
    integral equation is solved numerically for the two unit loads
    (closure only / remote K only); by linearity the steady-state conditions

        K_tip = sqrt(E'*Gamma0)          (tip at criticality, K0^2/E' = Gamma0)
        delta_t = delta0

    give a 2x2 linear system for the remote K at that l.  The zone size is
    found by bisection on delta_t - delta0.  A wrong bracket exponent in
    eq. (20) is rejected by this route (a +2 bracket is off by a factor
    [(sqrt(1+g)+sqrt(g))/(sqrt(1+g)-sqrt(g))]^2)."""
    from zone import ssb_zone_solve
    g = Gamma0_over_sigma0delta0
    K_req = np.sqrt(g)                    # K_tip required, units sigma0*sqrt(E'delta0/sigma0)

    def state(l):
        # basis (i): closure only; basis (ii): remote K only
        _, _, Kt_a, info_a = ssb_zone_solve(l, s=0.0, sigma0=1.0, K0=0.0, N=N)
        _, _, Kt_b, info_b = ssb_zone_solve(l, s=0.0, sigma0=0.0, K0=1.0, N=N)
        K0 = (K_req - Kt_a) / Kt_b        # linearity: K_tip = Kt_a + K0*Kt_b
        dt = info_a['delta_t'] + K0 * info_b['delta_t']
        return dt

    lo, hi = 1e-6, 100.0
    while state(hi) < 1.0:
        hi *= 10.0
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if state(mid) < 1.0:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)                # L0 in units E' delta0 / sigma0


# ---------------------------------------------------------------------------
# Strength, Section 3.2: eqs (44) and (47), Fig. 10
# ---------------------------------------------------------------------------

def strength_SSB(alpha):
    """eq. (44): Griffith/small-scale bridging limit."""
    return 1.0 / np.sqrt(np.pi * alpha)


def strength_LSB(alpha):
    """eq. (47): large-scale bridging, rectilinear law, K_tip = 0, delta_t = delta0.
    Equivalent closed form: (2/pi) arccos(exp(-pi/(8 alpha)))."""
    return (2.0 / np.pi) * np.arccos(np.exp(-np.pi / (8.0 * alpha)))


from functools import lru_cache


@lru_cache(maxsize=None)
def strength_numerical(alpha, s=0.0, N=700, nL=30):
    """Strength from the bridged-crack solution (infinite panel, unbridged
    half-length a = alpha, total c = a + L): for each damage extent L solve
    the K_tip = 0 state (unit sigma0, delta0, E'); failure when delta_t =
    delta0.  The strength is the applied stress at the limiting state; when
    sigma(L) is not monotonic in L (softening law) the paper's search
    applies: maximize sigma(L) subject to delta_t(L) <= delta0.

    s = 0    : rectilinear law  (compare with eq. 47)
    s = -1.0 : softening triangle sigma = sigma0*(1 - delta/delta0) (Fig. 11)
    """
    eta_max = 1.5 * np.exp(np.pi / (8.0 * alpha)) + 5.0
    etas = np.geomspace(1.0 + 1e-4, eta_max, nL)
    states = []
    for eta in etas:
        L = alpha * (eta - 1.0)
        c = alpha + L
        sol = solve_bridged(L, s=s, sigma0=1.0, free_remote=True, N=N, c=c)
        th_a = trailing_edge_theta(L, c=c)
        dt = sol['delta_at'](np.array([th_a]))[0]
        if dt > 0.0:
            states.append((eta, sol['remote'], dt))

    def state(eta):
        L = alpha * (eta - 1.0)
        c = alpha + L
        sol = solve_bridged(L, s=s, sigma0=1.0, free_remote=True, N=N, c=c)
        th_a = trailing_edge_theta(L, c=c)
        dt = sol['delta_at'](np.array([th_a]))[0]
        return sol['remote'], dt

    feas = [(e, sg, dt) for e, sg, dt in states if dt <= 1.0]
    infeas = [e for e, sg, dt in states if dt > 1.0]
    if not feas:
        return np.nan
    refined = list(feas)
    if infeas:
        lo = max(e for e, sg, dt in feas)
        hi_candidates = [e for e in infeas if e > lo]
        if hi_candidates:
            hi = min(hi_candidates)
            for _ in range(8):
                mid = 0.5 * (lo + hi)
                sg, dt = state(mid)
                if dt <= 1.0:
                    lo = mid
                    refined.append((lo, sg, dt))
                else:
                    hi = mid
        # local refinement around the sigma maximum (interior max for
        # softening laws)
        es = np.array([f[0] for f in feas])
        sgs = np.array([f[1] for f in feas])
        imax = int(np.argmax(sgs))
        el = np.geomspace(es[max(imax - 1, 0)], es[min(imax + 1, len(es) - 1)], 7)
        for e in el:
            if lo * 0.8 <= e <= hi * 1.2:
                sg, dt = state(e)
                if dt <= 1.0:
                    refined.append((e, sg, dt))
    return max(sg for e, sg, dt in refined)


# ---------------------------------------------------------------------------
# Notch, Section 3.3: eqs (48), (49), Fig. 13 limiting curves
# ---------------------------------------------------------------------------

def notch_brittle_limit(aw):
    """eq. (48) with the limiting stress concentration factor.
    Only the value C(0) = 3 is fixed by the paper (a/w -> 0, Fig. 13)."""
    if np.any(np.asarray(aw) != 0.0):
        raise ValueError("C(a/w) is not given analytically in the paper; "
                         "only the a/w -> 0 limit C = 3 is reproducible.")
    return 1.0 / 3.0


def notch_ductile_limit(aw):
    """eq. (49): net-section yield, sigma_bar_max/sigma_y = 1 - a/w."""
    return 1.0 - np.asarray(aw, dtype=float)


# ---------------------------------------------------------------------------
# Table 1
# ---------------------------------------------------------------------------

def table1():
    """Table 1 (p. 356): sigma0*delta0 [J/m^2] and delta0*E/sigma0 [m]
    for E = 1e11 N/m^2.  Returns the raw inputs and the products AS PRINTED
    in the paper.  Verified against the page image (300 dpi): the atomic-bond
    row is sigma0 = 1e10 N/m^2, delta0 = 1e-10 m, hence sigma0*delta0 = 1 J/m^2
    and delta0*E/sigma0 = 1e-9 m -- the table is internally consistent and
    there is NO erratum (an earlier draft misread delta0 as 1e-9 and invented
    a spurious 10x energy erratum)."""
    E = 1e11
    mech = {'atomic bond': (1e10, 1e-10),
            'fiber pull-out': (1e9, 1e-5),
            'fiber cross-over': (1e7, 1e-4)}
    printed = {'atomic bond': (1.0, 1e-9),
               'fiber pull-out': (1e4, 1e-3),
               'fiber cross-over': (1e3, 1.0)}
    out = {}
    for name, (s0, d0) in mech.items():
        pg, pl = printed[name]
        out[name] = dict(sigma0=s0, delta0=d0,
                         sigma0_delta0_computed=s0 * d0,
                         sigma0_delta0_printed=pg,
                         delta0_E_over_sigma0=d0 * E / s0,
                         delta0_E_over_sigma0_printed=pl)
    return out


# ---------------------------------------------------------------------------
# Fig. 8 R-curves
# ---------------------------------------------------------------------------

def rcurve_rectilinear_x(x):
    """Fig. 8, rectilinear law (same Gamma, sigma0 as itself):
    G/Gamma = L/L0 with L0 = pi/8 (units E' delta0/sigma0); kink at 0.393."""
    return np.minimum(x / (np.pi / 8.0), 1.0)


def rcurve_softening_points(n=25, N=200):
    """Fig. 8, softening triangle with delta0 doubled so that Gamma matches
    the rectilinear curve (paper's normalization).  x = L/(E' delta0^r/sigma0)
    with delta0^r the rectilinear limiting separation: x = 2*lam.
    Kink at x = 2*0.366 = 0.732, y = 2*c1^2*lam."""
    lam = np.linspace(1e-4, softening_plateau(N=N)[0], n)
    y = np.array([softening_rcurve_point(l, N=N)[1] for l in lam])
    return 2.0 * lam, y
