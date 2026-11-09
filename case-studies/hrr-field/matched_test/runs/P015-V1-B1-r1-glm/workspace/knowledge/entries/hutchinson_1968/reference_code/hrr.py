"""HRR-type near-tip fields for a tensile crack in a power-law hardening material.

Reproduces the main results of
  J. W. Hutchinson, "Singular behaviour at the end of a tensile crack in a
  hardening material", J. Mech. Phys. Solids 16 (1968) 13-31.

Only numpy and scipy are used.

Conventions (all stresses non-dimensionalised by the yield stress, strains by
the yield strain, lengths by the crack half-length):
  * total deformation theory, Ramberg-Osgood hardening  eps = sig + alpha sig^n
  * effective stress
        plane stress : se^2 = sr^2 + st^2 - sr*st + 3*srt^2
        plane strain : se^2 = (3/4)(sr-st)^2 + 3*srt^2
  * Airy stress function  phi = K r^s f(theta),
        sr  = K r^(s-2) (s f + f'')
        st  = K r^(s-2) s(s-1) f
        srt = K r^(s-2) (1-s) f'
  * the angular function f is normalised so that max_theta se_tilde = 1.
"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, least_squares

__all__ = [
    "singular_exponent",
    "effective_stress_coeff",
    "solve_angular",
    "integral_I",
    "energy_coefficient_En",
    "amplitude_K",
    "elastic_angular_stress_function",
    "elastic_angular_displacements",
    "elastic_J_integrand",
    "elastic_displacement_cartesian",
    "elastic_J_plane_stress",
    "elastic_J_numeric",
]


def singular_exponent(n):
    """Dominant singularity exponent s = (2n+1)/(n+1)   (eq. 15 / 26)."""
    return (2.0 * n + 1.0) / (n + 1.0)


# ---------------------------------------------------------------------------
# angular ODE
# ---------------------------------------------------------------------------
def _quadratic_coeff(plane_strain, n, s):
    """Coefficient matrix M with se_tilde^2 = z^T M z, z = (f, f', f'')."""
    M = np.zeros((3, 3))
    if plane_strain:
        M[0, 0] = 0.75 * s * s * (2 - s) ** 2
        M[1, 1] = 3 * (1 - s) ** 2
        M[2, 2] = 0.75
        M[0, 2] = M[2, 0] = 0.75 * s * (2 - s)
    else:
        M[0, 0] = s * s * (s * s - 3 * s + 3)
        M[1, 1] = 3 * (1 - s) ** 2
        M[2, 2] = 1.0
        M[0, 2] = M[2, 0] = s * (3 - s) / 2.0
    return M


def effective_stress_coeff(n, plane_strain):
    """Return ('se^2 = z M z') the 3x3 matrix M for the given n (uses s(n))."""
    return _quadratic_coeff(plane_strain, n, singular_exponent(n))


def _make_rhs(plane_strain, n, s=None):
    """f''''(f,f',f'',f''') obtained by solving the eigenvalue ODE (13)/(31)."""
    if s is None:
        s = singular_exponent(n)
    k = (n - 1.0) / 2.0
    a = n * (s - 2.0)
    b = a + 1.0
    M = _quadratic_coeff(plane_strain, n, s)

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
        # R0 = (d^2 rho/dtheta^2) with the f''''-term removed
        R0 = (k * Q ** (k - 1.0) * (2.0 * vTMv)
              + k * (k - 1.0) * Q ** (k - 2.0) * dQv ** 2
              + g0 * x2 + g1 * x3)
        P = g @ v
        if plane_strain:
            C = s * (2.0 - s) * x0 + x2
            Cp = s * (2.0 - s) * x1 + x3
            coef = rho + g2 * C
            rest = (R0 * C + 2.0 * P * Cp + rho * s * (2.0 - s) * x2
                    - a * (a + 2.0) * rho * C
                    + 4.0 * (s - 1.0) * b * (P * x1 + rho * x2))
        else:
            A = s * (s - 3.0) * x0 - 2.0 * x2
            Ap = s * (s - 3.0) * x1 - 2.0 * x3
            C = s * (2.0 * s - 3.0) * x0 - x2
            coef = 2.0 * rho - g2 * A
            rest = (a * rho * A - R0 * A - 2.0 * P * Ap
                    - rho * s * (s - 3.0) * x2 + a * b * rho * C
                    + 6.0 * b * (s - 1.0) * (P * x1 + rho * x2))
        return -rest / coef

    return rhs


def _integrate(plane_strain, n, p, s=None):
    rhs = _make_rhs(plane_strain, n, s)
    f = lambda t, y: [y[1], y[2], y[3], rhs(y[0], y[1], y[2], y[3])]
    return solve_ivp(f, [0.0, np.pi], [1.0, 0.0, p, 0.0],
                     rtol=1e-11, atol=1e-13, dense_output=True, max_step=0.005)


def _shoot(plane_strain, n, s=None):
    """Find phi''(0) such that phi(pi)=0; among roots pick the one that also
    makes phi'(pi) ~ 0 (the physical dominant mode)."""
    candidates = []
    ps = np.linspace(-2.0, -0.2, 19)
    vals = []
    for p in ps:
        sol = _integrate(plane_strain, n, p, s)
        vals.append(sol.y[0, -1])
    vals = np.asarray(vals)
    for i in range(len(ps) - 1):
        if vals[i] == 0 or vals[i] * vals[i + 1] < 0:
            try:
                root = brentq(lambda p: _integrate(plane_strain, n, p, s).y[0, -1],
                              ps[i], ps[i + 1], xtol=1e-13, rtol=1e-14)
            except Exception:
                continue
            sol = _integrate(plane_strain, n, root, s)
            candidates.append((abs(sol.y[1, -1]), root, sol))
    if not candidates:
        raise RuntimeError("no root bracketed")
    candidates.sort(key=lambda c: c[0])
    return candidates[0][1], candidates[0][2]


# ---------------------------------------------------------------------------
# public solution routine
# ---------------------------------------------------------------------------
def solve_angular(n, plane_strain=True, s=None, normalize=True, n_theta=4001):
    """Solve the dominant angular field.

    Returns a dict with theta grid, f,f',f'',f''', stress components se_tilde,
    sr, st, srt (normalised so max se_tilde = 1), the normalisation constant c
    and the shooting value p = f''(0).
    """
    p, sol = _shoot(plane_strain, n, s)
    th = np.linspace(0.0, np.pi, n_theta)
    Y = sol.sol(th)
    # raw effective stress
    M = _quadratic_coeff(plane_strain, n, s if s is not None else singular_exponent(n))
    Z = Y[:3].T
    Q = np.einsum('ij,jk,ik->i', Z, M, Z)
    se = np.sqrt(Q)
    c = 1.0 / se.max() if normalize else 1.0
    Y = Y * c
    ss = s if s is not None else singular_exponent(n)
    f, fp, fpp, fppp = Y
    sr = ss * f + fpp
    st = ss * (ss - 1.0) * f
    srt = (1.0 - ss) * fp
    Z = Y[:3].T
    se = np.sqrt(np.einsum('ij,jk,ik->i', Z, M, Z))
    return dict(theta=th, f=f, fp=fp, fpp=fpp, fppp=fppp,
                se=se, sr=sr, st=st, srt=srt, c=c, p=p, s=ss)


def _rho(se, n):
    return se ** (n - 1.0)


def integral_I(n, plane_strain=True, sol=None):
    """Path-integral normalisation I, eq. (24) (with the symmetry factor 2)."""
    if sol is None:
        sol = solve_angular(n, plane_strain)
    th = sol["theta"]
    s = sol["s"]
    f, fp, fpp, fppp = sol["f"], sol["fp"], sol["fpp"], sol["fppp"]
    se, sr, st, srt = sol["se"], sol["sr"], sol["st"], sol["srt"]
    M = _quadratic_coeff(plane_strain, n, s)
    m = n * (s - 2.0)
    b = n * (s - 2.0) + 1.0
    rho = _rho(se, n)
    # dQ = 2 M z
    z = np.stack([f, fp, fpp], -1)
    dQ = 2.0 * (z @ M.T)
    dQv = np.einsum('ij,ij->i', dQ, np.stack([fp, fpp, fppp], -1))
    rhop = (n - 1.0) * rho * (dQv / (2.0 * se ** 2))
    # tilde u_r
    if plane_strain:
        c1 = 0.75 * (sr - st)
        c1p = 0.75 * (s * (2.0 - s) * fp + fppp)
    else:
        c1 = sr - 0.5 * st
        c1p = (s * (3.0 - s) / 2.0) * fp + fppp
    ur = rho * c1 / b
    urp = (rhop * c1 + rho * c1p) / b
    # tilde u_theta and its derivative (strain-displacement relations)
    ut = (3.0 * rho * srt - urp) / m
    if plane_strain:
        utp = rho * (0.75 * (st - sr)) - ur
    else:
        utp = rho * (s * (s - 1.5) * f - fpp / 2.0) - ur
    integrand = ((n / (n + 1.0)) * se ** (n + 1.0) * np.cos(th)
                 - (np.sin(th) * (sr * (ut - urp) - srt * (ur + utp))
                    + b * np.cos(th) * (sr * ur + srt * ut)))
    I = 2.0 * np.trapezoid(integrand, th)
    return I


def energy_coefficient_En(n, plane_strain=True, sol=None, I=None):
    """Energy coefficient E_n, eq. (34)."""
    if sol is None:
        sol = solve_angular(n, plane_strain)
    if I is None:
        I = integral_I(n, plane_strain, sol)
    th = sol["theta"]
    shape = 2.0 * np.trapezoid(sol["se"] ** (n + 1.0), th)
    return (n / (n + 1.0)) * (np.pi / I) * shape


def amplitude_K(n, plane_strain=True, alpha=1.0, sigma_inf=1.0, nu=0.3, I=None):
    """Small-scale-yield amplitude K, eq. (27) plane stress / (32) plane strain."""
    if I is None:
        I = integral_I(n, plane_strain)
    amp = 1.0 if not plane_strain else (1.0 - nu ** 2)
    return ((amp / alpha) ** (1.0 / (n + 1.0))
            * (np.pi / I) ** (1.0 / (n + 1.0))
            * sigma_inf ** (2.0 / (n + 1.0)))


def elastic_angular_stress_function(theta):
    """Biharmonic mode-I angular function (eq. 40) with s = 3/2.

    phi = K r^(3/2) (cos(theta/2) + (1/3) cos(3 theta/2)).
    Returns (f, f', f'', f''') at theta.
    """
    t = np.asarray(theta, dtype=float)
    f = np.cos(t / 2.0) + np.cos(1.5 * t) / 3.0
    fp = -0.5 * np.sin(t / 2.0) - 0.5 * np.sin(1.5 * t)
    fpp = -0.25 * np.cos(t / 2.0) - 0.75 * np.cos(1.5 * t)
    fppp = 0.125 * np.sin(t / 2.0) + 1.125 * np.sin(1.5 * t)
    return f, fp, fpp, fppp


def elastic_J_plane_stress(K_I, E):
    """Exact linear-elastic energy release rate, J = K_I^2 / E (plane stress).

    Closed form.  It is cross-checked numerically against `elastic_J_numeric`
    (a first-principles circular J integral) in test_hrr.py.
    """
    return K_I ** 2 / E


def elastic_angular_displacements(theta, nu=0.3):
    """Angular displacement functions of the mode-I elastic field (eq. 40), E=1.

    With u_i = r^(1/2) utilde_i(theta) they follow from the plane-stress
    constitutive law and the strain-displacement relations, so they are
    consistent with the stress components of `elastic_angular_stress_function`.
    Returns (ur, ut, ur', ut') where the prime is d/dtheta.
    """
    t = np.asarray(theta, dtype=float)
    ur = 0.5 * ((5.0 - 3.0 * nu) * np.cos(t / 2.0)
                - (1.0 + nu) * np.cos(1.5 * t))
    ut = 0.5 * ((nu - 7.0) * np.sin(t / 2.0)
                + (1.0 + nu) * np.sin(1.5 * t))
    urp = (-0.25 * (5.0 - 3.0 * nu) * np.sin(t / 2.0)
           + 0.75 * (1.0 + nu) * np.sin(1.5 * t))
    utp = (0.25 * (nu - 7.0) * np.cos(t / 2.0)
           + 0.75 * (1.0 + nu) * np.cos(1.5 * t))
    return ur, ut, urp, utp


def elastic_J_integrand(theta, nu=0.3):
    """Dimensionless integrand of the elastic J integral (E = 1).

    Returns (integrand, terms).  `terms` exposes the angular stress components
    (sr, st, srt), the energy density W, the angular displacements
    (ur, ut, urp, utp) and a = s - 1 = 1/2, i.e. the intermediate quantities
    used by `elastic_J_numeric`, so they can be checked component-wise.
    """
    th = np.asarray(theta, dtype=float)
    f, fp, fpp, _ = elastic_angular_stress_function(th)
    s = 1.5
    sr = s * f + fpp
    st = s * (s - 1.0) * f
    srt = (1.0 - s) * fp
    W = 0.5 * (sr ** 2 + st ** 2 - 2.0 * nu * sr * st
               + 2.0 * (1.0 + nu) * srt ** 2)
    ur, ut, urp, utp = elastic_angular_displacements(th, nu)
    a = s - 1.0
    integrand = (W * np.cos(th)
                 - a * np.cos(th) * (sr * ur + srt * ut)
                 - np.sin(th) * (sr * (ut - urp) - srt * (ur + utp)))
    terms = dict(s=s, sr=sr, st=st, srt=srt, W=W,
                 ur=ur, ut=ut, urp=urp, utp=utp, a=a)
    return integrand, terms


def elastic_displacement_cartesian(x, y, nu=0.3, E=1.0):
    """Cartesian displacement of the unit-amplitude mode-I elastic field (eq. 40).

    The stress amplitude is the angular field used by `elastic_J_numeric`
    (conventionally K_I = sqrt(2*pi)), and u_i = (1/E) r^(1/2) utilde_i(theta).
    Used by the tests to check the displacement-gradient components.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    r = np.hypot(x, y)
    t = np.arctan2(y, x)
    ur, ut, _, _ = elastic_angular_displacements(t, nu)
    ct, stt = x / r, y / r
    scale = r ** 0.5 / E
    return scale * (ur * ct - ut * stt), scale * (ur * stt + ut * ct)


def elastic_J_numeric(K_I, E=1.0, nu=0.3, n_theta=200001):
    """First-principles Rice J integral of the mode-I elastic field (eq. 40).

    Evaluates  J = integral_Gamma (W dy - sigma_ij n_j u_{i,x} ds)  on a circular
    contour of unit radius around the tip, using the printed asymptotic fields
    and the *full* derivative of the displacement vector, including the
    derivatives of the polar basis vectors d(rhat)/dx, d(thetahat)/dx.  For the
    similarity fields u_i = A r^a utilde_i(theta), sigma_ij = B r^(s-2)
    sigma~_ij(theta) this reduces to

      J = (K_I^2/(2*pi*E)) integral_{-pi}^{pi} {
            W~ cos(theta)
            - a cos(theta) (sr~ ur~ + srth~ uth~)
            - sin(theta) [ sr~ (uth~ - ur~') - srth~ (ur~ + uth~') ] } dtheta

    with a = s - 1 = 1/2 and W~ the dimensionless plane-stress energy density.
    The dimensionless integral equals 2*pi for this normalisation (checked in
    test_hrr.py), so the result is J = K_I^2 / E for any K_I, E (E enters only
    through W and the displacement scale).  This is an integral, not the closed
    form returned by `elastic_J_plane_stress`; the two agree to quadrature
    accuracy.
    """
    th = np.linspace(-np.pi, np.pi, n_theta)
    integrand, _ = elastic_J_integrand(th, nu)
    N = np.trapezoid(integrand, th)
    return (K_I ** 2 / (2.0 * np.pi)) * (N / E)
