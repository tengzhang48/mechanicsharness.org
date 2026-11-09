"""Constraint checks for separable crack-tip and notch-tip fields (topic: crack-tip asymptotics).

Version 0.2 (logs P244, P249): revised after Grok's review and re-check (fixes in the topic's change log).

Any solver's field can be checked here, whatever code produced it. The checks need only the
angular stress functions, the stress exponent and the declared material law. They import no
reference solver.

Field form (2D, small strain, quasi-static, stationary tip, faces traction-free at theta = +-beta;
beta = pi is a crack, beta < pi a notch whose material occupies -beta < theta < beta):
    sigma_ij(r, theta) = r**(-p) * S_ij(theta)          S = (S_rr, S_tt, S_rt)
Material laws (all material constants set to 1; every pass/fail check is homogeneous in them):
    law="power":  eps_ij = 1.5 * Se**(n-1) * s_ij   (J2 deformation theory, pure power law;
                  n = 1 is incompressible linear elasticity)
    law="linear": Hooke's law with Poisson ratio nu (E = 1)
Plane strain uses eps_zz = 0; plane stress uses sigma_zz = 0. Shear strains are TENSOR
components (eps_rt), i.e. half the engineering shear.

Each check reports the conditions under which it applies and its tolerance. Passing shows that
the field is consistent with the declared model. It does not show that the model is adequate
for the problem (see TOPIC.md, "Leaving the scope").

Usage:
    from crack_tip_checks import check_field, format_report
    report = check_field(theta, S_rr, S_tt, S_rt, p=..., n=..., beta=np.pi, plane="strain")
    print(format_report(report))
theta must be increasing and cover [-beta, beta], or [0, beta] for mode I or mode II (the
other half is filled in by symmetry). Use leading=False for a higher-order or non-singular term
(e.g. a T-stress): the checks that fix the LEADING exponent are then not applied.
"""
import numpy as np
from scipy.interpolate import CubicSpline
from scipy.optimize import brentq

__all__ = ["check_field", "format_report", "williams_exponent", "hrr_crack_exponent"]


# ---------------------------------------------------------------- reference values (closed forms)
def hrr_crack_exponent(n):
    """Stress exponent p of the power-law crack-tip field: a finite, nonzero J requires p = 1/(n+1)."""
    return 1.0 / (n + 1.0)


def williams_exponent(beta, mode="I", lam_max=3.0):
    """Leading linear-elastic exponent p = 1 - lambda for a traction-free wedge, faces at +-beta.

    lambda is the smallest real root in (0, lam_max] of
        mode I:  sin(2*lambda*beta) + lambda*sin(2*beta) = 0
        mode II: sin(2*lambda*beta) - lambda*sin(2*beta) = 0, with the rigid-rotation root
                 lambda = 1 divided out; lambda = 1 counts only where it is a double root
                 (tan 2beta = 2beta, total angle 257.45 deg), which gives p = 0
    p > 0 is a singular field; p = 0 (lambda = 1, e.g. beta = pi/2) and p < 0 are bounded.
    Only real roots are searched.
    """
    if abs(beta - np.pi) < 1e-12:
        return 0.5
    sgn = 1.0 if mode == "I" else -1.0
    g0 = lambda lam: np.sin(2 * lam * beta) + sgn * lam * np.sin(2 * beta)
    if mode == "II":
        # lambda = 1 always solves the mode II equation (the rigid rotation). Divide it out:
        # h(lambda) = g(lambda)/(lambda - 1), with h(1) = g'(1) = 2 beta cos 2beta - sin 2beta.
        # A root of h at lambda = 1 is a genuine double root (tan 2beta = 2beta), giving p = 0.
        def g(lam):
            lam = np.asarray(lam, float)
            d = lam - 1.0
            with np.errstate(divide="ignore", invalid="ignore"):
                out = np.where(np.abs(d) > 1e-7, g0(lam) / np.where(d == 0, 1.0, d),
                               2 * beta * np.cos(2 * beta) - np.sin(2 * beta))
            return out if out.ndim else float(out)
    else:
        g = g0
    grid = np.linspace(1e-6, lam_max, 60001)
    v = g(grid)
    roots = []
    for i in range(len(grid) - 1):
        if v[i] == 0.0:
            roots.append(grid[i])
        elif v[i] * v[i + 1] < 0:
            roots.append(brentq(g, grid[i], grid[i + 1], xtol=1e-14))
    if not roots:
        raise ValueError(f"no real root in (0, {lam_max}] for beta = {beta}")
    return 1.0 - min(roots)


# ---------------------------------------------------------------- field construction
def _full_range(theta, S_rr, S_tt, S_rt, mode):
    theta = np.asarray(theta, float)
    S = [np.asarray(a, float) for a in (S_rr, S_tt, S_rt)]
    if theta[0] >= -1e-12:                       # given on [0, beta]: mirror
        if mode not in ("I", "II"):
            raise ValueError("a half range needs mode='I' or 'II'")
        sgn = (1, 1, -1) if mode == "I" else (-1, -1, 1)
        th = np.concatenate([-theta[:0:-1], theta])
        S = [np.concatenate([s_ * a[:0:-1], a]) for s_, a in zip(sgn, S)]
        return th, S
    return theta, S


class _Field:
    def __init__(self, theta, S, p, n, law, plane, nu):
        self.th, self.p, self.n, self.law, self.plane, self.nu = theta, p, n, law, plane, nu
        self.sp = [CubicSpline(theta, s) for s in S]
        self.q = -(n * p if law == "power" else p)   # strain exponent: eps ~ r**q

    def S(self, th):
        return [f(th) for f in self.sp]

    def dS(self, th, k=1):
        return [f(th, k) for f in self.sp]

    def strains(self, srr, stt, srt):
        """Angular strain functions (eps_rr, eps_tt, eps_rt[tensor], eps_zz, Se) from the declared law."""
        if self.law == "power":
            szz = 0.5 * (srr + stt) if self.plane == "strain" else 0.0 * srr   # plane strain: s_zz = 0
            m = (srr + stt + szz) / 3.0
            d = (srr - m, stt - m, srt, szz - m)
            se = np.sqrt(1.5 * (d[0]**2 + d[1]**2 + d[3]**2 + 2 * d[2]**2))
            f = 1.5 * np.where(se > 0, se, 1e-300) ** (self.n - 1)
            return f * d[0], f * d[1], f * d[2], f * d[3], se
        nu = self.nu                               # linear, E = 1
        szz = nu * (srr + stt) if self.plane == "strain" else 0.0 * srr
        tr = srr + stt + szz
        e = lambda s: (1 + nu) * s - nu * tr
        se = np.sqrt(0.5 * ((srr - stt)**2 + (stt - szz)**2 + (szz - srr)**2) + 3 * srt**2)
        return e(srr), e(stt), (1 + nu) * srt, e(szz), se

    def displacement_splines(self):
        """Mode I: u ~ r**(1+q) U(theta) from the strains: U_r = e_rr/(1+q), U_t' = e_tt - U_r, U_t(0) = 0."""
        th = np.linspace(self.th[0], self.th[-1], 8001)
        err, ett, ert, _, _ = self.strains(*self.S(th))
        Ur = err / (1.0 + self.q)
        Ut = CubicSpline(th, ett - Ur).antiderivative()
        Ut_vals = Ut(th) - Ut(0.0)
        return CubicSpline(th, Ur), CubicSpline(th, Ut_vals), (th, err, ett, ert)

    def cart_stress(self, x, y):
        r, th = np.hypot(x, y), np.arctan2(y, x)
        srr, stt, srt = self.S(th)
        c, s = np.cos(th), np.sin(th)
        sxx = srr * c**2 + stt * s**2 - 2 * srt * s * c
        syy = srr * s**2 + stt * c**2 + 2 * srt * s * c
        sxy = (srr - stt) * s * c + srt * (c**2 - s**2)
        return [r**(-self.p) * a for a in (sxx, syy, sxy)]

    def cart_strain(self, x, y):
        r, th = np.hypot(x, y), np.arctan2(y, x)
        err, ett, ert, _, _ = self.strains(*self.S(th))
        c, s = np.cos(th), np.sin(th)
        exx = err * c**2 + ett * s**2 - 2 * ert * s * c
        eyy = err * s**2 + ett * c**2 + 2 * ert * s * c
        exy = (err - ett) * s * c + ert * (c**2 - s**2)
        return [r**self.q * a for a in (exx, eyy, exy)]


# ---------------------------------------------------------------- checks
def _sample_angles(beta, k=41, margin=0.05):
    return np.linspace(-beta + margin, beta - margin, k)


def _res(name, value, tol, applies, conditions, note=""):
    return dict(check=name, value=float(value), tol=tol, applies=applies,
                passed=(bool(value <= tol) if applies else None), conditions=conditions, note=note)


def _na(name, conditions, note=""):
    return dict(check=name, value=float("nan"), tol=None, applies=False, passed=None,
                conditions=conditions, note=note)


def check_field(theta, S_rr, S_tt, S_rt, p, n=1.0, beta=np.pi, plane="strain", law="power",
                nu=0.3, mode="I", leading=True, tol=1e-3, tol_eq=1e-5, reported_J=None,
                h=1e-4, h2=1e-3):
    """Run every applicable check; returns a list of result dicts (see format_report).

    Tolerances assume angular functions resolved to about 1e-6 relative (typically a thousand or
    more points on smooth fields). A marginal failure on a coarse grid means: refine before
    blaming the model. tol_eq (angular equilibrium) is tight so that an exponent error of order
    1e-5 is detected even where no closed-form exponent exists (power-law notches).
    h: step for first differences; h2: step for the second differences of the Saint-Venant check
    (larger, because second differences amplify interpolation error).
    Residuals of the Cartesian checks are divided by the field's overall scale at r = 1, not by
    local term sizes, so correct fields whose terms vanish locally (T-stress, mode II at
    theta = 0) are not rejected.
    """
    th, S = _full_range(theta, S_rr, S_tt, S_rt, mode)
    if abs(th[0] + beta) > 1e-6 or abs(th[-1] - beta) > 1e-6:
        raise ValueError(f"theta must cover [-beta, beta]; got [{th[0]:.6f}, {th[-1]:.6f}]")
    F = _Field(th, S, p, n, law, plane, nu)
    scale = max(np.max(np.abs(a)) for a in S)
    if scale == 0:
        raise ValueError("the angular stress functions are identically zero")
    out = []
    if len(th) < 801:
        out.append(_na("resolution (report only)", "always",
                       f"{len(th)} points over [-beta, beta]; tolerances assume about 1e-6 accuracy; "
                       "refine before blaming the model if a residual check fails marginally"))

    # 1. traction-free faces
    face = max(abs(F.sp[1](a)) + abs(F.sp[2](a)) for a in (-beta, beta)) / scale
    out.append(_res("faces traction-free", face, tol, True,
                    "always (faces at theta = +-beta are traction-free)"))

    # 2. symmetry of the declared mode
    t = np.linspace(0, beta, 401)
    if mode in ("I", "II"):
        e, o = (1, -1) if mode == "I" else (-1, 1)
        par = max(np.max(np.abs(F.sp[0](t) - e * F.sp[0](-t))), np.max(np.abs(F.sp[1](t) - e * F.sp[1](-t))),
                  np.max(np.abs(F.sp[2](t) - o * F.sp[2](-t)))) / scale
        out.append(_res(f"mode {mode} symmetry", par, tol, True, f"declared mode {mode}"))

    # 3. equilibrium, angular form: (1-p)S_rr - S_tt + S_rt' = 0 ; (2-p)S_rt + S_tt' = 0
    tt = np.linspace(-beta, beta, 2001)
    srr, stt, srt = F.S(tt)
    _, dtt, drt = F.dS(tt)
    e1 = (1 - p) * srr - stt + drt
    e2 = (2 - p) * srt + dtt
    eq = max(np.max(np.abs(e1)), np.max(np.abs(e2))) / scale
    out.append(_res("equilibrium (angular ODEs, declared p)", eq, tol_eq, True,
                    "always; uses the declared p and spline derivatives of the given functions; "
                    "tight, so it also tests the exponent"))

    # 4. equilibrium in Cartesian components by finite differences (coordinate-independent)
    res, mags = [], []
    for a in _sample_angles(beta):
        x, y = np.cos(a), np.sin(a)
        dxx = (F.cart_stress(x + h, y)[0] - F.cart_stress(x - h, y)[0]) / (2 * h)
        dxy_y = (F.cart_stress(x, y + h)[2] - F.cart_stress(x, y - h)[2]) / (2 * h)
        dxy_x = (F.cart_stress(x + h, y)[2] - F.cart_stress(x - h, y)[2]) / (2 * h)
        dyy = (F.cart_stress(x, y + h)[1] - F.cart_stress(x, y - h)[1]) / (2 * h)
        res += [abs(dxx + dxy_y), abs(dxy_x + dyy)]
        mags += [abs(dxx) + abs(dxy_y), abs(dxy_x) + abs(dyy)]
    out.append(_res("equilibrium (Cartesian, finite differences)", max(res) / max(max(mags), scale), tol, True,
                    "always; samples away from the faces; residual / field scale at r = 1"))

    # 5. compatibility (Saint-Venant) on strains computed from the declared law
    estrain = max(np.max(np.abs(a)) for a in F.strains(*F.S(tt))[:3])
    res, mags = [], []
    for a in _sample_angles(beta):
        x, y = np.cos(a), np.sin(a)
        E = lambda dx, dy: F.cart_strain(x + dx, y + dy)
        c0 = E(0, 0)
        exx_yy = (E(0, h2)[0] - 2 * c0[0] + E(0, -h2)[0]) / h2**2
        eyy_xx = (E(h2, 0)[1] - 2 * c0[1] + E(-h2, 0)[1]) / h2**2
        exy_xy = (E(h2, h2)[2] - E(h2, -h2)[2] - E(-h2, h2)[2] + E(-h2, -h2)[2]) / (4 * h2**2)
        res.append(abs(exx_yy + eyy_xx - 2 * exy_xy))
        mags.append(abs(exx_yy) + abs(eyy_xx) + 2 * abs(exy_xy))
    out.append(_res("compatibility (Saint-Venant, strains from the declared law)",
                    max(res) / max(max(mags), estrain), 10 * tol, True,
                    f"law={law}, n={n}, plane {plane}; second differences, so tolerance is 10x; "
                    "residual / strain scale at r = 1"))

    # 6. finite energy near the tip: W ~ r**(q - p), so W r dr is integrable iff q - p + 2 > 0
    out.append(dict(check="finite energy near the tip", value=float(F.q - p + 2), tol="> 0",
                    applies=True, passed=bool(F.q - p + 2 > 0),
                    conditions="always; requires p*(n+1) < 2 for the power law (p < 1 linear)",
                    note=f"W*r ~ r^({F.q - p + 1:+.4f})"))

    # 7. the leading exponent: fixed by J for a crack; an eigenvalue for a notch
    crack = abs(beta - np.pi) < 1e-12
    if leading and p <= 0:
        out.append(_na("non-singular leading term (report)", "leading term with p <= 0: stresses stay bounded "
                       "at the tip", "no tip-field dominance; K- or J-type parameters do not apply; failure is "
                       "a stress (strength) question (see TOPIC.md, section 6)"))
    if not leading:
        out.append(_na("leading exponent", "leading=False: a higher-order or non-singular term; "
                       "exponent checks for the leading term not applied"))
    elif crack:
        target = hrr_crack_exponent(n) if law == "power" else 0.5
        out.append(_res("crack exponent fixed by J", abs(p - target), 1e-6, True,
                        "leading term of a crack (beta = pi), flanks parallel to x: finite nonzero J needs p = "
                        + ("1/(n+1)" if law == "power" else "1/2"),
                        "a crack's leading term is non-singular only if the singular amplitude vanishes; "
                        "check a T-stress or higher-order term with leading=False" if p <= 0 else ""))
    else:
        out.append(_na("crack exponent fixed by J", "not a crack: p is an eigenvalue, not fixed by J",
                       "compare instead with limits: beta -> pi, and n = 1 (Williams)"))
        if law == "linear" or n == 1:
            pw = williams_exponent(beta, mode if mode in ("I", "II") else "I")
            out.append(_res(f"linear notch exponent (Williams, mode {mode})", abs(p - pw), 1e-6,
                            mode in ("I", "II"), "leading term, linear (or n = 1), mode I or II",
                            f"Williams p = {pw:.9f}"))

    # 8. displacements from strains, the shear relation, and J (mode I only)
    if mode != "I":
        out.append(_na("displacement and J checks", "implemented for mode I only",
                       "mode II/mixed: use the equilibrium, compatibility and face checks"))
        out.append(_normalization(F, tt))
        return out
    if abs(1.0 + F.q) < 1e-9:
        out.append(_na("displacement and J checks", "not applicable: 1 + q = 0, the displacements are "
                       "logarithmic in r", "this field also fails the finite-energy check"))
        out.append(_normalization(F, tt))
        return out
    Ur, Ut, (tg, err, ett, ert) = F.displacement_splines()
    # admissibility of traction-free faces: the faces must separate. On the upper face the
    # displacement along the outward normal e_theta is U_t(beta) r**(1+q); opening needs U_t(beta) < 0.
    ut_face = float(Ut(beta))
    u_scale = max(float(np.max(np.abs(Ur(tg)))), float(np.max(np.abs(Ut(tg)))), 1e-300)
    if crack:
        out.append(dict(check="crack faces open (no interpenetration)", value=ut_face / u_scale,
                        tol=f"<= {tol} (U_theta(pi)/max|U|)", applies=True,
                        passed=bool(ut_face <= tol * u_scale),
                        conditions="crack, mode I, amplitude sign as given; if the faces would interpenetrate, "
                                   "traction-free faces are not admissible and contact is the missing physics",
                        note=f"U_theta(pi) = {ut_face:.6g}"))
    else:
        out.append(_na("notch faces open (report)", "notch: faces moving toward each other are allowed "
                       "unless they touch", f"U_theta(beta) = {ut_face:.6g} (< 0 means the notch widens)"))
    shear = Ur(tg, 1) + F.q * Ut(tg) - 2 * ert           # 2 e_rt = U_r' + (1+q)U_t - U_t (tensor e_rt)
    shear_res = np.max(np.abs(shear[50:-50])) / max(np.max(np.abs(ert)), estrain, 1e-300)
    out.append(_res("compatibility (displacements from strains, shear relation)", shear_res, tol * 10, True,
                    "mode I; integration constant fixed by symmetry (u_theta(0) = 0)"))

    def u_cart(x, y):
        r, a = np.hypot(x, y), np.arctan2(y, x)
        ur, ut = r**(1 + F.q) * Ur(a), r**(1 + F.q) * Ut(a)
        return ur * np.cos(a) - ut * np.sin(a), ur * np.sin(a) + ut * np.cos(a)

    def J_on_circle(R, k=4001):
        a = np.linspace(-beta, beta, k)
        x, y = R * np.cos(a), R * np.sin(a)
        sxx, syy, sxy = F.cart_stress(x, y)
        exx, eyy, exy = F.cart_strain(x, y)
        # sigma_zz*eps_zz = 0 in plane strain (eps_zz = 0) and in plane stress (sigma_zz = 0)
        W_fac = n / (n + 1.0) if law == "power" else 0.5
        W = W_fac * (sxx * exx + syy * eyy + 2 * sxy * exy)
        dux = (np.array(u_cart(x + h * R, y)) - np.array(u_cart(x - h * R, y))) / (2 * h * R)
        Nx, Ny = np.cos(a), np.sin(a)
        tx, ty = sxx * Nx + sxy * Ny, sxy * Nx + syy * Ny
        return np.trapezoid((W * Nx - (tx * dux[0] + ty * dux[1])) * R, a)

    # J from the field itself, in Cartesian components (no polar-basis algebra). For a separable
    # field J ~ r**(1 + q - p): its r-independence follows from p alone, so contour agreement is
    # not reported as a separate test.
    Js = [J_on_circle(R) for R in (0.5, 1.0, 2.0)]
    out.append(dict(check="J of the field (report)", value=float(Js[1]), tol=None, applies=False, passed=None,
                    conditions=("crack: J = I in the field's normalization, with all material constants = 1"
                                if crack else "not a crack: J on a circle depends on r and is not a crack parameter"),
                    note="J at r = 0.5, 1, 2: " + ", ".join(f"{j:.6g}" for j in Js)))
    if reported_J is not None:
        if crack:
            dev = abs(reported_J - Js[1]) / abs(Js[1])
            out.append(_res("reported J (or I) matches the field's own J", dev, 1e-2, True,
                            "crack, mode I; same normalization and material constants = 1",
                            f"field J = {Js[1]:.6g}, reported = {reported_J:.6g}"))
        else:
            out.append(_na("reported J (or I) matches the field's own J", "not compared: not a crack"))
    out.append(_normalization(F, tt))
    return out


def _normalization(F, tt):
    """Report only: the normalization of the given field (I, K and the field must share it)."""
    se = F.strains(*F.S(tt))[4]
    return dict(check="normalization (report only)", value=float(np.max(se)), tol=None, applies=False,
                passed=None, conditions="state your convention; I, K and the field must share it",
                note=f"max effective stress of the given angular field = {np.max(se):.6g} "
                     f"at theta = {tt[int(np.argmax(se))]:.4f}")


def format_report(report):
    lines = []
    for r in report:
        flag = {True: "PASS", False: "FAIL", None: "n/a "}[r["passed"]]
        v = r["value"]
        lines.append(f"{flag}  {r['check']}: {v:.3g} (tol {r['tol']})" if isinstance(v, float) and v == v
                     else f"{flag}  {r['check']}")
        lines.append(f"        applies: {r['conditions']}" + (f" | {r['note']}" if r.get("note") else ""))
    return "\n".join(lines)
