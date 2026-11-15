"""Near-tip field of an interface crack, derived numerically (no formulas copied from tables).

Geometry: material 1 occupies 0 < theta < pi (x2 > 0), material 2 occupies -pi < theta < 0,
crack faces at theta = +-pi, bonded interface at theta = 0.  Plane strain.

A complexified Airy function Phi = r^(lam+1) [a cos((lam+1)t) + b sin((lam+1)t)
                                              + c cos((lam-1)t) + d sin((lam-1)t)]
is taken in each material (8 complex unknowns).  Because the field equations, the constitutive
law and the boundary/interface conditions are real-linear, the real and imaginary parts of a
complex solution are both real solutions.  The 8 homogeneous conditions (traction-free faces,
traction and displacement continuity across theta = 0) give M(lam) v = 0; lam is a root of
det M near 1/2 + i*eps.  The physical field for a complex SIF Ka is Re(Ka * Fe), with Fe the
eigenfield normalised so that on theta = 0:  sigma_tt = r^(lam-1)/sqrt(2 pi).

The convention check: with that normalisation, sigma_rt(theta=0) must equal -i r^(lam-1)/sqrt(2pi)
for lam = 1/2 + i eps; then  sigma22 + i sigma12 = Ka r^(i eps)/sqrt(2 pi r)  on the interface,
exactly the task's definition.
"""
import numpy as np


def mu_kappa(E, nu):
    return E / (2 * (1 + nu)), 3 - 4 * nu


def bimat_eps(E1, nu1, E2, nu2):
    mu1, k1 = mu_kappa(E1, nu1)
    mu2, k2 = mu_kappa(E2, nu2)
    return np.log((k1 / mu1 + 1 / mu2) / (k2 / mu2 + 1 / mu1)) / (2 * np.pi)


def _rows_stress(lam, t):
    """Coefficients (for unknowns a,b,c,d) of r^(1-lam)*[s_rr, s_tt, s_rt] at angle t."""
    p, m = lam + 1, lam - 1
    cp, sp, cm, sm = np.cos(p * t), np.sin(p * t), np.cos(m * t), np.sin(m * t)
    srr = np.array([-lam * p * cp, -lam * p * sp, lam * (3 - lam) * cm, lam * (3 - lam) * sm])
    stt = np.array([lam * p * cp, lam * p * sp, lam * p * cm, lam * p * sm])
    srt = np.array([lam * p * sp, -lam * p * cp, lam * m * sm, -lam * m * cm])
    return srr, stt, srt


def _rows_disp(lam, t, mu, kap):
    """Coefficients of r^(-lam)*[u_r, u_t] at angle t."""
    p, m = lam + 1, lam - 1
    cp, sp, cm, sm = np.cos(p * t), np.sin(p * t), np.cos(m * t), np.sin(m * t)
    ur = np.array([-p * cp, -p * sp, (kap - lam) * cm, (kap - lam) * sm]) / (2 * mu)
    ut = np.array([p * sp, -p * cp, (kap + lam) * sm, -(kap + lam) * cm]) / (2 * mu)
    return ur, ut


def _rows_ddisp(lam, t, mu, kap):
    """Coefficients of r^(-lam)*[d u_r/dt, d u_t/dt]."""
    p, m = lam + 1, lam - 1
    cp, sp, cm, sm = np.cos(p * t), np.sin(p * t), np.cos(m * t), np.sin(m * t)
    dur = np.array([p * p * sp, -p * p * cp, -(kap - lam) * m * sm, (kap - lam) * m * cm]) / (2 * mu)
    dut = np.array([p * p * cp, p * p * sp, (kap + lam) * m * cm, (kap + lam) * m * sm]) / (2 * mu)
    return dur, dut


def eig_matrix(lam, mu1, k1, mu2, k2):
    M = np.zeros((8, 8), dtype=complex)
    _, stt, srt = _rows_stress(lam, np.pi)
    M[0, 0:4] = stt
    M[1, 0:4] = srt
    _, stt, srt = _rows_stress(lam, -np.pi)
    M[2, 4:8] = stt
    M[3, 4:8] = srt
    _, stt1, srt1 = _rows_stress(lam, 0.0)
    M[4, 0:4], M[4, 4:8] = stt1, -stt1
    M[5, 0:4], M[5, 4:8] = srt1, -srt1
    ur1, ut1 = _rows_disp(lam, 0.0, mu1, k1)
    ur2, ut2 = _rows_disp(lam, 0.0, mu2, k2)
    M[6, 0:4], M[6, 4:8] = ur1, -ur2
    M[7, 0:4], M[7, 4:8] = ut1, -ut2
    return M


def find_root(lam0, mu1, k1, mu2, k2, tol=1e-14):
    """Newton (secant) on det M(lam) = 0 starting at lam0."""
    f = lambda l: np.linalg.det(eig_matrix(l, mu1, k1, mu2, k2))
    l0, l1 = lam0, lam0 + 1e-3 * (1 + 1j)
    f0, f1 = f(l0), f(l1)
    for _ in range(100):
        if f1 == f0:
            break
        l2 = l1 - f1 * (l1 - l0) / (f1 - f0)
        l0, f0 = l1, f1
        l1, f1 = l2, f(l2)
        if abs(l1 - l0) < tol:
            break
    return l1


class TipField:
    def __init__(self, E1, nu1, E2, nu2, lam=None, verbose=False):
        self.mu1, self.k1 = mu_kappa(E1, nu1)
        self.mu2, self.k2 = mu_kappa(E2, nu2)
        self.eps = bimat_eps(E1, nu1, E2, nu2)
        lam_guess = 0.5 + 1j * self.eps if lam is None else lam
        self.lam = find_root(lam_guess, self.mu1, self.k1, self.mu2, self.k2)
        M = eig_matrix(self.lam, self.mu1, self.k1, self.mu2, self.k2)
        U, S, Vh = np.linalg.svd(M)
        self.sing = S
        lam = self.lam
        _, stt0, srt0 = _rows_stress(lam, 0.0)
        if S[-2] > 1e-8 * S[0]:
            # simple root (eps != 0): 1-D null space, the s_rt/s_tt ratio is fixed by the physics
            v = Vh[-1].conj()
            v = v / (stt0 @ v[0:4]) / np.sqrt(2 * np.pi)
        else:
            # eps = 0: double root (modes I and II decouple); pick the combination with
            # s_tt = 1/sqrt(2pi), s_rt = -i/sqrt(2pi) on theta = 0 (same convention as eps != 0)
            N = Vh[-2:].conj().T
            A = np.array([stt0 @ N[0:4], srt0 @ N[0:4]])
            alpha = np.linalg.solve(A, np.array([1.0, -1j]) / np.sqrt(2 * np.pi))
            v = N @ alpha
        self.v = v
        # ratio sigma_rt/sigma_tt on theta = 0 (expected -i for the r^(+i eps) convention)
        self.ratio0 = (srt0 @ v[0:4]) / (stt0 @ v[0:4])
        if verbose:
            print("eps =", self.eps, " lam =", self.lam, " lam-(1/2+i eps) =", self.lam - (0.5 + 1j * self.eps))
            print("singular values:", S[-3:], " ratio s_rt/s_tt at theta=0:", self.ratio0)

    def complex_field(self, x, y):
        """Complex eigenfield (normalised) at points; returns dict of complex arrays."""
        x = np.asarray(x, float)
        y = np.asarray(y, float)
        r = np.hypot(x, y)
        t = np.arctan2(y, x)
        lam = self.lam
        out = {k: np.zeros(x.shape, complex) for k in
               ("s11", "s22", "s12", "u1", "u2", "u11", "u12", "u21", "u22")}
        for mat, sel in ((1, t >= 0), (2, t < 0)):
            if not np.any(sel):
                continue
            tt, rr = t[sel], r[sel]
            if mat == 1:
                c, mu, kap = self.v[0:4], self.mu1, self.k1
            else:
                c, mu, kap = self.v[4:8], self.mu2, self.k2
            srr, stt, srt = (np.tensordot(c, R, axes=(0, 0)) for R in _rows_stress(lam, tt))
            ur, ut = (np.tensordot(c, R, axes=(0, 0)) for R in _rows_disp(lam, tt, mu, kap))
            dur, dut = (np.tensordot(c, R, axes=(0, 0)) for R in _rows_ddisp(lam, tt, mu, kap))
            rs = rr ** (lam - 1)
            rd = rr ** lam
            srr, stt, srt = srr * rs, stt * rs, srt * rs
            ur, ut, dur, dut = ur * rd, ut * rd, dur * rd, dut * rd
            C, S = np.cos(tt), np.sin(tt)
            out["s11"][sel] = srr * C * C + stt * S * S - 2 * srt * S * C
            out["s22"][sel] = srr * S * S + stt * C * C + 2 * srt * S * C
            out["s12"][sel] = (srr - stt) * S * C + srt * (C * C - S * S)
            u1 = ur * C - ut * S
            u2 = ur * S + ut * C
            out["u1"][sel], out["u2"][sel] = u1, u2
            dr_u1, dr_u2 = lam / rr * u1, lam / rr * u2
            dt_u1 = dur * C - ur * S - dut * S - ut * C
            dt_u2 = dur * S + ur * C + dut * C - ut * S
            out["u11"][sel] = C * dr_u1 - S / rr * dt_u1
            out["u12"][sel] = S * dr_u1 + C / rr * dt_u1
            out["u21"][sel] = C * dr_u2 - S / rr * dt_u2
            out["u22"][sel] = S * dr_u2 + C / rr * dt_u2
        return out

    def field(self, x, y, Ka):
        """Physical (real) auxiliary field for complex SIF Ka."""
        F = self.complex_field(x, y)
        return {k: np.real(Ka * v) for k, v in F.items()}


if __name__ == "__main__":
    # self-checks: root, convention, equilibrium, constitutive law, boundary/interface conditions
    rng = np.random.default_rng(0)
    for (E1, nu1, E2, nu2) in [(5.0, 0.22, 1.0, 0.35), (9.0, 1 / 3, 1.0, 1 / 3), (1.0, 0.3, 1.0, 0.3),
                               (1.0, 0.35, 5.0, 0.22)]:
        tf = TipField(E1, nu1, E2, nu2, verbose=True)
        Ka = 0.7 - 0.4j
        # equilibrium and compatibility by central differences
        pts = rng.uniform(-1, 1, (20, 2))
        pts = pts[np.abs(pts[:, 1]) > 0.05]
        hstep = 1e-5
        x, y = pts[:, 0], pts[:, 1]
        F = tf.field(x, y, Ka)
        Fxp, Fxm = tf.field(x + hstep, y, Ka), tf.field(x - hstep, y, Ka)
        Fyp, Fym = tf.field(x, y + hstep, Ka), tf.field(x, y - hstep, Ka)
        d = lambda key, A, B: (A[key] - B[key]) / (2 * hstep)
        eq1 = d("s11", Fxp, Fxm) + d("s12", Fyp, Fym)
        eq2 = d("s12", Fxp, Fxm) + d("s22", Fyp, Fym)
        g11, g12 = d("u1", Fxp, Fxm), d("u1", Fyp, Fym)
        g21, g22 = d("u2", Fxp, Fxm), d("u2", Fyp, Fym)
        smax = np.max(np.abs(F["s11"]))
        print("  equilibrium residual / max stress:", np.max(np.abs(eq1)) / smax, np.max(np.abs(eq2)) / smax)
        print("  grad u analytic vs FD:", max(np.max(np.abs(F["u11"] - g11)), np.max(np.abs(F["u12"] - g12)),
                                             np.max(np.abs(F["u21"] - g21)), np.max(np.abs(F["u22"] - g22))))
        # constitutive law
        mu = np.where(y > 0, tf.mu1, tf.mu2)
        kap = np.where(y > 0, tf.k1, tf.k2)
        nu = (3 - kap) / 4
        lmb = 2 * mu * nu / (1 - 2 * nu)
        tr = F["u11"] + F["u22"]
        e12 = 0.5 * (F["u12"] + F["u21"])
        r1 = F["s11"] - (2 * mu * F["u11"] + lmb * tr)
        r2 = F["s22"] - (2 * mu * F["u22"] + lmb * tr)
        r3 = F["s12"] - 2 * mu * e12
        print("  constitutive residual:", np.max(np.abs(np.r_[r1, r2, r3])) / smax)
        # faces and interface
        rr = np.array([0.1, 0.3, 0.9])
        Fu = tf.field(-rr, 0 * rr + 0.0, Ka)
        Fl = tf.field(-rr, -0.0 * np.ones_like(rr) - 1e-300, Ka)
        print("  face tractions upper:", np.max(np.abs(np.r_[Fu["s22"], Fu["s12"]])),
              " lower:", np.max(np.abs(np.r_[Fl["s22"], Fl["s12"]])))
        Fa = tf.field(rr, 0 * rr, Ka)
        Fb = tf.field(rr, -1e-300 * np.ones_like(rr), Ka)
        print("  interface jumps s22,s12,u1,u2:",
              [float(np.max(np.abs(Fa[k] - Fb[k]))) for k in ("s22", "s12", "u1", "u2")])
        trac = Fa["s22"] + 1j * Fa["s12"]
        target = Ka * rr ** (1j * tf.eps) / np.sqrt(2 * np.pi * rr)
        print("  interface traction vs Ka r^(i eps)/sqrt(2 pi r):", np.max(np.abs(trac - target)))
        target_m = Ka * rr ** (-1j * tf.eps) / np.sqrt(2 * np.pi * rr)
        print("  (same with r^(-i eps)):", np.max(np.abs(trac - target_m)))
