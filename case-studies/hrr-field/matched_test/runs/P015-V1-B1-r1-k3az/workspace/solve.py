#!/usr/bin/env python3
"""Mode III (anti-plane) HRR-type crack-tip field for a pure power-law shear material.

Law: gamma_vec = (gamma0/tau0^n) |tau|^{n-1} tau_vec   (engineering shear)
Crack along negative x-axis; faces theta = +-pi traction free; w(x,-y) = -w(x,y).

Separable field from stress function Phi = A r^lam f(theta):
    tau_rz = (1/r) Phi_,theta = A r^{-p} f'(theta)
    tau_tz = -Phi_,r          = -A lam r^{-p} f(theta),   p = 1 - lam.
Compatibility curl(gamma)=0 gives the eigenvalue ODE (g = f'^2 + lam^2 f^2):
    f'' = -[(n-1) lam^2 f f'^2 + lam q g f] / [n f'^2 + lam^2 f^2],  q = n(lam-1)+1.
BCs: f'(0) = 0 (tau_rz = 0 ahead, by antisymmetry), f(pi) = 0 (traction-free face).
Sign: f(0) = -1 so that w > 0 on the upper face.

J = c A^{n+1} I_n,  c = gamma0/tau0^n,
I_n = int_{-pi}^{pi} [(n/(n+1)) T^{n+1} cos th - T^{n-1} f'(f' cos th + lam f sin th)] dth,
T = sqrt(f'^2 + lam^2 f^2).  =>  A = [J/(c I_n)]^{1/(n+1)}.

w(r,th) = c A^n (r^p/p) T^{n-1} f'   (satisfies w_,r = gamma_rz, (1/r) w_,th = gamma_tz).
"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq


def solve_angular(n, n_th=200001):
    """Solve the eigenvalue problem. Returns lam, theta grid on [-pi,pi], f, fp, T, I_n."""
    def rhs(lam):
        q = n * (lam - 1.0) + 1.0
        def f_(th, y):
            f, fp = y
            g = fp * fp + lam * lam * f * f
            fpp = -((n - 1.0) * lam * lam * f * fp * fp + lam * q * g * f) / (
                n * fp * fp + lam * lam * f * f)
            return [fp, fpp]
        return f_

    def f_at_pi(lam):
        try:
            sol = solve_ivp(rhs(lam), (0.0, np.pi), [-1.0, 0.0],
                            rtol=3e-13, atol=1e-15)
            return sol.y[0, -1]
        except Exception:
            return np.nan

    lams = np.linspace(0.3, 0.999, 300)
    vals = [f_at_pi(l) for l in lams]
    lam_root = None
    for i in range(len(lams) - 1):
        if vals[i] * vals[i + 1] < 0:
            lam_root = brentq(f_at_pi, lams[i], lams[i + 1],
                              xtol=1e-14, rtol=1e-14)
            break
    if lam_root is None:
        raise RuntimeError("no eigenvalue bracket found")
    sol = solve_ivp(rhs(lam_root), (0.0, np.pi), [-1.0, 0.0],
                    rtol=3e-13, atol=1e-15, dense_output=True)
    th_h = np.linspace(0.0, np.pi, n_th)
    fh, fph = sol.sol(th_h)
    th = np.concatenate([-th_h[:0:-1], th_h])      # -pi .. pi
    f = np.concatenate([fh[:0:-1], fh])            # f even
    fp = np.concatenate([-fph[:0:-1], fph])        # f' odd
    T = np.hypot(fp, lam_root * f)
    integrand = ((n / (n + 1.0)) * T ** (n + 1) * np.cos(th)
                 - T ** (n - 1) * fp * (fp * np.cos(th) + lam_root * f * np.sin(th)))
    I_n = np.trapezoid(integrand, th)
    return lam_root, th, f, fp, T, I_n


def interp(th_grid, arr, thq):
    return float(np.interp(thq, th_grid, arr))



def main():
    # ---- given data ------------------------------------------------------
    tau0 = 80.0        # MPa
    gam0 = 0.002
    n = 3.5
    J = 0.120          # N/mm
    p_ref = 1.0 / (n + 1.0)

    # ---- n = 1 pipeline test against the exact elastic field -------------
    lam1, th1, f1, fp1, T1, I1 = solve_angular(1.0)
    f_exact = -np.cos(th1 / 2.0)
    fp_exact = 0.5 * np.sin(th1 / 2.0)
    print("=== n=1 check vs exact elastic solution ===")
    print(f"lam = {lam1:.12f}  (exact 0.5)")
    print(f"max|f - (-cos(th/2))|   = {np.max(np.abs(f1 - f_exact)):.3e}")
    print(f"max|fp - 0.5 sin(th/2)| = {np.max(np.abs(fp1 - fp_exact)):.3e}")
    print(f"I_1 = {I1:.12f}  (exact pi/4 = {np.pi/4:.12f})")

    # ---- actual material n = 3.5 ------------------------------------------
    lam, th, f, fp, T, I_n = solve_angular(n)
    p = 1.0 - lam
    c = gam0 / tau0 ** n
    A = (J / (c * I_n)) ** (1.0 / (n + 1.0))
    print("\n=== n = 3.5 eigenvalue / amplitude ===")
    print(f"lam (shot)      = {lam:.12f}")
    print(f"n/(n+1)         = {n/(n+1.0):.12f}")
    print(f"p = 1 - lam     = {p:.12f}   (1/(n+1) = {p_ref:.12f})")
    print(f"I_n             = {I_n:.10f}")
    print(f"c = gam0/tau0^n = {c:.10e} MPa^-n")
    print(f"A               = {A:.10f} MPa mm^p")
    print(f"f(pi) residual  = {f[-1]:.3e}   f'(pi) = {fp[-1]:.10f}")

    # ---- field quantities --------------------------------------------------
    def T_of(thq):
        return interp(th, T, thq)

    def f_of(thq):
        return interp(th, f, thq)

    def fp_of(thq):
        return interp(th, fp, thq)

    def tau_magnitude(r, thq):
        return A * r ** (-p) * T_of(thq)

    def tau_rz(r, thq):
        return A * r ** (-p) * fp_of(thq)

    def tau_tz(r, thq):
        return -A * lam * r ** (-p) * f_of(thq)

    def w(r, thq):
        return c * A ** n * r ** p / p * T_of(thq) ** (n - 1) * fp_of(thq)

    print("\n=== (1) relative face slip  Delta w = w(r,+pi) - w(r,-pi) ===")
    for r in (0.0200, 0.0800):
        dw = w(r, np.pi) - w(r, -np.pi)
        print(f"r = {r:.4f} mm :  w(r,+pi) = {w(r, np.pi):.8e} mm, "
              f"Delta w = {dw:.8e} mm")

    print("\n=== (2) signed stresses at r = 0.0200 mm, theta = 60 deg ===")
    r0, thq = 0.0200, np.pi / 3.0
    print(f"tau_rz = {tau_rz(r0, thq):.8f} MPa")
    print(f"tau_tz = {tau_tz(r0, thq):.8f} MPa")

    print("\n=== (3) radii where |gamma| = 0.0100 ===")
    gam_star = 0.0100
    tau_star = tau0 * (gam_star / gam0) ** (1.0 / n)
    print(f"tau_* (=tau0*(0.01/gam0)^(1/n)) = {tau_star:.6f} MPa")
    for deg, thq in ((0, 0.0), (90, np.pi / 2.0), (180, np.pi)):
        r_star = (A * T_of(thq) / tau_star) ** (1.0 / p)
        print(f"theta = {deg:3d} deg :  T = {T_of(thq):.8f},  r = {r_star:.8f} mm")

    # ================= verification ========================================
    print("\n=== verification ===")
    # (a) J from the field on circles, full Cartesian expressions
    def J_on_circle(R, k=200001):
        a = np.linspace(-np.pi, np.pi, k)
        trz = A * R ** (-p) * np.interp(a, th, fp)
        ttz = -A * lam * R ** (-p) * np.interp(a, th, f)
        tau = A * R ** (-p) * np.interp(a, th, T)
        grz = c * tau ** (n - 1) * trz
        gtz = c * tau ** (n - 1) * ttz
        gxz = grz * np.cos(a) - gtz * np.sin(a)      # w_,x (scalar: no basis trap)
        W = (n / (n + 1.0)) * c * tau ** (n + 1)     # stored work density
        tz = trz                                      # traction on circle (N = e_r)
        return np.trapezoid((W * np.cos(a) - tz * gxz) * R, a)

    for R in (1e-3, 0.02, 0.08, 0.25):
        Jc = J_on_circle(R)
        print(f"J on circle r = {R:7.4f} mm : {Jc:.10f} N/mm  "
              f"(rel. err {abs(Jc - J) / J:.2e})")

    # (b) compatibility: (1/r) w_,theta vs gamma_tz on a fine grid
    R = 0.02
    a = np.linspace(-np.pi + 1e-3, np.pi - 1e-3, 4001)
    wv = np.array([w(R, t) for t in a])
    gtz_fd = np.gradient(wv, a) / R
    trz = A * R ** (-p) * np.interp(a, th, fp)
    ttz = -A * lam * R ** (-p) * np.interp(a, th, f)
    tau = A * R ** (-p) * np.interp(a, th, T)
    gtz = c * tau ** (n - 1) * ttz
    rel = np.max(np.abs(gtz_fd[5:-5] - gtz[5:-5])) / np.max(np.abs(gtz))
    print(f"compat (1/r)w_,th vs gamma_tz: max rel diff = {rel:.2e}")

    # (c) equilibrium in Cartesian by finite differences: tau_xz,x + tau_yz,y = 0
    def tau_cart(x, y):
        r = np.hypot(x, y)
        a = np.arctan2(y, x)
        trz = A * r ** (-p) * np.interp(a, th, fp)
        ttz = -A * lam * r ** (-p) * np.interp(a, th, f)
        txz = trz * np.cos(a) - ttz * np.sin(a)
        tyz = trz * np.sin(a) + ttz * np.cos(a)
        return txz, tyz

    h = 1e-5
    pts = [(0.02, 0.01), (0.01, 0.03), (0.05, 0.02), (-0.03, 0.04)]
    for (x, y) in pts:
        dtxx = (tau_cart(x + h, y)[0] - tau_cart(x - h, y)[0]) / (2 * h)
        dtyy = (tau_cart(x, y + h)[1] - tau_cart(x, y - h)[1]) / (2 * h)
        scale = A * np.hypot(x, y) ** (-p) / np.hypot(x, y)
        print(f"equilibrium FD at ({x},{y}): div = {dtxx + dtyy:.3e} MPa/mm "
              f"(scale {scale:.3e}, rel {abs(dtxx+dtyy)/scale:.2e})")

    # (d) boundary residual: tau_tz on the faces
    print(f"face traction tau_tz(r=0.02, +pi) = {tau_tz(0.02, np.pi):.3e} MPa")
    print(f"face traction tau_tz(r=0.02, -pi) = {tau_tz(0.02, -np.pi):.3e} MPa")

    # (e) |gamma| definition consistency: |gamma| = gam0 (tau/tau0)^n
    r1, t1 = 0.02, np.pi / 3
    tau_v = tau_magnitude(r1, t1)
    gv = c * tau_v ** n
    print(f"|gamma| at r=0.02, th=60deg: {gv:.6f} vs gam0*(tau/tau0)^n = "
          f"{gam0 * (tau_v / tau0) ** n:.6f}")

    # (f) displacement finite-difference check: w_,r vs gamma_rz
    for (rr, tt) in ((0.02, np.pi / 3), (0.08, 2.5)):
        hr = 1e-6 * rr
        dwr = (w(rr + hr, tt) - w(rr - hr, tt)) / (2 * hr)
        grz = c * tau_magnitude(rr, tt) ** (n - 1) * tau_rz(rr, tt)
        print(f"w_,r vs gamma_rz at r={rr}, th={tt:.3f}: {dwr:.8e} vs {grz:.8e}")


if __name__ == "__main__":
    main()
