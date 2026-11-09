#!/usr/bin/env python3
"""Mode III (anti-plane shear) crack-tip field for a pure power-law material.

Constitutive law: gamma_vec parallel to tau_vec, |gamma| = gamma0*(tau/tau0)^n,
so tau = tau0*(|gamma|/gamma0)^(1/n). Anti-plane: w(r,theta), gamma = grad(w),
div(tau) = 0.

Separable HRR field: w = A r^q f(theta), q = 1/(n+1). The exponent q is fixed by
J path independence: W ~ tau*gamma ~ r^((q-1)(n+1)/n), J over a circle ~ W*r must
be r-independent -> (q-1)(n+1)/n = -1 -> q = 1/(n+1). Stress tau ~ r^s,
s = (q-1)/n = -1/(n+1).

Equilibrium (r tau_rz)_,r + tau_tz,th = 0 gives the angular ODE (p = f'):
    G = sqrt(q^2 f^2 + p^2), eta = p^2/G^2,
    p' = -q^2 f (n^2 - (n-1) eta) / (n - (n-1) eta)
BCs: f(0)=0 (antisymmetry), f'(pi)=0 (traction-free face). ODE homogeneous in f,
normalize f'(0)=1; amplitude A set by J:
    J = tau0 gamma0 (A/gamma0)^((n+1)/n) I,
    I = int_{-pi}^{pi} [(n/(n+1)) G^((n+1)/n) cos(th)
        - G^(1/n-1) q f (q f cos(th) - f' sin(th))] dth
"""
import numpy as np
from scipy.integrate import solve_ivp, trapezoid


def solve_angular(n, npts=40001, rtol=1e-12, atol=1e-14):
    """Angular functions f, p=f', G on [0, pi] with f(0)=0, f'(0)=1."""
    q = 1.0 / (n + 1.0)

    def rhs(th, y):
        f, p = y
        G2 = q * q * f * f + p * p
        eta = p * p / G2
        fp = -q * q * f * (n * n - (n - 1.0) * eta) / (n - (n - 1.0) * eta)
        return [p, fp]

    grid = np.linspace(0.0, np.pi, npts)
    sol = solve_ivp(rhs, (0.0, np.pi), [0.0, 1.0], t_eval=grid,
                    rtol=rtol, atol=atol, method="DOP853", dense_output=True)
    f, p = sol.y
    G = np.sqrt(q * q * f * f + p * p)
    return grid, f, p, G, q, sol.sol


def I_integral(n, grid, f, p, G, q):
    """Angular integral I (integrand even in theta: f odd, p and G even)."""
    a = 1.0 / n - 1.0
    integ = (n / (n + 1.0)) * G ** ((n + 1.0) / n) * np.cos(grid) \
        - G ** a * q * f * (q * f * np.cos(grid) - p * np.sin(grid))
    return 2.0 * trapezoid(integ, grid)


def amplitude(n, tau0, gam0, J, I):
    return gam0 * (J / (tau0 * gam0 * I)) ** (n / (n + 1.0))


def fields(n, tau0, gam0, A, grid, f, p, G, q, r, theta):
    """Angular functions interpolated to theta in [0, pi]; evaluate at radius r."""
    fq = np.interp(theta, grid, f)
    pq = np.interp(theta, grid, p)
    Gq = np.sqrt(q * q * fq * fq + pq * pq)
    s = -1.0 / (n + 1.0)
    C = tau0 * (A / gam0) ** (1.0 / n)
    tau_mag = C * r ** s * Gq ** (1.0 / n)
    gam_mag = A * r ** (q - 1.0) * Gq
    tau_rz = C * r ** s * Gq ** (1.0 / n - 1.0) * q * fq
    tau_tz = C * r ** s * Gq ** (1.0 / n - 1.0) * pq
    return tau_mag, gam_mag, tau_rz, tau_tz


def J_on_circle(n, tau0, gam0, A, grid, f, p, G, q, r):
    """Full J on a circle of radius r over [-pi, pi]; integrand even -> 2x upper half."""
    C = tau0 * (A / gam0) ** (1.0 / n)
    s = -1.0 / (n + 1.0)
    tau_rz = C * r ** s * G ** (1.0 / n - 1.0) * q * f
    gam_r = A * q * r ** (q - 1.0) * f
    gam_t = A * r ** (q - 1.0) * p
    gam_mag = A * r ** (q - 1.0) * G
    W = n / (n + 1.0) * tau0 * gam0 * (gam_mag / gam0) ** ((n + 1.0) / n)
    wx = gam_r * np.cos(grid) - gam_t * np.sin(grid)
    integ = W * np.cos(grid) - tau_rz * wx   # t_z = tau_rz (radial outward normal)
    return 2.0 * trapezoid(integ * r, grid)


def main():
    tau0 = 80.0      # MPa = N/mm^2
    gam0 = 0.002     # -
    n = 3.5          # -
    J_app = 0.120    # N/mm

    grid, f, p, G, q, solfun = solve_angular(n)
    print(f"n = {n}, q = 1/(n+1) = {q:.10f}, stress exponent s = {-1.0/(n+1.0):.10f}")
    print(f"traction-free check: f'(pi) = {p[-1]:.3e}  (should be ~0)")
    print(f"f(pi) = {f[-1]:.10f},  G(0) = {G[0]:.10f},  G(pi) = {G[-1]:.10f}")

    I = I_integral(n, grid, f, p, G, q)
    print(f"I (normalization f'(0)=1) = {I:.10f}")
    A = amplitude(n, tau0, gam0, J_app, I)
    print(f"Amplitude A = {A:.10e} mm^(n/(n+1))  [w = A r^q f]")

    # verification 1: n=1 limit vs exact elastic Mode III field
    # mu = tau0/gam0, J = K^2/(2 mu), w = (K/mu) sqrt(2 r/pi) sin(th/2);
    # in f'(0)=1 normalization: f = 2 sin(th/2), I = pi, A = sqrt(gam0 J/(pi tau0)).
    g1, f1, p1, G1, q1, _ = solve_angular(1.0)
    err_f = np.max(np.abs(f1 - 2.0 * np.sin(g1 / 2.0)))
    err_p = np.max(np.abs(p1 - np.cos(g1 / 2.0)))
    I1 = I_integral(1.0, g1, f1, p1, G1, q1)
    A1 = amplitude(1.0, tau0, gam0, J_app, I1)
    A1_exact = np.sqrt(gam0 * J_app / (np.pi * tau0))
    print("\n[n=1 sanity check vs exact elastic Mode III field]")
    print(f"max|f - 2 sin(th/2)| = {err_f:.3e}, max|f' - cos(th/2)| = {err_p:.3e}")
    print(f"I1 = {I1:.12f} vs pi, rel err = {abs(I1 - np.pi) / np.pi:.3e}")
    print(f"A1 rel err vs sqrt(gam0 J/(pi tau0)) = {abs(A1 - A1_exact) / A1_exact:.3e}")

    # verification 2: J path independence of the computed n=3.5 field
    print("\n[J on circles from the computed n=3.5 field]")
    for r in (0.005, 0.02, 0.1, 0.3):
        Jr = J_on_circle(n, tau0, gam0, A, grid, f, p, G, q, r)
        print(f"  r = {r:6.3f} mm:  J = {Jr:.10f} N/mm, rel dev = {(Jr - J_app) / J_app:+.2e}")

    print("\n================ REQUESTED NUMBERS (n = 3.5) ================")
    for r in (0.0200, 0.0800):
        dw = 2.0 * A * r ** q * f[-1]
        print(f"1) r = {r:.4f} mm:  Delta w = 2 A r^q f(pi) = {dw:.6e} mm")

    r2, th2 = 0.0200, np.deg2rad(60.0)
    tau_mag, gam_mag, tau_rz, tau_tz = fields(n, tau0, gam0, A, grid, f, p, G, q, r2, th2)
    print(f"2) r = 0.0200 mm, theta = 60 deg:  tau_rz = {tau_rz:.6f} MPa, "
          f"tau_tz = {tau_tz:.6f} MPa  (|tau| = {tau_mag:.6f}, |gamma| = {gam_mag:.3e})")

    gam_target = 0.0100
    for thdeg in (0.0, 90.0, 180.0):
        th = np.deg2rad(thdeg)
        Gt = np.interp(th, grid, G)
        rt = (A * Gt / gam_target) ** ((n + 1.0) / n)
        print(f"3) theta = {thdeg:5.1f} deg:  r(|gamma|=0.01) = {rt:.6e} mm  (G = {Gt:.8f})")

    # verification 3: residuals in Cartesian (equilibrium, constitutive, face traction)
    print("\n[Cartesian FD checks of the n=3.5 field]")
    rng = np.random.default_rng(0)
    hh = 1e-7
    C = tau0 * (A / gam0) ** (1.0 / n)
    s = -1.0 / (n + 1.0)

    def tau_cart(x, y):
        rr = np.hypot(x, y)
        th = np.arctan2(y, x)
        ath, sign = abs(th), np.sign(th)
        fq, pq = solfun(ath)   # dense output: smooth angular values
        Gq = np.sqrt(q * q * fq * fq + pq * pq)
        tr = C * rr ** s * Gq ** (1.0 / n - 1.0) * q * fq * sign  # tau_rz is odd in theta
        tt = C * rr ** s * Gq ** (1.0 / n - 1.0) * pq             # tau_tz even
        tx = tr * np.cos(th) - tt * np.sin(th)   # signed theta in the rotation
        ty = tr * np.sin(th) + tt * np.cos(th)
        return tx, ty

    eq_res, law_res, face_res = [], [], []
    for i in range(300):
        rr = rng.uniform(1e-4, 0.29)
        th = rng.uniform(-np.pi + 1e-3, np.pi - 1e-3)
        x, y = rr * np.cos(th), rr * np.sin(th)
        div = ((tau_cart(x + hh, y)[0] - tau_cart(x - hh, y)[0]) / (2 * hh)
               + (tau_cart(x, y + hh)[1] - tau_cart(x, y - hh)[1]) / (2 * hh))
        tx0, ty0 = tau_cart(x, y)
        eq_res.append(abs(div) / (np.hypot(tx0, ty0) / rr))
        # constitutive residual: |gamma| from law vs field
        tau_m = np.hypot(tx0, ty0)
        gam_law = gam0 * (tau_m / tau0) ** n
        ath = abs(np.arctan2(y, x))
        fq, pq = solfun(ath)
        gam_field = A * rr ** (q - 1.0) * np.sqrt(q * q * fq * fq + pq * pq)
        law_res.append(abs(gam_law - gam_field) / gam_field)
        # face traction residual, upper face (outward normal (0,-1))
        xf = -rr
        txf, tyf = tau_cart(xf, 1e-9)
        face_res.append(abs(-tyf) / np.hypot(txf, tyf))
    print(f"  max |div tau|/(|tau|/r) over 300 pts = {max(eq_res):.3e}")
    print(f"  max constitutive residual            = {max(law_res):.3e}")
    print(f"  max face traction |t_z|/|tau|        = {max(face_res):.3e}")


if __name__ == "__main__":
    main()

