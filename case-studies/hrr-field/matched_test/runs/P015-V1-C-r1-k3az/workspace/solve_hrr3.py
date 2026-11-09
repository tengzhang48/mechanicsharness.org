"""Mode III stationary crack in power-law hardening material (HRR-type field).

Constitutive law (deformation plasticity, anti-plane):
    gamma_vec = gamma0 * (tau/tau0)^(n-1) * (tau_vec/tau0)   [engineering shear]
Separable crack-tip field (stress ~ r^(-1/(n+1)), strain ~ r^(-n/(n+1))):
    tau_r   = tau0 * C * r^t * V(theta),   t = -1/(n+1)
    tau_th  = tau0 * C * r^t * U(theta)
    w       = (n+1) * gamma0 * C^n * S^(n-1) * V * r^s,   s = 1/(n+1)
where  V = -U'/c,  c = n/(n+1),  S = sqrt(U^2+V^2),  and U solves
    d/dtheta[ M U' ] = -k M U,   k = n/(n+1)^2,  M = (U^2 + U'^2/c^2)^((n-1)/2)
    U'(0) = 0 (antisymmetry),  U(pi) = 0 (traction-free face).
Amplitude C from J:  J = gamma0 * tau0 * C^(n+1) * I_n,
    I_n = int_{-pi}^{pi} [ n/(n+1) S^(n+1) cos - S^(n-1) V^2 cos + S^(n-1) V U sin ] dtheta
"""
import numpy as np
from scipy.integrate import solve_ivp

PI = np.pi
OUT = "<matched-runs>/runs/P015-V1-C-r1-k3az/results.txt"


def solve_shape(n, rtol=1e-12, atol=1e-14):
    c = n / (n + 1.0)
    k = n / (n + 1.0) ** 2
    a = (n - 1.0) / 2.0

    def rhs(th, y):
        U1, U2 = y
        Q = U1 * U1 + (U2 / c) ** 2
        # U2' = -U1 (k Q + 2a U2^2) / (Q + 2a U2^2 / c^2)
        U2p = -U1 * (k * Q + 2.0 * a * U2 * U2) / (Q + 2.0 * a * U2 * U2 / (c * c))
        return [U2, U2p]

    sol = solve_ivp(rhs, (0.0, PI), [1.0, 0.0], dense_output=True,
                    rtol=rtol, atol=atol, max_step=0.01)
    return sol, c, k, a


def main():
    lines = []
    log = lines.append

    # ---------- verification 0: elastic limit n = 1 has exact K-field ----------
    sol1, c1, k1, a1 = solve_shape(1.0)
    tht = np.linspace(0, PI, 2001)
    y1 = sol1.sol(tht)
    errU = np.max(np.abs(y1[0] - np.cos(tht / 2.0)))
    errV = np.max(np.abs(-y1[1] / c1 - np.sin(tht / 2.0)))
    log(f"[check n=1] max|U-cos(th/2)| = {errU:.3e}, max|V-sin(th/2)| = {errV:.3e}, "
        f"U(pi) = {sol1.sol(PI)[0]:.3e}")

    # ---------- given data ----------
    n = 3.5
    tau0 = 80.0       # MPa
    gam0 = 0.002
    J = 0.120         # N/mm

    t = -1.0 / (n + 1.0)   # stress exponent
    s = 1.0 / (n + 1.0)    # displacement exponent

    sol, c, k, a = solve_shape(n)
    Upi = sol.sol(PI)[0]
    log(f"[check n=3.5] U(pi) = {Upi:.3e}  (traction-free face residual, normalized U(0)=1)")

    # dense angular grid on [0, pi]
    th = np.linspace(0, PI, 40001)
    y = sol.sol(th)
    U = y[0]
    V = -y[1] / c
    S = np.hypot(U, V)

    # ---------- J-integral coefficient ----------
    integrand = (n / (n + 1.0)) * S ** (n + 1) * np.cos(th) \
        - S ** (n - 1) * V ** 2 * np.cos(th) \
        + S ** (n - 1) * V * U * np.sin(th)
    I_n = 2.0 * np.trapezoid(integrand, th)   # integrand is even in theta
    C = (J / (gam0 * tau0 * I_n)) ** (1.0 / (n + 1.0))
    log(f"I_n = {I_n:.10f},  C = {C:.10f} mm^(1/(n+1))")

    # ---------- symmetry extension to [-pi, pi] ----------
    th_full = np.concatenate((-th[:0:-1], th))          # -pi..pi
    U_full = np.concatenate((U[:0:-1], U))              # U even
    V_full = np.concatenate((-V[:0:-1], V))             # V odd
    S_full = np.concatenate((S[:0:-1], S))              # S even

    def Uf(x): return np.interp(x, th_full, U_full)
    def Vf(x): return np.interp(x, th_full, V_full)
    def Sf(x): return np.interp(x, th_full, S_full)

    def tau_r(r, x):  return C * r ** t * tau0 * Vf(x)
    def tau_th(r, x): return C * r ** t * tau0 * Uf(x)
    def gam_mag(r, x): return gam0 * (C * r ** t * Sf(x)) ** n
    def g_fun(x): return (n + 1.0) * gam0 * C ** n * Sf(x) ** (n - 1) * Vf(x)
    def w(r, x): return g_fun(x) * r ** s

    # ---------- independent J recomputation on circle (verification) ----------
    def J_circle(r):
        x = np.linspace(-PI, PI, 80001)
        tt = np.hypot(tau_r(r, x), tau_th(r, x))            # |tau|
        gg = gam0 * (tt / tau0) ** n                        # |gamma|
        W = n / (n + 1.0) * tt * gg                         # stored work density
        g_r = gg * Vf(x) / Sf(x)                            # gamma_r
        g_th = gg * Uf(x) / Sf(x)                           # gamma_theta
        wx = g_r * np.cos(x) - g_th * np.sin(x)             # w_,x
        tz = tau_r(r, x)                                    # t_z = tau . N, N = e_r
        return np.trapezoid(W * np.cos(x) - tz * wx, x) * r

    for rr in (0.02, 0.08, 0.29):
        log(f"[check J] r = {rr:.2f} mm : J = {J_circle(rr):.8f} N/mm (input 0.120)")

    # ---------- compatibility residual: g'(theta) ?= gamma0 C^n S^(n-1) U ----------
    gg_num = g_fun(th)
    dg = np.gradient(gg_num, th)
    rhs_c = gam0 * C ** n * S ** (n - 1) * U
    compat = np.max(np.abs(dg - rhs_c)) / np.max(np.abs(rhs_c))
    log(f"[check compat] max rel. residual of dg/dtheta vs gamma_theta*r... = {compat:.3e}")

    # ---------- 1) face slip ----------
    for rr in (0.0200, 0.0800):
        dw = w(rr, PI) - w(rr, -PI)
        log(f"slip  Dw({rr:.4f} mm) = {dw:.8e} mm")

    # ---------- 2) signed stresses at r=0.02, theta=60 deg ----------
    th60 = np.deg2rad(60.0)
    tr = tau_r(0.0200, th60)
    tt = tau_th(0.0200, th60)
    log(f"tau_rz(0.02, 60deg)  = {tr:.8e} MPa")
    log(f"tau_thz(0.02, 60deg) = {tt:.8e} MPa")

    # ---------- 3) radii where |gamma| = 0.0100 ----------
    gtar = 0.0100
    for name, x in (("0 deg", 0.0), ("90 deg", PI / 2.0), ("180 deg (upper face)", PI)):
        Sx = Sf(x)
        # gam0 * (C r^t Sx)^n = gtar, t = -1/(n+1)  ->  r = [ gam0 (C Sx)^n / gtar ]^((n+1)/n)
        rr = (gam0 * (C * Sx) ** n / gtar) ** ((n + 1.0) / n)
        chk = gam_mag(rr, x)
        log(f"r(|gamma|={gtar}) at theta={name}: r = {rr:.8e} mm  (check gamma={chk:.6f}, S={Sx:.6f})")

    # ---------- extra data for the report ----------
    log(f"[data] C^n = {C**n:.8f} mm^(n/(n+1));  U(60deg) = {Uf(PI/3):.6f}, V(60deg) = {Vf(PI/3):.6f}, "
        f"S(60deg) = {Sf(PI/3):.6f}")
    log(f"[data] V(pi) = {Vf(PI):.6f}, S(90deg) = {Sf(PI/2):.6f}, S(pi) = {Sf(PI):.6f}, "
        f"max S = {S.max():.6f} at theta = {np.degrees(th[np.argmax(S)]):.2f} deg")
    log(f"[data] r^t at 0.02 mm = {0.02**t:.6f}; w(r,pi) amplitude 2*g(pi) = {2*g_fun(PI):.8e} mm^(1-1/(n+1))")

    with open(OUT, "w") as f:
        f.write("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
