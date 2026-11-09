"""Mode III (anti-plane) crack-tip field in a pure power-law shear material.

Material: gamma = (g0/t0) (tau/t0)^(n-1) tau  (vector law), tau0=80 MPa, g0=0.002, n=3.5.
Crack along negative x-axis, faces traction-free at theta=+-pi, antisymmetric w.

Model: separable leading field  tau = A r^-p T(theta),  w = r^s F(theta).
 - equilibrium in polar: (1-p) T_r + T_theta' = 0
 - compatibility (gamma = grad w, collinear with tau):
       d/dtheta [ m^(n-1) T_r ] = s m^(n-1) T_theta,   m = |T|,  s = 1 - n p
 - BCs: T_r(0)=0 (antisymmetry), T_theta(pi)=0 (traction-free face)
 - p is an eigenvalue (shoot); for a crack, finite nonzero J also requires p=1/(n+1).
 - Amplitude fixed by J = 0.120 N/mm.
"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

TAU0, G0 = 80.0, 0.002
J_APP = 0.120


def rhs_factory(n, p):
    s = 1.0 - n * p

    def rhs(th, y):
        Tr, Tt = y
        m = np.hypot(Tr, Tt)
        dTt = -(1.0 - p) * Tr
        mn1 = m ** (n - 1.0)
        mn3 = m ** (n - 3.0)
        A = mn1 + (n - 1.0) * mn3 * Tr * Tr
        b = s * mn1 * Tt - (n - 1.0) * mn3 * Tr * Tt * dTt
        return [b / A, dTt]

    return rhs, s


def shoot_Ttheta_pi(n, p):
    rhs, s = rhs_factory(n, p)
    sol = solve_ivp(rhs, (0.0, np.pi), [0.0, 1.0], rtol=1e-12, atol=1e-14,
                    dense_output=True, max_step=0.01)
    return sol.y[1, -1], sol


def solve_field(n, p_bracket=(0.12, 0.30)):
    f = lambda p: shoot_Ttheta_pi(n, p)[0]
    p = brentq(f, *p_bracket, xtol=1e-14, rtol=1e-15)
    _, sol = shoot_Ttheta_pi(n, p)
    return p, sol


def angular(n, p, sol, nth=20001):
    """Normalized angular functions on [0, pi] with T_theta(0)=1."""
    rhs, s = rhs_factory(n, p)
    th = np.linspace(0.0, np.pi, nth)
    Tr, Tt = sol.sol(th)
    m = np.hypot(Tr, Tt)
    # strain angular functions Gamma = (g0/t0^n) m^(n-1) T  (r^(s-1) factored out)
    Gr = (G0 / TAU0**n) * m**(n - 1.0) * Tr
    Gt = (G0 / TAU0**n) * m**(n - 1.0) * Tt
    F = Gr / s                       # w = r^s F(theta)
    return dict(th=th, Tr=Tr, Tt=Tt, m=m, Gr=Gr, Gt=Gt, F=F, s=s)


def In_integral(n, a):
    """J = (g0/t0^n) A^(n+1) * I_n  with the normalization max_theta ... T_theta(0)=1."""
    th, Tr, Tt, m = a["th"], a["Tr"], a["Tt"], a["m"]
    h = ((n / (n + 1.0)) * m**(n + 1.0) * np.cos(th)
         - m**(n - 1.0) * Tr**2 * np.cos(th)
         + m**(n - 1.0) * Tr * Tt * np.sin(th))
    return 2.0 * np.trapezoid(h, th)          # integrand even in theta


# ---------------------------------------------------------------- n=1 regression vs closed form
def regression_n1():
    n = 1.0
    p, sol = solve_field(n, p_bracket=(0.30, 0.60))
    a = angular(n, p, sol)
    I1 = In_integral(n, a)
    s = a["s"]
    # closed form: T_r = sin(th/2), T_t = cos(th/2), I_1 = pi, p = 1/2
    err_p = abs(p - 0.5)
    err_T = max(np.max(np.abs(a["Tr"] - np.sin(a["th"] / 2))),
                np.max(np.abs(a["Tt"] - np.cos(a["th"] / 2))))
    err_I = abs(I1 - np.pi)
    # end-to-end displacement check against K_III solution, mu = tau0/g0
    A = np.sqrt(J_APP * TAU0 / (G0 * I1))          # A_T = K/sqrt(2 pi)
    F_pi = A**n * (G0 / TAU0**n) * a["m"][-1]**(n - 1) * a["Tr"][-1] / s
    dw_num = 2 * np.sqrt(0.05) * F_pi
    K = np.sqrt(2 * (TAU0 / G0) * J_APP)           # J = K^2/(2 mu)
    dw_exact = 2 * (K * G0 / TAU0) * np.sqrt(2 * 0.05 / np.pi)
    print(f"[n=1] p*={p:.12f} (err {err_p:.1e}), T err {err_T:.2e}, I_1={I1:.12f} (err {err_I:.1e})")
    print(f"[n=1] dw(0.05): numeric {dw_num:.9f}, closed form {dw_exact:.9f}, rel {abs(dw_num/dw_exact-1):.2e}")
    return err_p, err_T, err_I, abs(dw_num / dw_exact - 1)


# ---------------------------------------------------------------- main case n=3.5
def main():
    n = 3.5
    p, sol = solve_field(n)
    a = angular(n, p, sol)
    s, th = a["s"], a["th"]
    print(f"n={n}: eigenvalue p* = {p:.12f}, 1/(n+1) = {1/(n+1):.12f}, "
          f"|diff| = {abs(p - 1/(n+1)):.2e}")
    print(f"face residual |T_theta(pi)| = {abs(a['Tt'][-1]):.2e} (normalization T_theta(0)=1)")
    # ODE residual: compatibility  F' = G_t
    dF = np.gradient(a["F"], th)
    print(f"compatibility residual max|F' - Gamma_theta|/max|Gamma_theta| = "
          f"{np.max(np.abs(dF - a['Gt'])) / np.max(np.abs(a['Gt'])):.2e} (finite diff)")

    I_n = In_integral(n, a)
    A = (J_APP * TAU0**n / (G0 * I_n)) ** (1.0 / (n + 1.0))   # tau = A r^-p T, in MPa mm^p
    print(f"I_n = {I_n:.9f},  A_T = {A:.6f} MPa mm^{p:.4f}")

    # ---- 1. face slip
    F = a["F"] * A**n                    # full displacement angular function
    fpi = F[-1]
    for r in (0.02, 0.08):
        print(f"Dw(r={r}) = {2 * r**s * fpi:.6g} mm")
    print(f"   (F(pi) = {fpi:.6g} mm^(1-s); check ratio Dw.08/Dw.02 = "
          f"{(0.08/0.02)**s:.6f} = (4)^s)")

    # ---- 2. stresses at r=0.02, theta=60 deg
    r0, t0d = 0.02, np.pi / 3
    i60 = np.searchsorted(th, t0d)
    Tr60 = np.interp(t0d, th, a["Tr"]); Tt60 = np.interp(t0d, th, a["Tt"])
    trz = A * r0**(-p) * Tr60
    ttz = A * r0**(-p) * Tt60
    print(f"tau_rz(0.02,60deg) = {trz:.6g} MPa, tau_thz(0.02,60deg) = {ttz:.6g} MPa")

    # ---- 3. |gamma| = 0.01 contours
    gstar = 0.01
    for name, ang in [("0 deg", 0.0), ("90 deg", np.pi / 2), ("180 deg (face)", np.pi)]:
        mv = np.interp(ang, th, a["m"])
        gam_ang = (G0 / TAU0**n) * (A * mv)**n        # |gamma| = r^-np * this
        r = (gam_ang / gstar) ** (1.0 / (n * p))
        print(f"|gamma|=0.01 at theta={name}: r = {r:.6g} mm  (m={mv:.6g})")

    # ---- verification: J on Cartesian FD contours
    def w_xy(x, y):
        r_ = np.hypot(x, y)
        t_ = np.arctan2(y, x)              # in (-pi, pi]
        T_r_ = np.interp(np.abs(t_), th, a["Tr"]) * np.sign(t_) if False else None
        Tr_ = np.where(t_ >= 0,
                       np.interp(np.abs(t_), th, a["Tr"]),
                       -np.interp(np.abs(t_), th, a["Tr"]))
        return r_**s * A**n * (G0 / TAU0**n) * np.interp(np.abs(t_), th, a["m"])**(n - 1) * Tr_ / s

    def tau_xy(x, y):
        r_ = np.hypot(x, y); t_ = np.arctan2(y, x)
        at = np.abs(t_)
        Tr_ = np.interp(at, th, a["Tr"]) * np.where(t_ >= 0, 1, -1)
        Tt_ = np.interp(at, th, a["Tt"])
        tr = A * r_**(-p) * Tr_; tt = A * r_**(-p) * Tt_
        return tr * np.cos(t_) - tt * np.sin(t_), tr * np.sin(t_) + tt * np.cos(t_)

    def J_cart(R, k=200001, h=1e-5):
        t_ = np.linspace(-np.pi + 1e-9, np.pi - 1e-9, k)
        x, y = R * np.cos(t_), R * np.sin(t_)
        tx, ty = tau_xy(x, y)
        m_ = np.hypot(*tau_xy(x, y)) / (A * R**(-p))   # = m(theta)
        W = (n / (n + 1.0)) * (G0 / TAU0**n) * (A * m_)**(n + 1.0) * R**(-p * (n + 1.0))
        dwx = (w_xy(x + h, y) - w_xy(x - h, y)) / (2 * h)
        Nx, Ny = np.cos(t_), np.sin(t_)
        tz = tx * Nx + ty * Ny
        return np.trapezoid((W * Nx - tz * dwx) * R, t_)

    for R in (0.005, 0.02, 0.05, 0.15):
        print(f"J on circle R={R}: {J_cart(R):.6f} N/mm  (target 0.120)")

    # equilibrium FD check in Cartesian
    x, y = 0.03, 0.02
    h = 1e-6
    dtx = (tau_xy(x + h, y)[0] - tau_xy(x - h, y)[0]) / (2 * h)
    dty = (tau_xy(x, y + h)[1] - tau_xy(x, y - h)[1]) / (2 * h)
    print(f"equilibrium residual d(tx)/dx + d(ty)/dy at (0.03,0.02): {dtx + dty:.3e} MPa/mm "
          f"(terms ~{abs(dtx):.3f})")


if __name__ == "__main__":
    regression_n1()
    print("-" * 70)
    main()
