"""Mode III (anti-plane) crack-tip field for a pure power-law material.

Model: separable leading field  w(r,theta) = a * r^p * h(theta),  p = 1/(n+1)
(stress ~ r^(-1/(n+1)), strain ~ r^(-n/(n+1)); exponent fixed by finite, nonzero,
path-independent J). Constitutive law: gamma_vec = (gamma0/tau0)(tau/tau0)^(n-1) tau_vec,
inverted: tau_vec = tau0 (|gamma|/gamma0)^(1/n) gamma_vec/|gamma|.

Equilibrium div(tau)=0 gives the 2nd-order angular ODE
    d/dtheta[ q^nu h' ] + beta q^nu h = 0,
    q = sqrt(p^2 h^2 + h'^2),  nu = (1-n)/n,  beta = n/(n+1)^2.
BCs: h(0)=0 (antisymmetry), h'(pi)=0 (traction-free face).

Eshelby conservation behind J gives the pointwise first integral
    W sin(theta) + tau_theta_z w_,x = 0   (constant 0 from the free face),
i.e.  sin(t) h'^2 - cos(t) h h' - beta sin(t) h^2 = 0, so
    h' = h [cos + sqrt(cos^2 + 4 beta sin^2)] / (2 sin)   (+ branch, h' > 0),
which implies both BCs. Normalize h ~ theta as theta -> 0 (h'(0)=1).

Amplitude: J = K_H * I_J, K_H = tau0*gamma0*(a/gamma0)^((n+1)/n),
    I_J = int_{-pi}^{pi} [ (n/(n+1)) q^((n+1)/n) cos(t)
                           - q^nu p h (p h cos(t) - h' sin(t)) ] dt.
"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline

TAU0 = 80.0      # MPa = N/mm^2
GAM0 = 0.002     # -
N = 3.5          # -
J_APPLIED = 0.120  # N/mm


def angular_solution(n, ntheta=200001):
    """Integrate the first-integral ODE; return theta, h, h', q on [0, pi], h ~ theta at 0."""
    p = 1.0 / (n + 1.0)
    beta = n / (n + 1.0) ** 2

    def u(t):
        s, c = np.sin(t), np.cos(t)
        return (c + np.sqrt(c * c + 4.0 * beta * s * s)) / (2.0 * s)

    th0 = 1e-9
    sol = solve_ivp(lambda t, y: [y[0] * u(t)], (th0, np.pi - 1e-12), [th0],
                    method="DOP853", dense_output=True, rtol=1e-12, atol=1e-16)
    assert sol.success, sol.message
    tg = np.linspace(0.0, np.pi, ntheta)
    hg = tg.copy()                      # h ~ theta below th0
    hg[tg >= th0] = sol.sol(tg[tg >= th0])[0]
    hpg = np.ones_like(tg)              # normalization h'(0) = 1
    pos = tg > 0.0
    hpg[pos] = hg[pos] * u(tg[pos])
    q = np.sqrt(p * p * hg * hg + hpg * hpg)
    return tg, hg, hpg, q


def J_angular_integral(tg, hg, hpg, q, n):
    """I_J = J / K_H with the even extension h(-t) = -h(t); full-range quadrature."""
    p = 1.0 / (n + 1.0)
    nu = (1.0 - n) / n
    tF = np.concatenate([-tg[:0:-1], tg])
    hF = np.concatenate([-hg[:0:-1], hg])
    hpF = np.concatenate([hpg[:0:-1], hpg])
    qF = np.concatenate([q[:0:-1], q])
    A = ((n / (n + 1.0)) * qF ** ((n + 1.0) / n) * np.cos(tF)
         - qF ** nu * p * hF * (p * hF * np.cos(tF) - hpF * np.sin(tF)))
    return np.trapezoid(A, tF)


def solve_field(n, J):
    tg, hg, hpg, q = angular_solution(n)
    I = J_angular_integral(tg, hg, hpg, q, n)
    K_H = J / I                                  # = tau0*gam0*(a/gam0)^((n+1)/n)
    a_over_g0 = (K_H / (TAU0 * GAM0)) ** (n / (n + 1.0))
    a = GAM0 * a_over_g0                          # units mm^(n/(n+1))
    return dict(tg=tg, hg=hg, hpg=hpg, q=q, I=I, K_H=K_H, a=a, n=n, J=J)


def interp(F, t):
    return (np.interp(t, F["tg"], F["hg"]), np.interp(t, F["tg"], F["hpg"]),
            np.interp(t, F["tg"], F["q"]))


def stresses(F, r, t):
    n = F["n"]; p = 1.0 / (n + 1.0); nu = (1.0 - n) / n
    h, hp, q = interp(F, t)
    S = TAU0 * (F["a"] / GAM0) ** (1.0 / n) * r ** (-1.0 / (n + 1.0))
    return S * q ** nu * p * h, S * q ** nu * hp     # tau_rz, tau_theta_z


def slip(F, r):
    n = F["n"]; p = 1.0 / (n + 1.0)
    return 2.0 * F["a"] * r ** p * F["hg"][-1]


def strain_radius(F, t, gam_target=0.0100):
    n = F["n"]
    _, _, q = interp(F, t)
    return (F["a"] * q / gam_target) ** ((n + 1.0) / n)


# ---------------------------------------------------------------------- checks
def check_n1_limit():
    """n=1 must reduce to the linear K_III field. With normalization h'(0)=1:
    h = 2 sin(t/2); I_J = pi (integrand ~ h^2, so 4x the h=sin(t/2) value pi/4);
    a = (K/mu)/sqrt(2 pi) from w = a r^1/2 2 sin(t/2), J = K^2/(2 mu), mu = tau0/gam0."""
    F = solve_field(1.0, J_APPLIED)
    h_ref = 2.0 * np.sin(F["tg"] / 2.0)
    dh = np.max(np.abs(F["hg"] - h_ref)) / np.max(np.abs(h_ref))
    dI = abs(F["I"] - np.pi) / np.pi
    mu = TAU0 / GAM0
    a_ref = (np.sqrt(2.0 * mu * J_APPLIED) / mu) / np.sqrt(2.0 * np.pi)
    da = abs(F["a"] - a_ref) / a_ref
    print(f"[check n=1] max|h-2sin(t/2)|/max = {dh:.3e}; |I-pi/4|/(pi/4) = {dI:.3e}; "
          f"|a - K-field a|/a = {da:.3e}")
    assert dh < 1e-8 and dI < 1e-10 and da < 1e-8


def check_equilibrium_residual(F):
    """Residual of the 2nd-order equilibrium ODE d/dt[q^nu h'] + beta q^nu h = 0."""
    n = F["n"]; tg = F["tg"]
    nu = (1.0 - n) / n; beta = n / (n + 1.0) ** 2
    g = F["q"] ** nu * F["hpg"]
    gp = CubicSpline(tg, g).derivative()(tg)
    res = gp + beta * F["q"] ** nu * F["hg"]
    scale = np.max(np.abs(gp)) + np.max(np.abs(beta * F["q"] ** nu * F["hg"]))
    r = np.max(np.abs(res[200:-200])) / scale
    print(f"[check equilibrium] max |ODE residual|/scale (interior) = {r:.3e}")
    assert r < 1e-6


def check_faces(F):
    hp_face = abs(F["hpg"][-1])
    print(f"[check faces] h'(pi) = {hp_face:.3e} (traction-free face, should be ~0)")
    assert hp_face < 1e-8


def check_J_cartesian(F, R=0.05, k=20002):
    """Independent J: build w(x,y) on a circle, Cartesian FD for w_,x, w_,y, then
    stresses from the vector law and W = n/(n+1) tau.gamma; compare to applied J.
    Midpoint grid: no node on the slit faces (w jumps there); FD step << node gap."""
    n = F["n"]; p = 1.0 / (n + 1.0); a = F["a"]
    dt = 2.0 * np.pi / k
    t = -np.pi + (np.arange(k) + 0.5) * dt
    x, y = R * np.cos(t), R * np.sin(t)

    def wfield(xx, yy):
        rr = np.hypot(xx, yy)
        tt = np.arctan2(yy, xx)
        return a * rr ** p * (np.sign(tt) * np.interp(np.abs(tt), F["tg"], F["hg"]))

    d = 1e-6 * R   # << node-face gap (pi/(k-1) ~ 1.6e-4 rad), FD error ~(d/R)^2 = 1e-12
    wx = (wfield(x + d, y) - wfield(x - d, y)) / (2.0 * d)
    wy = (wfield(x, y + d) - wfield(x, y - d)) / (2.0 * d)
    gm = np.hypot(wx, wy)
    taum = TAU0 * (gm / GAM0) ** (1.0 / n)
    txz, tyz = taum * wx / gm, taum * wy / gm
    W = (n / (n + 1.0)) * (txz * wx + tyz * wy)
    tz = txz * np.cos(t) + tyz * np.sin(t)
    Jnum = np.sum((W * np.cos(t) - tz * wx) * R) * dt
    rel = abs(Jnum - F["J"]) / F["J"]
    print(f"[check J cartesian FD] J(R={R}) = {Jnum:.6f} N/mm vs applied {F['J']:.6f}; rel = {rel:.3e}")
    return rel


def main():
    check_n1_limit()
    F = solve_field(N, J_APPLIED)
    p = 1.0 / (N + 1.0)
    print(f"\nn = {N}, p = 1/(n+1) = {p:.6f}, beta = {N/(N+1)**2:.6f}")
    print(f"I_J = {F['I']:.8f};  h(pi) = {F['hg'][-1]:.8f};  a = {F['a']:.6e} mm^(n/(n+1))")
    check_equilibrium_residual(F)
    check_faces(F)
    for R in (0.02, 0.05, 0.10):
        check_J_cartesian(F, R=R)

    print("\n=== Requested numbers ===")
    for r in (0.0200, 0.0800):
        print(f"1) Delta w(r={r} mm) = {slip(F, r):.6e} mm")
    tr, tt = stresses(F, 0.0200, np.pi / 3.0)
    print(f"2) tau_rz(0.02, 60deg) = {tr:.6f} MPa ; tau_theta_z = {tt:.6f} MPa")
    for name, t in (("0 deg", 0.0), ("90 deg", np.pi / 2.0), ("180 deg (upper face)", np.pi)):
        print(f"3) |gamma|=0.0100 at theta={name}: r = {strain_radius(F, t):.6e} mm")


if __name__ == "__main__":
    main()

