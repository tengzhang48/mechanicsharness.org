#!/usr/bin/env python3
"""
Mode III (anti-plane shear) asymptotic field at a stationary slit tip in a
pure power-law material.

Constitutive law (given):  gamma_vec = (gamma0/tau0^n) * tau^(n-1) * tau_vec
Inverse:                   tau_vec  = tau0 * (gamma/gamma0)^(1/n) * gamma_vec/gamma
Equilibrium:  div tau = 0        Compatibility:  curl gamma = 0,  gamma = grad w
Stress function:  tau_r = (1/r) psi_theta,  tau_th = -psi_r,  psi = B r^(1-s) Phi(th)
Eigen ODE (s from J-dominance, s = 1/(n+1)):
    d/dth[ T^(n-1) Phi' ] + (1-s*n)(1-s) Phi T^(n-1) = 0,
    T^2 = Phi'^2 + (1-s)^2 Phi^2,   Phi(0) = -1, Phi'(0) = 0, Phi(pi) = 0
(displacement chosen so upper face moves +z).
Displacement:  w(r,th) = (c/s) B^n r^s G(th),  G = T^(n-1) Phi',  c = gamma0/tau0^n.
Amplitude B fixed by J = integral over circle (radius-independent).
Units: mm, N, MPa.
"""
import numpy as np
from scipy.integrate import solve_ivp, quad
from scipy.optimize import brentq

TAU0, GAMMA0, N, JAPP = 80.0, 0.002, 3.5, 0.120


def shoot(n, s):
    k = 1.0 - s
    mu = (1.0 - s * n) * k

    def ode(th, y):
        Phi, P = y
        T2 = P * P + k * k * Phi * Phi
        return [P, -Phi * (mu * T2 + (n - 1.0) * k * k * P * P) / (T2 + (n - 1.0) * P * P)]

    return solve_ivp(ode, (0.0, np.pi), (-1.0, 0.0), method='DOP853',
                     rtol=1e-12, atol=1e-13, dense_output=True, max_step=np.pi / 200.0)


def phi_pi(n, s):
    return shoot(n, s).y[0, -1]


def solve(n):
    s0 = 1.0 / (n + 1.0)
    sol = shoot(n, s0)
    res0 = sol.y[0, -1]
    if abs(res0) > 1e-8:  # verify s=1/(n+1) is the BVP eigenvalue
        sstar = brentq(lambda s: phi_pi(n, s), 0.7 / (n + 1.0), 1.5 / (n + 1.0), xtol=1e-14)
    else:
        sstar = s0
    s = sstar
    sol = shoot(n, s)
    k = 1.0 - s
    c = GAMMA0 / TAU0 ** n

    # solution defined on [0, pi]; extended by symmetry: Phi even, P odd, T even
    def Phi(th):
        return sol.sol(np.abs(th))[0]

    def P(th):
        return np.sign(th) * sol.sol(np.abs(th))[1]

    def T(th):
        return np.sqrt(sol.sol(np.abs(th))[1] ** 2 + k * k * sol.sol(np.abs(th))[0] ** 2)

    # J integral:  J = (gamma0/tau0^n) B^(n+1) * I , integrand even in theta
    def Ifun(th):
        p, ph = P(th), Phi(th)
        t = np.sqrt(p * p + k * k * ph * ph)
        return t ** (n - 1) * ((n / (n + 1.0)) * t * t * np.cos(th)
                               - p * p * np.cos(th) - k * ph * p * np.sin(th))

    Ival = 2.0 * quad(Ifun, 0.0, np.pi, limit=400)[0]
    B = (JAPP * TAU0 ** n / (GAMMA0 * Ival)) ** (1.0 / (n + 1.0))
    return dict(n=n, s=s, k=k, c=c, B=B, I=Ival, sol=sol, res=res0, sstar=sstar,
                Phi=Phi, P=P, T=T)


def w_field(r, th, f):
    G = f['T'](th) ** (f['n'] - 1) * f['P'](th)
    return (f['c'] / f['s']) * f['B'] ** f['n'] * r ** f['s'] * G


def tau_components(r, th, f):
    tr = f['B'] * r ** (-f['s']) * f['P'](th)
    tt = -f['k'] * f['B'] * r ** (-f['s']) * f['Phi'](th)
    return tr, tt


def J_direct(r, f):
    """J from direct contour integral over circle of radius r."""
    n, s, k, B, c = f['n'], f['s'], f['k'], f['B'], f['c']

    def integrand(th):
        Phi = f['Phi'](th)
        P = f['P'](th)
        T = f['T'](th)
        tau_r = B * r ** (-s) * P
        gam_r = c * B ** n * r ** (-n * s) * T ** (n - 1) * P
        gam_th = -k * c * B ** n * r ** (-n * s) * T ** (n - 1) * Phi
        gx = gam_r * np.cos(th) - gam_th * np.sin(th)
        Wd = (n / (n + 1.0)) * TAU0 * GAMMA0 * (B * r ** (-s) * T / TAU0) ** (n + 1.0)
        return (Wd * np.cos(th) - tau_r * gx) * r

    return quad(integrand, -np.pi, np.pi, limit=400)[0]


def report(tag, f):
    n, s, k, B = f['n'], f['s'], f['k'], f['B']
    print(f"\n=== {tag}: n={n}, s={s:.12f}, eigen-s*={f['sstar']:.12f}, "
          f"Phi(pi) residual={f['sol'].y[0,-1]:.3e}, I={f['I']:.8f}, B={B:.8f}")
    if abs(n - 3.5) < 1e-12:
        print("-- requested quantities --")
        for r in (0.0200, 0.0800):
            print(f"  Dw(r={r:.4f}) = {w_field(r, np.pi, f) - w_field(r, -np.pi, f):.6e} mm"
                  f"   (w(r,pi)={w_field(r, np.pi, f):.4e})")
        th = np.deg2rad(60.0)
        tr, tt = tau_components(0.0200, th, f)
        print(f"  tau_rz(0.02,60deg) = {tr:.6f} MPa,  tau_thz(0.02,60deg) = {tt:.6f} MPa")
        gstar = 0.0100
        for thd in (0.0, 90.0, 180.0):
            th = np.deg2rad(thd)
            Tval = f['T'](th)
            rstar = (B * Tval / (TAU0 * (gstar / GAMMA0) ** (1.0 / n))) ** (1.0 / s)
            print(f"  r*(|gamma|=0.01, theta={thd:5.1f} deg) = {rstar:.6f} mm  (T={Tval:.6f})")
        print("-- verifications --")
        for r in (0.0200, 0.0800, 0.2000):
            print(f"  J_direct(r={r:.3f}) = {J_direct(r, f):.9f} N/mm  (target {JAPP})")
        ttp = -k * B * 0.0200 ** (-s) * f['Phi'](np.pi)
        print(f"  face traction residual tau_thz(r=0.02, theta=pi) = {ttp:.3e} MPa")
        print(f"  check w(r,-pi)=-w(r,pi): {w_field(0.02,-np.pi,f):.4e}")


def elastic_check():
    """n -> 1 limit vs exact LEFM K-field, mu = tau0/gamma0, K = sqrt(2 mu J)."""
    f = solve(1.0)
    mu = TAU0 / GAMMA0
    K = np.sqrt(2.0 * mu * JAPP)
    print("\n=== elastic limit check (n=1), mu=%.1f MPa, K=%.6f MPa*sqrt(mm)" % (mu, K))
    th = np.linspace(0.0, np.pi, 21)
    print("  max |Phi + cos(th/2)| =", np.max(np.abs(f['sol'].sol(th)[0] + np.cos(th / 2))))
    for r in (0.02, 0.08):
        dw_num = w_field(r, np.pi, f) - w_field(r, -np.pi, f)
        dw_ex = 4.0 * np.sqrt(JAPP * r / (np.pi * mu))
        print(f"  Dw({r}): num={dw_num:.8e} exact={dw_ex:.8e} rel={abs(dw_num/dw_ex-1):.2e}")
    th = np.deg2rad(60.0)
    trn, ttn = tau_components(0.02, th, f)
    tre = K * np.sin(th / 2) / np.sqrt(2 * np.pi * 0.02)
    tte = K * np.cos(th / 2) / np.sqrt(2 * np.pi * 0.02)
    print(f"  tau_rz: num={trn:.6f} exact={tre:.6f} rel={abs(trn/tre-1):.2e}")
    print(f"  tau_thz: num={ttn:.6f} exact={tte:.6f} rel={abs(ttn/tte-1):.2e}")
    r_ex = K ** 2 / (2 * np.pi * (0.01 * mu) ** 2)  # |tau|=0.01*mu -> r (theta-indep)
    for thd in (0.0, 90.0, 180.0):
        Tval = f['T'](np.deg2rad(thd))
        r_num = (f['B'] * Tval / (TAU0 * (0.01 / GAMMA0))) ** 2
        print(f"  r*(theta={thd:5.1f}): num={r_num:.8e} exact={r_ex:.8e} rel={abs(r_num/r_ex-1):.2e}")


if __name__ == "__main__":
    f35 = solve(3.5)
    report("power law n=3.5", f35)
    elastic_check()
