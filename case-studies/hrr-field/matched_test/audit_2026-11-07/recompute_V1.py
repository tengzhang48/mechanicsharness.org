"""Independent recomputation of the V1 key (antiplane power-law crack), audit 2026-11-07.

Written from V1/TASK_PROBLEM.md only. Differs from the builder's solve_reference.py (and from the
case page's hrr_check.py, which repeats the builder's route A) in the steps at risk:
  * formulation: Airy-type stress function (equilibrium exact), COMPATIBILITY as the governing
    equation, unknown g(theta) of phi = A r^mu g(theta); the builder uses the displacement w and
    equilibrium with the exponent fixed in advance;
  * exponent: mu found by SHOOTING on the face condition g(pi) = 0, not assumed; the radial
    exponent of w, lam = n(mu-1)+1, is then compared with 1/(n+1);
  * work density: W = int_0^|gamma| tau(s) ds by numerical quadrature of the stated law;
  * J: evaluated on a SQUARE contour and on a circle, in Cartesian components, with w_,x from
    central finite differences of the reconstructed displacement w(x, y); no polar J formula.
Run: python -I recompute_V1.py   (writes recompute_V1_results.json next to this file)
"""
import json
import os
import numpy as np
from scipy.integrate import solve_ivp, quad
from scipy.optimize import brentq

HERE = os.path.dirname(os.path.abspath(__file__))
n, tau0, gam0, J = 3.5, 80.0, 0.002, 0.120
c = gam0 / tau0**n                      # gamma = c |tau|^(n-1) tau   (engineering shear)


def tau_of_gamma(gmag):                 # inverse law, scalar
    return (gmag / c) ** (1.0 / n)


def solve_angular(mu, dense=False):
    """phi = r^mu g(theta); tau_r = r^(mu-1) g', tau_th = -mu r^(mu-1) g.
    Compatibility d(r gam_th)/dr = d gam_r/dtheta  ->  h' = -mu*lam*Q^(n-1) g,
    with h = Q^(n-1) g', Q^2 = g'^2 + mu^2 g^2, lam = n(mu-1)+1.  Symmetric (tau_th even):
    g'(0) = 0 -> h(0) = 0; normalization g(0) = -1 (tau_y > 0 ahead of the tip)."""
    lam = n * (mu - 1.0) + 1.0

    def gp_of(g, h):                    # invert h = (gp^2 + mu^2 g^2)^((n-1)/2) gp, monotone in gp
        if h == 0.0:
            return 0.0
        F = lambda x: (x * x + mu * mu * g * g) ** ((n - 1) / 2) * x - abs(h)
        hi = max(abs(h) ** (1.0 / n), 1e-300)
        while F(hi) < 0:
            hi *= 2.0
        return np.copysign(brentq(F, 0.0, hi, xtol=1e-300, rtol=1e-15, maxiter=500), h)

    def rhs(th, y):
        g, h = y
        gp = gp_of(g, h)
        Q = np.hypot(gp, mu * g)
        return [gp, -mu * lam * Q ** (n - 1) * g]

    sol = solve_ivp(rhs, (0.0, np.pi), [-1.0, 0.0], method="DOP853", rtol=1e-13, atol=1e-15,
                    dense_output=dense)
    return sol, gp_of, lam


# 1. shoot on mu: g(pi; mu) = 0 ----------------------------------------------------------------
mus = np.linspace(0.55, 0.98, 44)
gpi = [solve_angular(m)[0].y[0, -1] for m in mus]
brackets = [(mus[i], mus[i + 1]) for i in range(len(mus) - 1) if gpi[i] * gpi[i + 1] < 0]
mu = brentq(lambda m: solve_angular(m)[0].y[0, -1], *brackets[0], xtol=1e-15, rtol=1e-15)
sol, gp_of, lam = solve_angular(mu, dense=True)


def ang(th):
    """g, g', h at |th| <= pi with symmetry: g, g' even/odd; h odd."""
    s = np.sign(th) if th != 0 else 1.0
    g, h = sol.sol(abs(th))
    return g, s * gp_of(g, h), s * h


# 2. fields for amplitude A (phi = A r^mu g) ---------------------------------------------------
def w_field(x, y, A):
    r = np.hypot(x, y); th = np.arctan2(y, x)
    g, gp, h = ang(th)
    return c * A**n * r**lam * h / lam      # from w_,r = gamma_r (no rigid translation)


def grad_w_fd(x, y, A, d=1e-6):
    s = d * np.hypot(x, y)
    return np.array([(w_field(x + s, y, A) - w_field(x - s, y, A)) / (2 * s),
                     (w_field(x, y + s, A) - w_field(x, y - s, A)) / (2 * s)])


def W_quad(gmag):                        # stored work by quadrature of the stated law
    return quad(tau_of_gamma, 0.0, gmag, epsabs=0, epsrel=1e-13, limit=200)[0]


def J_integrand(x, y, N, A, Wfun=W_quad, drop=None):
    gr = grad_w_fd(x, y, A)
    gm = np.hypot(*gr)
    tau = tau_of_gamma(gm) * gr / gm
    Wv = Wfun(gm)
    t1 = Wv * N[0]
    t2 = (tau @ N) * gr[0]
    if drop == "traction":
        return t1
    if drop == "W":
        return -t2
    return t1 - t2


def J_square(A, h=1.0, **kw):
    """Counter-clockwise from the lower face (-h, 0-) to the upper face (-h, 0+)."""
    eps = 1e-12
    segs = [  # (param range, point(t), outward normal)
        ((-eps, -h), lambda t: (-h, t), (-1.0, 0.0)),   # left, lower: y from 0- down to -h
        ((-h, h), lambda t: (t, -h), (0.0, -1.0)),       # bottom
        ((-h, h), lambda t: (h, t), (1.0, 0.0)),         # right
        ((h, -h), lambda t: (t, h), (0.0, 1.0)),         # top, x from h to -h
        ((h, eps), lambda t: (-h, t), (-1.0, 0.0)),      # left, upper: y from h down to 0+
    ]
    tot = 0.0
    for (a, b), P, N in segs:          # ds is positive arc length: integrate over [min, max]
        f = lambda t: J_integrand(*P(t), np.array(N), A, **kw)
        tot += quad(f, min(a, b), max(a, b), epsabs=1e-13, epsrel=1e-11, limit=400)[0]
    return tot


def J_circle(A, r=1.0, **kw):
    f = lambda th: r * J_integrand(r * np.cos(th), r * np.sin(th),
                                   np.array([np.cos(th), np.sin(th)]), A, **kw)
    return quad(f, -np.pi + 1e-12, np.pi - 1e-12, epsabs=1e-13, epsrel=1e-11, limit=400)[0]


def outputs(A, cfac=1.0):
    """Requested numbers for amplitude A (cfac scales the compliance for a mutant)."""
    cc = c * cfac
    g, gp, h = ang(np.pi)
    slip = lambda r: 2 * cc * A**n * r**lam * h / lam
    g60, gp60, h60 = ang(np.pi / 3)
    r1 = 0.02
    tau_rz = A * r1 ** (mu - 1) * gp60
    tau_tz = -mu * A * r1 ** (mu - 1) * g60

    def rcont(th, gtarget=0.01):
        g_, gp_, h_ = ang(th)
        Q = np.hypot(gp_, mu * g_)
        # |gamma| = cc (A r^(mu-1) Q)^n = gtarget
        return ((gtarget / cc) ** (1 / n) / (A * Q)) ** (1 / (mu - 1))
    return {"slip_0.02": slip(0.02), "slip_0.08": slip(0.08), "tau_rz": tau_rz, "tau_tz": tau_tz,
            "r_0": rcont(0.0), "r_90": rcont(np.pi / 2), "r_180": rcont(np.pi)}


if __name__ == "__main__":
    out = {}
    out["shooting"] = {"mu_brackets": [list(map(float, b)) for b in brackets], "mu": mu,
                       "mu_expected_n_over_n+1": n / (n + 1), "lam_w_exponent": lam,
                       "lam_expected_1_over_n+1": 1 / (n + 1)}
    # W check: quadrature vs n/(n+1) tau*gamma at one point
    gt = 0.0137
    out["W_check"] = {"W_quad": W_quad(gt), "n_over_n+1_tau_gamma": n / (n + 1) * tau_of_gamma(gt) * gt}
    # J at A = 1 on two contours; J scales as A^(n+1)
    J1_sq = J_square(1.0, h=1.0)
    J1_sq2 = J_square(1.0, h=0.37)
    J1_ci = J_circle(1.0, r=0.6)
    out["J_unit_amplitude"] = {"square_h1": J1_sq, "square_h0.37": J1_sq2, "circle_r0.6": J1_ci}
    A = (J / J1_sq) ** (1 / (n + 1))
    vals = outputs(A)
    out["A_stress_function_amplitude"] = A
    out["values"] = vals
    key = {"slip_0.02": 0.001318126777760562, "slip_0.08": 0.0017936937381386464,
           "tau_rz": 38.586596862313044, "tau_tz": 138.01061663235828, "r_0": 0.046894423054945324,
           "r_90": 0.025066123488862637, "r_180": 0.013398406587127234}   # reference_output.json route A
    out["rel_diff_vs_key"] = {k: vals[k] / key[k] - 1 for k in key}
    # J of the final field on a third contour (closure check)
    out["J_final_square_h0.05"] = J_square(A, h=0.05)
    # face traction residual and n = 1 limit
    g, gp, h = ang(np.pi)
    out["face_tau_theta_over_max"] = abs(mu * g) / max(abs(np.hypot(*[ang(t)[1], mu * ang(t)[0]])) for t in np.linspace(0, np.pi, 50))
    json.dump(out, open(os.path.join(HERE, "recompute_V1_results.json"), "w"), indent=1)
    print(json.dumps(out, indent=1))
