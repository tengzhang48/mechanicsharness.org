import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline

t0, g0 = 80.0, 0.002


def solve(n, q=None):
    if q is None:
        q = 1.0 / (n + 1.0)

    def rhs(t, y):
        F, Fp = y
        num = -q * (1 - q) * F * (q * q * F * F + Fp * Fp) \
              + (n - 1.0) / n * q * q * F * Fp * Fp
        den = q * q * F * F + Fp * Fp / n
        return [Fp, num / den]

    sol = solve_ivp(rhs, [0, np.pi], [0, 1], rtol=1e-13, atol=1e-15,
                    dense_output=True, max_step=2e-4)
    th = np.linspace(0, np.pi, 600001)
    Y = sol.sol(th)
    return q, th, Y[0], Y[1]


def build(n, J=0.120):
    q, th, F, Fp = solve(n)
    g = np.sqrt(q * q * F * F + Fp * Fp)
    thf = np.concatenate([-th[::-1], th[1:]])
    Ff = np.concatenate([-F[::-1], F[1:]])
    Fpf = np.concatenate([Fp[::-1], Fp[1:]])
    gf = np.sqrt(q * q * Ff * Ff + Fpf * Fpf)
    integ = (n / (n + 1.0) * gf ** (1 + 1.0 / n) * np.cos(thf)
             - q * Ff * gf ** (1.0 / n - 1)
             * (q * Ff * np.cos(thf) - Fpf * np.sin(thf)))
    Int = np.trapezoid(integ, thf)
    C = (J / (t0 * g0 ** (-1.0 / n) * Int)) ** (n / (n + 1.0))
    return dict(n=n, q=q, th=th, F=F, Fp=Fp, C=C, Int=Int,
                Fs=CubicSpline(th, F), Fps=CubicSpline(th, Fp))


def polar_tau(d, r, t):
    q, n, C = d['q'], d['n'], d['C']
    at = np.abs(t)
    F = np.where(t < 0, -1.0, 1.0) * d['Fs'](at)   # F odd
    Fp = d['Fps'](at)                                # F' even
    g = np.sqrt(q * q * F * F + Fp * Fp)
    A = t0 * g0 ** (-1.0 / n) * C ** (1.0 / n)
    pref = A * r ** (-q)
    tr = pref * q * F * g ** (1.0 / n - 1)
    tt = pref * Fp * g ** (1.0 / n - 1)
    return tr, tt, F, Fp, g, C


# ---- verification 1: J on circular contours, analytic gamma ----
print("Verification 1: J = R * int (W cos - tau_r * w_,x) dtheta (analytic)")
for n in (1.0, 2.0, 3.5, 6.0):
    d = build(n)
    q = d['q']
    for R in (0.01, 0.02, 0.05, 0.2):
        t = np.linspace(-np.pi, np.pi, 2000001)
        tr, tt, F, Fp, g, C = polar_tau(d, R, t)
        # gamma in polar
        gr = C * q * R ** (q - 1) * F
        gt = C * R ** (q - 1) * Fp
        W = n / (n + 1.0) * (tr * gr + tt * gt)
        w_x = C * R ** (q - 1) * (q * F * np.cos(t) - Fp * np.sin(t))
        Jnum = R * np.trapezoid(W * np.cos(t) - tr * w_x, t)
        print("  n=%.1f R=%.3f J=%.8f" % (n, R, Jnum))

# ---- verification 2: equilibrium residual ----
print("Verification 2: polar equilibrium residual (max over r,theta), n=3.5")
d = build(3.5)
for R in (0.02, 0.05):
    t = np.linspace(-np.pi + 1e-3, np.pi - 1e-3, 4001)
    h = 1e-6 * R
    trp = polar_tau(d, R + h, t)[0]
    trm = polar_tau(d, R - h, t)[0]
    dtr = (trp - trm) / (2 * h)
    # d tau_theta/dtheta
    tt_p = polar_tau(d, R, t + 1e-6)[1]
    tt_m = polar_tau(d, R, t - 1e-6)[1]
    dtt = (tt_p - tt_m) / 2e-6
    tr, tt, *_ = polar_tau(d, R, t)
    res = dtr + dtt / R + tr / R
    print("  R=%.3f max|res|/scale = %.3e" % (R, np.max(np.abs(res)) / np.max(np.abs(tr))))

# ---- verification 3: constitutive law consistency ----
print("Verification 3: gamma = (g0/t0)(tau/t0)^(n-1) tau, n=3.5")
n = 3.5
d = build(n)
q, C = d['q'], d['C']
t = np.linspace(-np.pi, np.pi, 4001)
R = 0.03
tr, tt, F, Fp, g, _ = polar_tau(d, R, t)
tau_mag = np.hypot(tr, tt)
gamma_pred_mag = g0 * (tau_mag / t0) ** n
gr = C * q * R ** (q - 1) * F
gt = C * R ** (q - 1) * Fp
gamma_mag = np.hypot(gr, gt)
print("  max rel err in |gamma| = %.3e" % np.max(np.abs(gamma_mag - gamma_pred_mag) / gamma_pred_mag))
# direction
print("  max angle err (rad) = %.3e" % np.max(np.abs(np.arctan2(gt, gr) - np.arctan2(tt, tr))))
