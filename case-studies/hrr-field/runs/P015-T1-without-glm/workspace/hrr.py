"""
HRR plane-strain crack-tip field, n = 6.5.
sigma from Airy Phi = r^(2-p) phi(theta):
  sig_rr = (2-p)phi + phi'', sig_tt = (2-p)(1-p)phi, sig_rt = -(1-p)phi'
plane strain: sig_zz=(sig_rr+sig_tt)/2, s_rr=(sig_rr-sig_tt)/2=-s_tt, s_rt=sig_rt,
  sigma_e = sqrt(3(s_rr^2+s_rt^2)).
separable eps ~ r^(-pn) with p(n+1)=1; kinematic compatibility for u_r=r^p a, u_th=r^p b:
  e~_rr = p a, e~_tt = a+b', 2 e~_rt = a' + (p-1)b  and e~_tt = -e~_rr (incompressible)
  -> ODE: (e~_rr)'' + (1-p^2) e~_rr - 2p (e~_rt)' = 0,
     e~_rr = (3/2) se^(n-1) s_rr, e~_rt = (3/2) se^(n-1) s_rt.
BVP on [0,pi], p free eigenvalue: phi'(0)=0, phi'''(0)=0 (symmetry),
  phi(pi)=phi'(pi)=0 (traction-free), phi(0)=1 (scale).
J = alpha e0 sigma0 I_n (B/sigma0)^(n+1); sigma = sigma0 (J/(alpha e0 sigma0 I_n r))^p sig~
H = sig~_tt(0)/I_n^p ;  C = U(R)/(J R) = (n/(n+1)) * int sig~_e^(n+1) dtheta / I_n.
"""
import numpy as np
from scipy.integrate import solve_bvp
from scipy.interpolate import CubicSpline

PI = np.pi
n_t = 6.5
alpha, E, s0 = 0.02, 200e9, 400e6
e0 = s0 / E

def F_phi4(y1, y2, y3, y4, n, p):
    """phi'''' from compatibility ODE, linear solve in phi''''."""
    s1 = (y3 + p * (2.0 - p) * y1) / 2.0
    s2 = -(1.0 - p) * y2
    R = np.sqrt(s1 * s1 + s2 * s2)
    se = np.sqrt(3.0) * R
    s1p = (y4 + p * (2.0 - p) * y2) / 2.0
    s2p = -(1.0 - p) * y3
    Rp = (s1 * s1p + s2 * s2p) / R
    AA = (s1p * s1p + s2p * s2p - Rp * Rp + s2 * (-(1.0 - p) * y4)) / R  # R'' w/o s1'' term
    BB = s1 / (2.0 * R)                                                   # coeff of phi'''' in R''
    wcoef = (n - 1.0) * se ** (n - 2.0) * np.sqrt(3.0) * BB * s1 + 0.5 * se ** (n - 1.0)
    c0 = ((n - 1.0) * (n - 2.0) * se ** (n - 3.0) * 3.0 * Rp * Rp * s1
          + (n - 1.0) * se ** (n - 2.0) * np.sqrt(3.0) * AA * s1
          + 2.0 * (n - 1.0) * se ** (n - 2.0) * np.sqrt(3.0) * Rp * s1p
          + se ** (n - 1.0) * (p * (2.0 - p) * y3 / 2.0)
          + (1.0 - p * p) * se ** (n - 1.0) * s1
          - 2.0 * p * ((n - 1.0) * se ** (n - 2.0) * np.sqrt(3.0) * Rp * s2
                       + se ** (n - 1.0) * s2p))
    return -c0 / wcoef

def make_bvp(nv):
    def ode(th, y, pp):
        p = pp[0]
        return np.vstack([y[1], y[2], y[3], F_phi4(y[0], y[1], y[2], y[3], nv, p)])
    def bc(ya, yb, pp):
        return np.array([ya[1], ya[3], yb[0], yb[1], ya[0] - 1.0])
    return ode, bc

# ---- guess: elastic-like shape + perturbation so phi''(pi) != 0 (avoids sigma_e = 0) ----
def guess(x, A=0.15):
    y1 = np.cos(x / 2) ** 3 + A * (1.0 + np.cos(x))
    y1 /= y1[0]
    y2 = -1.5 * np.cos(x / 2) ** 2 * np.sin(x / 2) - A * np.sin(x)
    y3 = 1.5 * np.cos(x / 2) * np.sin(x / 2) ** 2 - 0.75 * np.cos(x / 2) ** 3 + A * np.cos(x)
    y4 = 1.5 * (np.cos(x / 2) ** 2 * np.sin(x / 2) - 0.5 * np.sin(x / 2) ** 3) \
         + 1.125 * np.cos(x / 2) ** 2 * np.sin(x / 2) - A * np.sin(x)
    s = 1.0 / (np.cos(0) ** 3 + A * 2.0)
    return s * np.vstack([y1, y2, y3, y4])

xs = np.linspace(0.0, PI, 241)
yg = guess(xs)
sol = None
for tol in [1e-6, 1e-9, 1e-11]:
    ode, bc = make_bvp(6.5)
    s2_ = solve_bvp(ode, bc, xs, yg, p=[sol.p[0] if sol else 2.0 / 7.5],
                    tol=tol, max_nodes=100000)
    print(f"tol={tol:.0e} status={s2_.status} p={s2_.p[0]:.12f} "
          f"maxres={s2_.rms_residuals.max():.2e} nodes={s2_.x.size}")
    if s2_.status != 0:
        break
    sol = s2_
    yg = sol.sol(xs)

if sol is None or sol.status != 0:
    raise SystemExit("bvp failed")

p_num = sol.p[0]
p_th = 1.0 / (n_t + 1)
print(f"eigenvalue p = {p_num:.12f}, theory 1/(n+1) = {p_th:.12f}, diff = {abs(p_num-p_th):.2e}")

# ---------------- angular functions ----------------
def angular(s, p, n=n_t, M=200001):
    th = np.linspace(0.0, PI, M)
    phi, d1, d2, d3 = s.sol(th)
    s1 = (d2 + p * (2 - p) * phi) / 2
    s2 = -(1 - p) * d1
    se = np.sqrt(3 * (s1 ** 2 + s2 ** 2))
    s1p = (d3 + p * (2 - p) * d1) / 2
    s2p = -(1 - p) * d2
    sep = 3 * (s1 * s1p + s2 * s2p) / se
    A = dict(th=th, se=se,
             sig_rr=(2 - p) * phi + d2, sig_tt=(2 - p) * (1 - p) * phi,
             sig_rt=-(1 - p) * d1,
             e_rr=1.5 * se ** (n - 1) * s1, e_rt=1.5 * se ** (n - 1) * s2)
    A['e_rrp'] = 1.5 * ((n - 1) * se ** (n - 2) * sep * s1 + se ** (n - 1) * s1p)
    A['a'] = A['e_rr'] / p
    A['ap'] = A['e_rrp'] / p
    A['b'] = (2 * A['e_rt'] - A['e_rrp']) / (p - 1)
    A['bp'] = -A['e_rr'] * (1.0 + 1.0 / p)
    return A

def full_range(A):
    th = A['th']
    out = {'TH': np.concatenate([-th[::-1], th[1:]])}
    odd = ['sig_rt', 'e_rt', 'b']
    for k in ['se', 'sig_rr', 'sig_tt', 'e_rr', 'e_rrp', 'a', 'ap'] + odd:
        f = A[k]
        g = np.concatenate([f[::-1], f[1:]])
        if k in odd:
            g = np.concatenate([-f[::-1], f[1:]])
        out[k] = g
    return out

def In_integral(F, n, p):
    TH = F['TH']
    fI = ((n / (n + 1)) * F['se'] ** (n + 1) * np.cos(TH)
          - F['sig_rr'] * (p * F['a'] * np.cos(TH) - np.sin(TH) * (F['ap'] - F['b']))
          - F['sig_rt'] * (p * F['b'] * np.cos(TH) - np.sin(TH) * (F['bp'] + F['a'])))
    return np.trapz(fI, TH), np.trapz(F['se'] ** (n + 1), TH)

def results(s, p, label):
    A = angular(s, p)
    F = full_range(A)
    In, Mint = In_integral(F, n_t, p)
    H = A['sig_tt'][0] / In ** p
    C = (n_t / (n_t + 1)) * Mint / In
    print(f"[{label}] p={p:.12f} I_n={In:.8f} M={Mint:.8f} H={H:.8f} C={C:.8f} "
          f"sig~tt(0)={A['sig_tt'][0]:.6f} min_se={A['se'].min():.4f} "
          f"theta_at_min={A['th'][np.argmin(A['se'])]:.3f}")
    return A, F, In, Mint, H, C

A, F, In, Mint, H, C = results(sol, p_th, "p=1/(n+1)")
_, _, _, _, H2, C2 = results(sol, p_num, "p=eigenvalue")   # sensitivity to p

# normalization invariance
lam = 1.0 / A['se'].max()
Hn = (lam * A['sig_tt'][0]) / (lam ** (n_t + 1) * In) ** p_th
print(f"normalization invariance: H'={Hn:.10f} vs H={H:.10f}")

# kinematic compatibility check: FD derivative of b vs analytic b'
db = np.gradient(F['b'], F['TH'])
print(f"compatibility: max|b'(FD)-b'p|/max|b'p| = "
      f"{np.max(np.abs(db - F['bp'])) / np.max(np.abs(F['bp'])):.3e}")

# ---------------- verification: direct numerical J on two circles ----------------
cs_a = CubicSpline(A['th'], A['a']); cs_b = CubicSpline(A['th'], A['b'])
cs_rr = CubicSpline(A['th'], A['sig_rr']); cs_rt = CubicSpline(A['th'], A['sig_rt'])
cs_se = CubicSpline(A['th'], A['se'])
aof = lambda t: cs_a(np.abs(t))
bof = lambda t: np.sign(t) * cs_b(np.abs(t))
rrof = lambda t: cs_rr(np.abs(t))
rtof = lambda t: np.sign(t) * cs_rt(np.abs(t))
seof = lambda t: cs_se(np.abs(t))

B = s0
Jst = alpha * e0 * s0 * In
def ufield(X, Y):
    r = np.sqrt(X * X + Y * Y); t = np.arctan2(Y, X)
    rp = r ** p_th
    ur = alpha * e0 * (B / s0) ** n_t * rp * aof(t)
    ut = alpha * e0 * (B / s0) ** n_t * rp * bof(t)
    return ur * np.cos(t) - ut * np.sin(t), ur * np.sin(t) + ut * np.cos(t)

for r0 in [0.3, 0.05]:
    h = 1e-6 * r0
    TT = np.linspace(-PI, PI, 40001)
    X, Y = r0 * np.cos(TT), r0 * np.sin(TT)
    _, uyp2 = ufield(X + h, Y); _, uym2 = ufield(X - h, Y)
    uxp2, _ = ufield(X + h, Y); uxm2, _ = ufield(X - h, Y)
    dux_dx = (uxp2 - uxm2) / (2 * h)
    duy_dx = (uyp2 - uym2) / (2 * h)
    tr = B * r0 ** (-p_th) * rrof(TT); tth = B * r0 ** (-p_th) * rtof(TT)
    tx = tr * np.cos(TT) - tth * np.sin(TT)
    ty = tr * np.sin(TT) + tth * np.cos(TT)
    Wr = (n_t / (n_t + 1)) * alpha * e0 * s0 * (B / s0) ** (n_t + 1) * seof(TT) ** (n_t + 1)
    Jn = np.trapz(Wr * np.cos(TT) - (tx * dux_dx + ty * duy_dx), TT) * r0
    print(f"direct J, r={r0}: J_num={Jn:.6e} J*={Jst:.6e} rel.err={Jn / Jst - 1:.3e}")

# ---------------- verification: n=1 elastic limit, I_1 = 3*pi/2 ----------------
M = 400001
TH = np.linspace(-PI, PI, M)
c, s = np.cos(TH / 2), np.sin(TH / 2)
sig_rr_e = c * (1 + s ** 2); sig_tt_e = c ** 3; sig_rt_e = s * c ** 2
s1_e = (sig_rr_e - sig_tt_e) / 2; s2_e = sig_rt_e
se_e = np.sqrt(3 * (s1_e ** 2 + s2_e ** 2))
e_rr_e = 1.5 * s1_e; e_rt_e = 1.5 * s2_e
e_rrp_e = np.gradient(e_rr_e, TH)
p1 = 0.5
a1 = e_rr_e / p1; ap1 = np.gradient(a1, TH)
b1 = (2 * e_rt_e - e_rrp_e) / (p1 - 1); bp1 = -e_rr_e * (1 + 1 / p1)
fI1 = (0.5 * se_e ** 2 * np.cos(TH)
       - sig_rr_e * (p1 * a1 * np.cos(TH) - np.sin(TH) * (ap1 - b1))
       - sig_rt_e * (p1 * b1 * np.cos(TH) - np.sin(TH) * (bp1 + a1)))
I1 = np.trapz(fI1, TH)
print(f"n=1 elastic: I_1={I1:.10f} vs 3pi/2={1.5 * PI:.10f}  "
      f"H_1={sig_tt_e[M // 2] / I1 ** 0.5:.8f} vs 1/sqrt(3pi/2)={1 / np.sqrt(1.5 * PI):.8f}")

print("\n================ FINAL ================")
print(f"p = {p_th:.10f} (eigenvalue solve: {p_num:.10f})")
print(f"H = {H:.6f}")
print(f"C = {C:.6f}")
semax = A['se'].max()
print(f"Hutchinson norm (max sig~_e=1): I_n={In / semax ** (n_t + 1):.6f}, "
      f"sig~_tt(0)={A['sig_tt'][0] / semax:.6f}, sig~_rr(0)={A['sig_rr'][0] / semax:.6f}, "
      f"theta(max se)={A['th'][np.argmax(A['se'])]:.4f}")
np.save("results.npy", np.array([p_th, H, C]))
