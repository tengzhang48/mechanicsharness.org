"""Solve the P015-T1 task: HRR near-tip field, plane strain, fully plastic law.

Compute:
  1. p   in sigma_tt(r,0) ~ r^{-p}
  2. H   in sigma_tt(r,0) = H sigma0 [J/(alpha eps0 sigma0 r)]^p
  3. C   = U(R)/(J R),  U(R) = int_0^R int_{-pi}^{pi} W r dtheta dr
"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

# ---- material ----
n = 6.5
alpha = 0.02
sig0 = 400e6          # Pa
E = 200e9             # Pa
eps0 = sig0 / E


def Mcoeff(s):
    """se_tilde^2 = z^T M z, z=(f,f',f''), plane strain."""
    M = np.zeros((3, 3))
    M[0, 0] = 0.75 * s * s * (2 - s) ** 2
    M[1, 1] = 3 * (1 - s) ** 2
    M[2, 2] = 0.75
    M[0, 2] = M[2, 0] = 0.75 * s * (2 - s)
    return M


def make_rhs(n, s):
    """f'''' as a function of (f,f',f'',f''') from the plane-strain ODE."""
    k = (n - 1.0) / 2.0
    a = n * (s - 2.0)
    b = a + 1.0
    M = Mcoeff(s)

    def rhs(x0, x1, x2, x3):
        z = np.array([x0, x1, x2])
        v = np.array([x1, x2, x3])
        Q = z @ M @ z
        dQ = 2.0 * (M @ z)
        dQv = dQ @ v
        vTMv = v @ M @ v
        rho = Q ** k
        g = k * Q ** (k - 1.0) * dQ
        g0, g1, g2 = g
        R0 = (k * Q ** (k - 1.0) * (2.0 * vTMv)
              + k * (k - 1.0) * Q ** (k - 2.0) * dQv ** 2
              + g0 * x2 + g1 * x3)
        P = g @ v
        C = s * (2.0 - s) * x0 + x2
        Cp = s * (2.0 - s) * x1 + x3
        coef = rho + g2 * C
        rest = (R0 * C + 2.0 * P * Cp + rho * s * (2.0 - s) * x2
                - a * (a + 2.0) * rho * C
                + 4.0 * (s - 1.0) * b * (P * x1 + rho * x2))
        return -rest / coef

    return rhs


def integrate(n, s, p):
    rhs = make_rhs(n, s)
    f = lambda t, y: [y[1], y[2], y[3], rhs(y[0], y[1], y[2], y[3])]
    return solve_ivp(f, [0.0, np.pi], [1.0, 0.0, p, 0.0],
                     rtol=1e-12, atol=1e-14, dense_output=True, max_step=0.002)


def shoot(n, s_guess=None):
    """Find (s, f''(0)) with f(pi)=0 and f'(pi)=0 (the dominant HRR mode).

    Strategy: at the analytic guess s0 = (2n+1)/(n+1) scan p = f''(0) on a
    wide grid, bracket the root of f(pi)=0 that lies nearest to the known
    physical value, then polish both s and p with a 2-parameter Newton.
    """
    if s_guess is None:
        s_guess = (2 * n + 1) / (n + 1)

    def fend(s, p, rtol=1e-10, mstep=0.005):
        rhs = make_rhs(n, s)
        ff = lambda t, y: [y[1], y[2], y[3], rhs(y[0], y[1], y[2], y[3])]
        sol = solve_ivp(ff, [0.0, np.pi], [1.0, 0.0, p, 0.0],
                        rtol=rtol, atol=1e-12, dense_output=False, max_step=mstep)
        return sol.y[0, -1], sol.y[1, -1]

    # bracket f(pi)=0 in p at the analytic s
    ps = np.linspace(-3.0, -0.05, 120)
    vals = [fend(s_guess, p)[0] for p in ps]
    roots = []
    for i in range(len(ps) - 1):
        if vals[i] * vals[i + 1] < 0:
            r = brentq(lambda p: fend(s_guess, p)[0], ps[i], ps[i + 1], xtol=1e-12)
            roots.append(r)
    if not roots:
        raise RuntimeError("no root bracketed for f(pi)=0")
    # physical root: the one whose f'(pi) is smallest in magnitude
    p0 = min(roots, key=lambda r: abs(fend(s_guess, r)[1]))

    # 2-parameter Newton from (s_guess, p0)
    x = np.array([s_guess, p0])
    for _ in range(60):
        F = np.array(fend(x[0], x[1], rtol=1e-11, mstep=0.004))
        J = np.zeros((2, 2))
        for j in range(2):
            dx = np.zeros(2)
            dx[j] = 1e-7
            J[:, j] = (np.array(fend(*(x + dx), rtol=1e-11, mstep=0.004)) - F) / 1e-7
        try:
            step = np.linalg.solve(J, -F)
        except np.linalg.LinAlgError:
            break
        x = x + step
        if np.linalg.norm(step) < 1e-13:
            break
    return x[0], x[1]


S_FORMULA = (2 * n + 1) / (n + 1)
s_num, p_shoot = shoot(n)
print(f"[eigenvalue] s = {s_num:.12f}, formula (2n+1)/(n+1) = {S_FORMULA:.12f}, "
      f"diff = {abs(s_num-S_FORMULA):.2e}")
print(f"[eigenvalue] f''(0) = {p_shoot:.12f}")

# ---- angular field on a grid ----
s = S_FORMULA
NT = 40001
theta = np.linspace(0.0, np.pi, NT)
sol = integrate(n, s, p_shoot)
Y = sol.sol(theta)
f, fp, fpp, fppp = Y
M = Mcoeff(s)
Z = Y[:3].T
se = np.sqrt(np.einsum('ij,jk,ik->i', Z, M, Z))
se_max = se.max()
c = 1.0 / se_max          # normalise max se_tilde = 1
Y = Y * c
f, fp, fpp, fppp = Y
sr = s * f + fpp
st = s * (s - 1.0) * f
srt = (1.0 - s) * fp
se = np.sqrt(np.einsum('ij,jk,ik->i', Y[:3].T, M, Y[:3].T))

# =====================================================================
# J and I (eq. 24 with symmetry factor 2), plane strain
# =====================================================================
m = n * (s - 2.0)          # = -(n/(n+1))
bb = n * (s - 2.0) + 1.0   # = 1/(n+1)
rho = se ** (n - 1.0)
z = np.stack([f, fp, fpp], -1)
dQ = 2.0 * (z @ M.T)
dQv = np.einsum('ij,ij->i', dQ, np.stack([fp, fpp, fppp], -1))
rhop = (n - 1.0) * rho * (dQv / (2.0 * se ** 2))
c1 = 0.75 * (sr - st)
c1p = 0.75 * (s * (2.0 - s) * fp + fppp)
ur = rho * c1 / bb
urp = (rhop * c1 + rho * c1p) / bb
ut = (3.0 * rho * srt - urp) / m
utp = rho * (0.75 * (st - sr)) - ur

integrand = ((n / (n + 1.0)) * se ** (n + 1.0) * np.cos(theta)
             - (np.sin(theta) * (sr * (ut - urp) - srt * (ur + utp))
                + bb * np.cos(theta) * (sr * ur + srt * ut)))
I = 2.0 * np.trapezoid(integrand, theta)
print(f"[I] I = {I:.10f}")

# =====================================================================
# Dimensionless setup
# ---------------------------------------------------------------------
# Non-dimensionalise lengths by some r0, stresses by sig0, strains by eps0=alpha?? 
# The natural normalisation in the task:
#   constitutive: eps_ij = (3/2) alpha eps0 (se/sig0)^{n-1} s_ij/sig0
#   phi = K r^s f(theta)  =>  sigma = K r^{s-2} * (angular)
# Let us keep everything dimensional and just trace the algebra.
#
# W = n/(n+1) sigma_ij eps_ij
# sigma_ij eps_ij = (3/2) alpha eps0 (se/sig0)^{n-1} (s_ij s_ij)/sig0
#                 = (3/2) alpha eps0 (se/sig0)^{n-1} * (2/3) se^2 / sig0
#                 = alpha eps0 se^{n+1} / sig0^n
# => W = (n/(n+1)) * alpha eps0 * se^{n+1} / sig0^n
#
# se = K r^{s-2} se_tilde(theta), (n+1)(s-2) = -1
# se^{n+1} = K^{n+1} r^{-1} se_tilde^{n+1}
# W = (n/(n+1)) alpha eps0/sig0^n * K^{n+1} r^{-1} se_tilde^{n+1}
#
# U(R) = int_0^R int_{-pi}^{pi} W r dtheta dr
#      = (n/(n+1)) alpha eps0/sig0^n * K^{n+1} * [int_0^R dr] * shape(2x)
#      = (n/(n+1)) alpha eps0/sig0^n * K^{n+1} * R * shape2
#   shape2 = int_{-pi}^{pi} se_tilde^{n+1} dtheta = 2*int_0^{pi} se_tilde^{n+1} dtheta
#
# J = alpha eps0 sig0 * K^{n+1} * I   (I as computed, J/m^2; see below)
# Actually: J has units of stress*length. The inner-contour identity in the
# paper is J = alpha K^{n+1} r^{(n+1)(s-2)+1} I with K carrying
# sigma^(-1)*length^{2-s}... i.e. the paper is non dimensional. Dimensionally:
# se is a stress ~ K_nondim * (sig0?) ... Let us express K dimensionful, Kdim,
# so sigma = Kdim r^{s-2} angular. Then
#   W = (n/(n+1)) alpha eps0/sig0^n Kdim^{n+1} r^{-1} se_tilde^{n+1}
#   J = alpha eps0/sig0^n * Kdim^{n+1} * I     (the eps0 sig0 from the traction
#                                               term combine to the same units)
# Check: the J integral = int (W r dtheta) - (t d u/dx) r dtheta. Both terms
# scale as W*r ~ (stress)*length and W has units of stress.
# W ~ Kdim^{n+1} * (alpha eps0/sig0^n) * r^{-1}
# So J ~ Kdim^{n+1}*(alpha eps0/sig0^n)*I  [units: (stress^{n+1})*(eps0/sig0^n)]
#     = stress^{n+1} * (stress/(E*sig0^n)) = stress^{n+2}/(E sig0^n). Hmm units.
# This is getting unwieldy; do it cleanly with the substitution below.
# =====================================================================

# Cleaner: use the task's target form directly.
#   sigma_tt(r,0) = H sig0 [J/(alpha eps0 sig0 r)]^p
# The bracket J/(alpha eps0 sig0 r) must be dimensionless:
#   J [J/m^2 = N/m] / (eps0 [1] * sig0 [Pa=N/m^2] * r [m]) = (N/m)/(N/m) = 1. good.
#
# From the field:
#   sigma_tt(r,0) = Kdim * st0 * r^{s-2} = Kdim * st0 * r^{-p},  p=1/(n+1)
#   J = alpha *(eps0/sig0^n)* Kdim^{n+1} * I   [derive from J ~ W*r scaling]
# Solve for Kdim:
#   Kdim^{n+1} = J * sig0^n / (alpha eps0 I)
#   Kdim = [J sig0^n/(alpha eps0 I)]^{1/(n+1)}
# Then
#   sigma_tt(r,0) = st0 [J sig0^n/(alpha eps0 I)]^{p} r^{-p}
#   target        = H sig0 [1/(alpha eps0 sig0)]^p J^p r^{-p}
# Match:
#   st0 [sig0^n/(alpha eps0 I)]^p = H sig0 [1/(alpha eps0 sig0)]^p
#   st0 sig0^{np} (alpha eps0 I)^{-p} = H sig0 (alpha eps0 sig0)^{-p}
#   H = st0 sig0^{np-1} (alpha eps0 I)^{-p} (alpha eps0 sig0)^{p}
#     = st0 I^{-p} sig0^{np-1} (sig0)^{p}
#     = st0 I^{-p} sig0^{np-1+p} = st0 I^{-p} sig0^{p(n+1)-1}
#   with p(n+1)=1:  H = st0 I^{-p} sig0^{0} = st0 * I^{-p}
# Beautiful: H = st0 * I^{-p}, independent of sig0, eps0, alpha magnitudes!
# =====================================================================

p = 1.0 / (n + 1.0)
st0 = st[0]           # angular sigma_tt at theta=0, normalised max se_tilde=1
print(f"[field] st~(0) = {st0:.12f}, f(0) = {f[0]:.12f}, se_max_raw={se_max:.12f}")

H = st0 * I ** (-p)
print(f"[p] p = {p:.12f}")
print(f"[H] H = st0 * I^(-p) = {H:.12f}")

# --- direct numeric verification of H through full dimensional chain ---
def sigma_tt_dim(r, Kdim):
    return Kdim * st0 * r ** (s - 2.0)

def J_from_Kdim(Kdim):
    return alpha * eps0 / sig0 ** n * Kdim ** (n + 1.0) * I

def Kdim_from_J(J):
    return (J * sig0 ** n / (alpha * eps0 * I)) ** (1.0 / (n + 1.0))

Jin = 12345.678       # arbitrary J in J/m^2
Kdim = Kdim_from_J(Jin)
print(f"[check] J_from_Kdim(Kdim_from_J(J)) / J = {J_from_Kdim(Kdim)/Jin:.12f}")
r0 = 1e-4
sig_val = sigma_tt_dim(r0, Kdim)
target_factor = sig0 * (Jin / (alpha * eps0 * sig0 * r0)) ** p
H_check = sig_val / target_factor
print(f"[check] H via dimensional chain = {H_check:.12f}")

# =====================================================================
# C = U(R)/(J R)
#   U(R) = (n/(n+1)) (alpha eps0/sig0^n) Kdim^{n+1} R * shape2
#   shape2 = int_{-pi}^{pi} se_tilde^{n+1} dtheta
#   J     = alpha eps0/sig0^n Kdim^{n+1} I
#   C = (n/(n+1)) R shape2 / (I R) = (n/(n+1)) * shape2 / I
# =====================================================================
shape2 = 2.0 * np.trapezoid(se ** (n + 1.0), theta)
C = (n / (n + 1.0)) * shape2 / I
print(f"[shape2] int se_tilde^(n+1) dtheta over [-pi,pi] = {shape2:.12f}")
print(f"[C] C = (n/(n+1)) shape2 / I = {C:.12f}")

print("\n==== RESULTS ====")
print(f"p = {p:.10f}")
print(f"H = {H:.10f}")
print(f"C = {C:.10f}")