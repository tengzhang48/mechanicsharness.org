"""Separation between unequal strips -- beam (DCB) cohesive-zone model.

Geometry (test specimen), from the free ends on the left moving right:
    free ends (left) -- starter cut (a0 = 20 mm, no connectors) --
    root (x = L) -- bridging zone (length L, connectors) -- front (x = 0) --
    bonded region (long, straight).

Coordinate x runs from the front (x = 0) toward the root (x = L) and the
starter cut (L < x < L + a0); the free ends are at x = L + a0.

The two strips are Euler-Bernoulli beams clamped at the front (the long,
straight bonded region does not rotate).  Opening delta(x) = w1 + w2.

    D~ delta'''' = -p(x),   1/D~ = 1/D1 + 1/D2,
    p(x) = b*t(delta) on 0 < x < L  (connectors), 0 for x > L (starter cut),
    t(delta) = sb*(1 - delta/dc) on 0 < delta < dc, else 0 (inactive).

BC:  delta(0) = delta'(0) = 0          (clamped at front)
     D~ delta''(L+a0) = M, delta'''(L+a0) = 0   (pure opening couple M per strip)

Fracture (J-integral, path independence):
     G_far = M^2/(2 b D~)                      (pure-moment DCB, length independent)
     G_tip = G_far - int_0^{delta_t} t dd    =  Gamma_0   at front advance
State of interest: trailing-edge separation delta_t = delta(L) = dc.
"""
import numpy as np

# ---- material / geometry ---------------------------------------------------
E1, h1 = 110e9, 0.60e-3      # upper strip
E2, h2 = 190e9, 0.90e-3      # lower strip
b = 25.0e-3                  # width (in-plane)
a0 = 20.0e-3                 # starter cut
sb = 0.20e6                  # bridging strength
dc = 0.20e-3                 # limiting separation
Gam0 = 1.40                  # intrinsic toughness J/m^2

I1 = b * h1**3 / 12.0
I2 = b * h2**3 / 12.0
D1 = E1 * I1
D2 = E2 * I2
Dtil = D1 * D2 / (D1 + D2)

# ---- length scale ----------------------------------------------------------
Lc = (Dtil * dc / (b * sb))**0.25        # bending-cohesive length
k = 1.0 / Lc                             # wave number (linear softening law)

print("D1 =", D1, " D2 =", D2, " D~ =", Dtil)
print("Lc (mm) =", Lc * 1e3)

# ---- cohesive energy for the softening law --------------------------------
def bridging_work(dt):
    """int_0^{dt} t(delta) d delta  (J/m^2)."""
    if dt <= 0:
        return 0.0
    dt = min(dt, dc)
    return sb * dt - sb * dt**2 / (2 * dc)

# ---- analytic linear BVP ---------------------------------------------------
def solve_delta(L, M):
    """Solve the linear BVP; return coefficients and delta(x) evaluator.

    Bridging zone: delta = dc + A e^{kx} + B e^{-kx} + C cos kx + D sin kx
    Starter cut : delta = P x^3 + Q x^2 + R x + S
    """
    x1 = L
    x2 = L + a0

    def basis_zone(x):
        return np.array([np.exp(k * x), np.exp(-k * x),
                         np.cos(k * x), np.sin(k * x)])

    def dbasis_zone(x):
        return np.array([k * np.exp(k * x), -k * np.exp(-k * x),
                         -k * np.sin(k * x), k * np.cos(k * x)])

    def ddbasis_zone(x):
        return np.array([k**2 * np.exp(k * x), k**2 * np.exp(-k * x),
                         -k**2 * np.cos(k * x), -k**2 * np.sin(k * x)])

    def dddbasis_zone(x):
        return np.array([k**3 * np.exp(k * x), -k**3 * np.exp(-k * x),
                         k**3 * np.sin(k * x), -k**3 * np.cos(k * x)])

    def basis3(x):
        return np.array([x**3, x**2, x, 1.0])

    def dbasis3(x):
        return np.array([3 * x**2, 2 * x, 1.0, 0.0])

    def ddbasis3(x):
        return np.array([6 * x, 2.0, 0.0, 0.0])

    def dddbasis3(x):
        return np.array([6.0, 0.0, 0.0, 0.0])

    # unknowns u = [A,B,C,D,P,Q,R,S]
    mat = np.zeros((8, 8))
    rhs = np.zeros(8)
    # x = 0: delta(0)=0 -> dc + A + B + C = 0
    mat[0, :4] = [1, 1, 1, 0]
    rhs[0] = -dc
    # x = 0: delta'(0)=0 -> A - B + D = 0
    mat[1, :4] = [k, -k, 0, k]
    rhs[1] = 0.0
    # x = L continuity delta
    mat[2, :4] = basis_zone(x1)
    mat[2, 4:] = -basis3(x1)
    rhs[2] = -dc
    # delta'
    mat[3, :4] = dbasis_zone(x1)
    mat[3, 4:] = -dbasis3(x1)
    # delta''
    mat[4, :4] = ddbasis_zone(x1)
    mat[4, 4:] = -ddbasis3(x1)
    # delta'''
    mat[5, :4] = dddbasis_zone(x1)
    mat[5, 4:] = -dddbasis3(x1)
    # free end x2: delta'' = M/D~
    mat[6, 4:] = ddbasis3(x2)
    rhs[6] = M / Dtil
    # free end x2: delta''' = 0
    mat[7, 4:] = dddbasis3(x2)
    rhs[7] = 0.0

    u = np.linalg.solve(mat, rhs)
    A, B, C, D, P, Q, R, S = u

    def delta(x):
        x = np.atleast_1d(np.asarray(x, dtype=float))
        out = np.empty_like(x)
        for i, xi in enumerate(x):
            if xi <= L:
                out[i] = dc + A * np.exp(k * xi) + B * np.exp(-k * xi) \
                    + C * np.cos(k * xi) + D * np.sin(k * xi)
            else:
                out[i] = P * xi**3 + Q * xi**2 + R * xi + S
        return out

    return delta, u

# ---- M from the tip condition at delta_t = dc ------------------------------
def M_from_tip(dt=dc):
    G = Gam0 + bridging_work(dt)
    return np.sqrt(2.0 * b * Dtil * G), G

# ---- find L such that delta(L) = dc ----------------------------------------
def L_for_M(M, dt_target=dc):
    from scipy.optimize import brentq

    def f(L):
        delta, _ = solve_delta(L, M)
        return delta(L + 0.0)[0] - dt_target   # delta at the root x=L

    # bracket
    lo, hi = 1e-6, 1e-2
    while f(hi) < 0:
        hi *= 1.5
        if hi > 1.0:
            raise RuntimeError("no bracket")
    return brentq(f, lo, hi, xtol=1e-14)

if __name__ == "__main__":
    M, G = M_from_tip()
    print("G_far required (J/m^2) =", G)
    print("M (N m) =", M)
    L = L_for_M(M)
    print("L (mm) =", L * 1e3)
    delta, u = solve_delta(L, M)
    print("delta(L) (um) =", delta(L)[0] * 1e6)
    print("delta(L/2) (um) =", delta(L / 2.0)[0] * 1e6)
    print("delta(free end) (um) =", delta(L + a0)[0] * 1e6)