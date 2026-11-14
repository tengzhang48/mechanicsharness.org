"""
A crack front trapped by tough bands under mixed-mode loading.

First-order (in front-slope) perturbation theory for a planar crack in an
infinite isotropic elastic body.  All quantities are computed analytically;
Fourier-series checks are included at the end.

Part A: prescribed front a(z) = abar + A cos(2*pi*z/lam), A/lam << 1.
Part B: equilibrium front a(z) = abar + eps*lam*f(z/lam) in a periodic
        toughness Gamma(z) = Gamma0 [1 + eps*g(z/lam)].

Usage:  python3 solve.py
Writes: results.json
"""

import json
import numpy as np

# ----------------------------------------------------------------------
# Material and loading
# ----------------------------------------------------------------------
nu = 0.3
K2r = 0.6          # K2^0 / K1^0
K3r = 0.8          # K3^0 / K1^0
q = 0.5            # toughness-band parameter

# Irwin relation: G = (1-nu^2)(K1^2+K2^2)/E + (1+nu)K3^2/E
Gcoef = (1.0 - nu**2) * (1.0 + K2r**2) + (1.0 + nu) * K3r**2

# weights w_a = (elastic factor * (K_a^0)^2) / Gcoef
w1 = (1.0 - nu**2) * 1.0 / Gcoef
w2 = (1.0 - nu**2) * K2r**2 / Gcoef
w3 = (1.0 + nu) * K3r**2 / Gcoef

# Geometric (mode-rotation) coefficient in the energy release rate:
#   dG/G0 = 2*F + beta*f'   (F is the common relative nonlocal K perturbation)
beta = 2.0 * (-(4.0 / 3.0) * w2 + 0.75 * w3)

# ----------------------------------------------------------------------
# Part A.  Prescribed cosine front
#
# Nonlocal Rice kernel for a Fourier mode:
#   dK_a/K_a^0 = -(1/2) |k| da          (|k| = 2*pi/lam)
# With da = A cos(kz) this gives the relative perturbation
#   F = -(pi A/lam) cos(kz).
# The rotation of the local (n,t) frame adds
#   dK2 = -a' K3^0,  dK3 = +a' K2^0,  a' = -(2*pi A/lam) sin(kz).
# ----------------------------------------------------------------------
C1 = -np.pi
S1 = 0.0
C2 = K2r * (-np.pi)
S2 = 2.0 * np.pi * K3r          # = 0.6 * (8*pi/3)
C3 = K3r * (-np.pi)
S3 = -2.0 * np.pi * K2r         # = 0.8 * (-3*pi/2)
CG = -2.0 * np.pi
SG = -2.0 * np.pi * beta

partA = {
    "C1": float(C1), "S1": float(S1),
    "C2": float(C2), "S2": float(S2),
    "C3": float(C3), "S3": float(S3),
    "CG": float(CG), "SG": float(SG),
}

# ----------------------------------------------------------------------
# Part B.  Equilibrium front in periodic tough bands
#
# Toughness:  g(s) = 2*sum_{n>=1} q^n cos(2*pi*n*s)
# Equilibrium: G = Gamma  =>  -H[f'] + beta f' = g
# Fourier solution:
#   f_n = q^n (-1 - i beta) / (2*pi*n*(1+beta^2))
#   f(s) = [ (1/2) ln D(s) - beta*phi(s) ] / (pi*(1+beta^2))
# with D = 1 - 2q cos(2*pi*s) + q^2, phi = arg(1 - q exp(2*pi*i*s)).
#
# The reported kappa_a are defined by K_a = K_a^0 + eps*K1^0*kappa_a:
#   kappa1 = -H[f']/2
#   kappa2 = 0.6*kappa1 - 0.8*f'
#   kappa3 = 0.8*kappa1 + 0.6*f'
# ----------------------------------------------------------------------
def theta(s):
    return 2.0 * np.pi * s


def Dfun(s):
    t = theta(s)
    return 1.0 - 2.0 * q * np.cos(t) + q * q


def phifun(s):
    t = theta(s)
    return np.arctan2(-q * np.sin(t), 1.0 - q * np.cos(t))


def gfun(s):
    t = theta(s)
    return 2.0 * q * (np.cos(t) - q) / Dfun(s)


def ffun(s):
    return (0.5 * np.log(Dfun(s)) - beta * phifun(s)) / (np.pi * (1.0 + beta**2))


def fpfun(s):
    t = theta(s)
    return ((2.0 * q / (1.0 + beta**2)) *
            (np.sin(t) - beta * (q - np.cos(t))) / Dfun(s))


def fppfun(s):
    t = theta(s)
    dd = Dfun(s)
    N = np.sin(t) - beta * (q - np.cos(t))
    Np = np.cos(t) - beta * np.sin(t)
    return (2.0 * q / (1.0 + beta**2)) * 2.0 * np.pi * \
        (Np * dd - N * 2.0 * q * np.sin(t)) / dd**2


def gpfun(s):
    t = theta(s)
    dd = Dfun(s)
    return -4.0 * np.pi * q * (1.0 - q * q) * np.sin(t) / dd**2


def kappas(s):
    fp = fpfun(s)
    k1 = 0.5 * (gfun(s) - beta * fp)
    k2 = K2r * k1 - K3r * fp
    k3 = K3r * k1 + K2r * fp
    return k1, k2, k3


# high-precision extrema of f and kappa3
from scipy.optimize import brentq, minimize_scalar

sgrid = np.linspace(-0.5, 0.5, 20001)
froots = []
for i in range(len(sgrid) - 1):
    a, b = sgrid[i], sgrid[i + 1]
    if fpfun(a) * fpfun(b) < 0.0:
        froots.append(brentq(fpfun, a, b, xtol=1e-15, rtol=1e-15))
fvals = [ffun(r) for r in froots]
f_min, f_max = min(fvals), max(fvals)
front_peak_to_peak = f_max - f_min

# kappa3'(s) = 0.4*g'(s) + (0.6 - 0.4*beta)*f''(s)
def k3p(s):
    return 0.4 * gpfun(s) + (0.6 - 0.4 * beta) * fppfun(s)


k3roots = []
for i in range(len(sgrid) - 1):
    a, b = sgrid[i], sgrid[i + 1]
    if k3p(a) * k3p(b) < 0.0:
        k3roots.append(brentq(k3p, a, b, xtol=1e-15, rtol=1e-15))
k3cand = [(kappas(r)[2], r) for r in k3roots]
kappa3_max, s_kappa3_max = max(k3cand, key=lambda x: x[0])
# include endpoints of the scanned interval in case the maximum sits there
for r in (-0.5, 0.5):
    val = kappas(r)[2]
    if val > kappa3_max:
        kappa3_max, s_kappa3_max = val, r

ss = np.arange(16, dtype=float) / 16.0
f16 = [float(ffun(s)) for s in ss]
kappa1_16 = [float(kappas(s)[0]) for s in ss]
kappa2_16 = [float(kappas(s)[1]) for s in ss]
kappa3_16 = [float(kappas(s)[2]) for s in ss]

K10_scaled = float(np.sqrt((1.0 - nu**2) / Gcoef))

partB = {
    "K10_scaled": K10_scaled,
    "front_peak_to_peak": float(front_peak_to_peak),
    "f": f16,
    "kappa1": kappa1_16,
    "kappa2": kappa2_16,
    "kappa3": kappa3_16,
    "kappa3_max": float(kappa3_max),
    "s_kappa3_max": float(s_kappa3_max),
}

results = {"partA": partA, "partB": partB}

with open("results.json", "w") as fh:
    json.dump(results, fh, indent=2)
    fh.write("\n")

# ----------------------------------------------------------------------
# Checks
# ----------------------------------------------------------------------
def check_fourier(M=2048):
    """Check -H[f'] + beta f' = g to machine precision."""
    s = np.arange(M) / M
    fv = np.array([ffun(x) for x in s])
    gv = np.array([gfun(x) for x in s])
    F = np.fft.rfft(fv)                 # unnormalised: M * fhat_n
    n = np.arange(len(F))
    Hf = np.fft.irfft(2.0 * np.pi * n * F, n=M)
    fp = np.fft.irfft((2j * np.pi * n) * F, n=M)
    res = -Hf + beta * fp - gv
    return float(np.max(np.abs(res)))


def check_equilibrium_points():
    """First-order G=Gamma residual at the 16 reported points for small eps."""
    eps = 1e-3
    worst = 0.0
    for s in ss:
        k1, k2, k3 = kappas(s)
        # relative perturbations
        d1 = eps * k1
        d2 = eps * k2 / K2r
        d3 = eps * k3 / K3r
        dG = 2.0 * (w1 * d1 + w2 * d2 + w3 * d3)
        worst = max(worst, abs(dG - eps * gfun(s)))
    return worst


if __name__ == "__main__":
    print("Part A:", partA)
    print("Part B: K10_scaled =", K10_scaled)
    print("        front_peak_to_peak =", front_peak_to_peak)
    print("        kappa3_max =", kappa3_max, "at s =", s_kappa3_max)
    print("Fourier check (max |residual|):", check_fourier())
    print("Equilibrium check (max residual, eps=1e-3):", check_equilibrium_points())
