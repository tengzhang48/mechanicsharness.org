#!/usr/bin/env python3
"""
MP-22: A suddenly pressurized crack in a fluid-saturated porous solid.

Computes the quantities requested in TASK.md part (c) and writes results.json.
Units: K in sqrt(pi*L), V in L^2/G, per unit loading magnitude (sigma0 = 1
or p_f = 1).  G, L, c are set to 1 in the numerical checks.

Physics used (see ANSWER.md for derivations)
---------------------------------------------
* Biot poroelasticity. Instantaneous response (t = 0+) = undrained elastic
  (G, nu_u); long-time response (t -> inf) = drained elastic (G, nu).
* Griffith crack (half-length L) under unit opening traction:
    K_I = sqrt(pi L)                    (independent of elastic moduli)
    V   = pi (1 - nu') L^2 / G          (nu' = operative Poisson ratio)
* Loading 1 (sigma_n = -H(t), p = 0 on faces):
    t=0+  : undrained crack. Undrained overpressure at faces
            p_b = (nu_u - nu) sigma0 / (alpha (1 - 2 nu))  [from zeta = 0,
            sigma_kk = -sigma0 at faces].
    t=inf : p diffused away -> drained crack.
    short : draining of p_b through a 1-D oedometric diffusion boundary
            layer of thickness sqrt(c t) contracts the near-face skeleton
            and pulls the faces apart:
            V ~ V1_0 + A sqrt(c t/L^2), A = 4 (nu_u - nu)/((1-nu) sqrt(pi)).
* Loading 2 (sigma_n = 0, p = H(t) on faces):
    t=0+  : p = 0 in the bulk; the effective face traction alpha*p_f is
            cancelled by the boundary body force in the zero-thickness
            layer -> no crack loading: K2_0 = 0.
    t=inf : 2-D exterior diffusion with maintained face value -> p -> p_f
            everywhere; uniform effective stress, uniform strain, faces do
            not separate: K2_inf = V2_inf = 0.
"""
import json
import math

import numpy as np

# ---------------- given parameters (part c) ----------------
nu = 0.2      # drained Poisson ratio
nuu = 0.4     # undrained Poisson ratio
alp = 0.89    # Biot coefficient
pi = math.pi

# ---------------- analytic results (normalized) ----------------
K1_0 = 1.0                     # sqrt(pi*L) per unit traction
K1_inf = 1.0
V1_0 = pi * (1.0 - nuu)        # L^2/G per unit traction
V1_inf = pi * (1.0 - nu)
A1 = 4.0 * (nuu - nu) / ((1.0 - nu) * math.sqrt(pi))
pf0 = (nuu - nu) / (alp * (1.0 - 2.0 * nu))   # per unit traction
K2_0 = 0.0
K2_inf = 0.0
V2_inf = 0.0

results = {
    "K1_0": K1_0,
    "K1_inf": K1_inf,
    "V1_0": V1_0,
    "V1_inf": V1_inf,
    "V1_short_coeff": A1,
    "p_faces_0": pf0,
    "K2_0": K2_0,
    "K2_inf": K2_inf,
    "V2_inf": V2_inf,
    "formulas": {
        "K1_0": "K = sigma0*sqrt(pi*L) (undrained elastic crack; K set by traction BC, independent of moduli); per unit loading in sqrt(pi*L): 1",
        "K1_inf": "K = sigma0*sqrt(pi*L) (drained elastic crack, p -> 0 everywhere); per unit loading in sqrt(pi*L): 1",
        "V1_0": "V = pi*(1-nu_u)*sigma0*L^2/G; in L^2/G per unit loading: pi*(1-nu_u)",
        "V1_inf": "V = pi*(1-nu)*sigma0*L^2/G; in L^2/G per unit loading: pi*(1-nu)",
        "V1_short_coeff": "V ~= V1_0 + A*sqrt(c*t/L^2)*sigma0*L^2/G with A = 4*(nu_u-nu)/((1-nu)*sqrt(pi)) (1-D oedometric drainage layer at the faces)",
        "p_faces_0": "p = (nu_u-nu)*sigma0/(alpha*(1-2*nu)) (from zeta=0 with sigma_kk=-sigma0 at faces); per unit loading: (nu_u-nu)/(alpha*(1-2*nu))",
        "K2_0": "0 (no bulk pressure yet; effective face traction alpha*p_f cancelled by boundary-layer body force)",
        "K2_inf": "0 (p -> p_f everywhere in 2-D at long time; uniform effective stress gives no crack-tip singularity)",
        "V2_inf": "0 (uniform strain state u_y = e_yy*y leaves the crack faces unseparated)",
    },
    "reproduce_command": "python3 solve.py",
}

print("=== results.json values (K in sqrt(pi L), V in L^2/G, per unit loading) ===")
for key in ("K1_0", "K1_inf", "V1_0", "V1_inf", "V1_short_coeff",
            "p_faces_0", "K2_0", "K2_inf", "V2_inf"):
    print(f"{key:15s} = {results[key]:.12g}")

with open("results.json", "w") as fh:
    json.dump(results, fh, indent=2)
print("results.json written")

# ======================================================================
# CHECKS
# ======================================================================
print("\n=== checks ===")
trap = getattr(np, "trapezoid", None) or np.trapz
from scipy.special import erf, erfc
fails = 0

def report(name, got, want, tol):
    global fails
    err = abs(got - want) / (abs(want) if want else 1.0)
    ok = err <= tol
    if not ok:
        fails += 1
    print(f"[{'PASS' if ok else 'FAIL'}] {name}: got {got:.12g}, want {want:.12g}, rel err {err:.2e}")

# Check 1: crack volume by numerical integration of the COD ellipse
#   V = int_{-L}^{L} 2(1-nu') sigma0/G sqrt(L^2-x^2) dx  vs  pi(1-nu')L^2/G
xs = np.linspace(-1.0, 1.0, 2_000_001)
for nup, tag in ((nu, "drained"), (nuu, "undrained")):
    cod = 2.0 * (1.0 - nup) * np.sqrt(np.clip(1.0 - xs**2, 0.0, None))
    report(f"COD quadrature for V ({tag})", trap(cod, xs), pi * (1.0 - nup), 1e-9)

# Check 2: integral of erfc used in the boundary-layer analysis
#   int_0^inf erfc(s) ds = 1/sqrt(pi)
s = np.linspace(0.0, 20.0, 4_000_001)
report("int_0^inf erfc(s) ds", trap(erfc(s), s), 1.0 / math.sqrt(pi), 1e-9)

# Check 3: oedometric (uniaxial-strain) response of the 2-D constitutive law.
# Impose e_xx = 0, sigma_yy fixed; verify de_yy/dp = alpha(1-2nu)/(2G(1-nu)).
G = 1.0
syy0 = -1.0
def eyy_oedo(p):
    # e_xx = [sxx - nu*(sxx+syy0) + alpha(1-2nu)p]/2G = 0  ->  sxx
    sxx = (nu * syy0 - alp * (1.0 - 2.0 * nu) * p) / (1.0 - nu)
    return (syy0 - nu * (sxx + syy0) + alp * (1.0 - 2.0 * nu) * p) / (2.0 * G)
dp = 1e-6
slope = (eyy_oedo(0.37 + dp) - eyy_oedo(0.37)) / dp
report("oedometric de_yy/dp", slope,
       alp * (1.0 - 2.0 * nu) / (2.0 * G * (1.0 - nu)), 1e-6)

# Check 4: given diffusivity c vs 1-D mass-balance derivation.
# Along the oedometric path with sigma_yy = 0: sigma_kk = -alpha(1-2nu)p/(1-nu),
# 2G zeta = alpha(1-2nu) sigma_kk + alpha^2(1-2nu)^2 p/(nu_u-nu),
# c = kappa / (dzeta/dp).  Compare with the formula quoted in TASK.md.
kap = 1.0
def zeta_oedo(p):
    skk = -alp * (1.0 - 2.0 * nu) * p / (1.0 - nu)
    return (alp * (1.0 - 2.0 * nu) * skk
            + alp**2 * (1.0 - 2.0 * nu)**2 * p / (nuu - nu)) / (2.0 * G)
dzdp_num = (zeta_oedo(0.53 + dp) - zeta_oedo(0.53)) / dp
c_derived = kap / dzdp_num
c_given = (2.0 * kap * G * (1.0 - nu) * (nuu - nu)
           / (alp**2 * (1.0 - 2.0 * nu)**2 * (1.0 - nuu)))
report("diffusivity c (derived vs given)", c_derived, c_given, 1e-6)

# Check 5: finite-difference 1-D diffusion simulation of the face boundary
# layer (loading 1, short time).  p(y,0)=p_b, p(0,t)=0, p(Y,t)=p_b.
#   dV = -2L alpha(1-2nu)/(G(1-nu)) * int (p - p_b) dy
#   A_num(t) = (dV / (L^2/G)) / sqrt(c t / L^2)  should equal A1.
L, c = 1.0, 1.0
pb = pf0
Y, dy = 2.0, 1.0e-3
y = np.arange(0.0, Y + 0.5 * dy, dy)
p = np.full(y.shape, pb)
p[0] = 0.0
dt = 0.4 * dy * dy / c
r = c * dt / dy**2
targets = [1e-3, 3e-3, 1e-2]
rec = {}
t = 0.0
nsteps = int(targets[-1] / dt) + 1
ti = 0
for n in range(1, nsteps + 1):
    p[1:-1] += r * (p[2:] - 2.0 * p[1:-1] + p[:-2])
    p[-1] = pb
    t = n * dt
    while ti < len(targets) and t >= targets[ti]:
        rec[targets[ti]] = p.copy()
        ti += 1
for tt in targets:
    pp = rec[tt]
    dV = (-2.0 * L * alp * (1.0 - 2.0 * nu) / (G * (1.0 - nu))
          * trap(pp - pb, y))
    A_num = (dV / (L**2 / G)) / math.sqrt(c * tt / L**2)
    report(f"FD diffusion short-time A at ct/L^2={tt:g}", A_num, A1, 2e-3)

# Check 6: same short-time amplitude with the analytic erf profile
# (quadrature check of int (p - p_b) dy = -p_b * 2 sqrt(c t / pi)).
ct = 3e-3
yy = np.linspace(0.0, 1.0, 2_000_001)
pprof = pb * erf(yy / (2.0 * math.sqrt(ct)))
dV_erf = (-2.0 * L * alp * (1.0 - 2.0 * nu) / (G * (1.0 - nu))
          * trap(pprof - pb, yy))
A_erf = (dV_erf / (L**2 / G)) / math.sqrt(ct / L**2)
report("analytic-profile short-time A", A_erf, A1, 1e-9)

print(f"\nchecks finished: {fails} failure(s)")
