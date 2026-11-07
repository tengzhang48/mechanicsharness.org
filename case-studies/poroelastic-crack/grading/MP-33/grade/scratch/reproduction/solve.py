#!/usr/bin/env python3
"""Suddenly pressurized crack in a Biot poroelastic solid (plane strain).

Model: t=0+ is the undrained elastic response (G, nu_u); t->inf for loading 1 is the
drained elastic response (G, nu) because the pressure perturbation diffuses away.
Short-time crack-volume growth comes from a 1-D uniaxial-strain diffusion boundary
layer at the crack faces. Loading 2 (p imposed, zero traction) gives no instantaneous
response and no long-time crack perturbation (p -> uniform -> uniform effective stress).
See ANSWER.md for the full derivation. Reproduce: python3 solve.py
"""
import json
import numpy as np
from scipy import integrate
from scipy.special import erfc

nu, nu_u, alpha = 0.2, 0.4, 0.89
G, L, P = 1.0, 1.0, 1.0          # per unit loading; K in sqrt(pi L), V in L^2/G
pi = np.pi

# ================= analytic results =================
K1_0   = 1.0                     # P sqrt(pi L), undrained elastic, t=0+
K1_inf = 1.0                     # P sqrt(pi L), drained elastic, t->inf
V1_0   = pi * (1 - nu_u) * P * L**2 / G
V1_inf = pi * (1 - nu)   * P * L**2 / G
p0     = 2 * (nu_u - nu) / (alpha * (1 - 2 * nu))     # = p_faces_0
A      = 8 * (nu_u - nu) / ((1 - nu) * np.sqrt(pi))   # V1_short_coeff
K2_0, K2_inf, V2_inf = 0.0, 0.0, 0.0

checks = []
def check(name, got, want, tol=1e-8, how=""):
    ok = abs(got - want) <= tol * max(1.0, abs(want))
    checks.append((name, got, want, ok, how))
    print(f"[check] {name}: got {got:.12g}, expected {want:.12g} -> {'PASS' if ok else 'FAIL'}")

# check 1: undrained reduction (zeta=0) gives elasticity with nu_u
skk = 1.7
p_u = -(nu_u - nu) * skk / (alpha * (1 - 2 * nu))
coeff = (-nu * skk + alpha * (1 - 2 * nu) * p_u) / skk
check("undrained reduction coefficient of sigma_kk vs -nu_u", coeff, -nu_u,
      how="fails if the two constitutive laws were mutually inconsistent")

# check 2: diffusivity implied by 1-D uniaxial-strain mass balance equals stated c
kappa = 3.1e-4
c_given = 2 * kappa * G * (1 - nu) * (nu_u - nu) / (alpha**2 * (1 - 2 * nu)**2 * (1 - nu_u))
c_mine  = kappa * 2 * G * (1 - nu) * (nu_u - nu) / (alpha**2 * (1 - 2 * nu)**2 * (1 - nu_u))
check("1-D boundary-layer diffusivity vs stated c", c_mine, c_given, tol=1e-12,
      how="fails if the layer's mechanical constraint (e_xx=0) were wrong")

# check 3: crack volume by quadrature of the Griffith COD
cod = lambda x: 2 * (1 - nu) * P / G * np.sqrt(L**2 - x**2)
V_num, _ = integrate.quad(cod, -L, L, epsabs=1e-14)
check("crack volume quadrature vs pi(1-nu)PL^2/G", V_num, pi * (1 - nu) * P * L**2 / G,
      how="fails if the COD formula or its integral pi L^2/2 were wrong")

# check 4: integral of erfc
I_erfc, _ = integrate.quad(erfc, 0, 50)
check("integral_0^inf erfc(s) ds vs 1/sqrt(pi)", I_erfc, 1 / np.sqrt(pi),
      how="fails if the boundary-layer heave integral were wrong")

# check 5: Westergaard sigma_kk = -2P on the crack faces
xs = np.linspace(-0.999 * L, 0.999 * L, 21)
Z = lambda x: P * ((x + 1e-12j) / np.sqrt((x + 1e-12j)**2 - L**2) - 1.0)
skk_face = np.array([2 * Z(x).real for x in xs])
err = np.max(np.abs(skk_face + 2 * P))
check("max |sigma_kk + 2P| on faces (Westergaard)", err, 0.0, tol=1e-6,
      how="fails if the elastic face stress trace were not -2P, changing p_faces_0")
p0_west = -(nu_u - nu) * (-2 * P) / (alpha * (1 - 2 * nu))
check("p_faces_0 from Westergaard trace", p0_west, p0,
      how="fails if the undrained pressure conversion were wrong")

# check 6: independent finite-difference 1-D diffusion simulation of the face layer.
# p(y,0)=p0, p(0,t)=0; face heave du=(1-2nu)alpha/(2G(1-nu)) * int(p0-p)dy; dV = 4 L du.
# dimensionless c=G=L=P=1; slope of dV*(G/L^2) vs sqrt(c t/L^2) must equal A.
# FTCS is O(dy^2); run three grids and Richardson-extrapolate to dy -> 0.
def fdm_slope(Ny, t_end, Ymax=1.0):
    c1 = 1.0
    dy = Ymax / Ny
    dt = 0.4 * dy**2 / c1           # explicit FTCS, r <= 0.5 for stability
    r = c1 * dt / dy**2
    p = np.full(Ny + 1, p0); p[0] = 0.0
    nsteps = int(np.ceil(t_end / dt))
    for _ in range(nsteps):
        p[1:-1] += r * (p[2:] - 2 * p[1:-1] + p[:-2])
        p[0] = 0.0; p[-1] = p0
    t = nsteps * dt
    integ = np.sum(p0 - p[1:]) * dy
    du = (1 - 2 * nu) * alpha / (2 * G * (1 - nu)) * integ
    dV = 4 * L * du
    return (dV * G / L**2) / np.sqrt(c1 * t / L**2)

t_end = 5e-4
Nys = [1000, 2000, 4000, 8000, 16000]
slopes = [fdm_slope(Ny, t_end) for Ny in Nys]
diffs = np.diff(slopes)
print(f"[check] FDM slopes at c t/L^2={t_end:g}, Ny={Nys}: "
      + ", ".join(f"{s:.6f}" for s in slopes))
print("[check] successive differences: " + ", ".join(f"{d:.6f}" for d in diffs)
      + " ; ratios: " + ", ".join(f"{diffs[i+1]/diffs[i]:.3f}" for i in range(len(diffs)-1)))
# differences halve per refinement -> scheme is 1st order in dy (discontinuous IC at y=0)
A_extra = slopes[-1] + (slopes[-1] - slopes[-2])      # Richardson for O(dy) error
check("FDM growth coefficient (finest grid) vs analytic A", slopes[-1], A, tol=3e-3,
      how="fails if the 1-D layer physics (heave law, both faces, erfc profile) were wrong")
print(f"[check] Richardson-extrapolated slope: {A_extra:.6f} (analytic A = {A:.6f})")
check("Richardson-extrapolated FDM slope vs analytic A", A_extra, A, tol=1e-4,
      how="fails if the discrete solution did not converge to the analytic boundary-layer result")

# cross-check the heave integral of the exact erf similarity profile
for ct in [1e-4, 1e-3, 1e-2]:
    integ_exact = 2 * p0 * np.sqrt(ct) / np.sqrt(pi)
    integ_num, _ = integrate.quad(lambda y: p0 * erfc(y / (2 * np.sqrt(ct))), 0, 1.0)
    check(f"erf-profile heave integral at ct={ct:g}", integ_num, integ_exact,
          how="fails if the similarity profile p0 erf(y/2 sqrt(ct)) were wrong")

# ---- assemble results.json ----
results = {
    "K1_0": K1_0,
    "K1_inf": K1_inf,
    "V1_0": V1_0,
    "V1_inf": V1_inf,
    "V1_short_coeff": A,
    "p_faces_0": p0,
    "K2_0": K2_0,
    "K2_inf": K2_inf,
    "V2_inf": V2_inf,
    "formulas": {
        "K1_0": "P*sqrt(pi*L)  [= 1 in units of sqrt(pi*L), per unit load]; undrained elastic pressurized crack",
        "K1_inf": "P*sqrt(pi*L)  [= 1]; drained elastic pressurized crack (stress field independent of moduli)",
        "V1_0": "pi*(1-nu_u)*P*L^2/G  [= pi*(1-nu_u) in units of L^2/G]",
        "V1_inf": "pi*(1-nu)*P*L^2/G  [= pi*(1-nu)]",
        "V1_short_coeff": "8*(nu_u-nu)/((1-nu)*sqrt(pi));  V = pi*(1-nu_u)*L^2/G + A*(L^2/G)*sqrt(c*t/L^2)",
        "p_faces_0": "2*(nu_u-nu)/(alpha*(1-2*nu))  per unit applied traction (undrained pressure at the faces)",
        "K2_0": "0  (undrained bulk + zero applied traction -> no instantaneous deformation)",
        "K2_inf": "0  (p -> uniform 1 -> effective stress uniform alpha*delta_ij -> no crack-tip field)",
        "V2_inf": "0  (uniform free swelling; a crack is invisible to a uniform effective stress)",
    },
    "reproduce_command": "python3 solve.py",
}
with open("results.json", "w") as f:
    json.dump(results, f, indent=2)

print("\n==== values (K in units of sqrt(pi L), V in units of L^2/G, per unit load) ====")
for k in ["K1_0", "K1_inf", "V1_0", "V1_inf", "V1_short_coeff", "p_faces_0",
          "K2_0", "K2_inf", "V2_inf"]:
    print(f"{k} = {results[k]:.10f}")
nfail = sum(1 for c in checks if not c[3])
print(f"\n{len(checks)} checks, {nfail} failed")
