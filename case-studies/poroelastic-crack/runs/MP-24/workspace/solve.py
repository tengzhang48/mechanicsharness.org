#!/usr/bin/env python3
"""
MP-24: Suddenly pressurized crack in a fluid-saturated poroelastic solid.

Model: Biot linear poroelasticity, plane strain (x-y), crack |x|<=L, y=0.
Loading 1: sigma_n = -H(t) (unit internal pressure), p = 0 on faces (drained).
Loading 2: sigma_n = 0, p = H(t) on faces.

Analytic results (derived in ANSWER.md), evaluated for part (c):
    nu = 0.2, nu_u = 0.4, alpha = 0.89
plus two independent numerical checks:
  CHECK 1: reciprocity consistency  int(sigma_kk^2 dA) = 2 pi L^2
  CHECK 2: 2-D slit-diffusion simulation reproducing the short-time
           coefficient A = 16 (nu_u-nu)/sqrt(pi)  [units L^2/G].
"""
import json
import numpy as np

# ---------------- part (c) parameters ----------------
nu, nu_u, alpha = 0.2, 0.4, 0.89

# ---------------- analytic report values ----------------
# K in units of sqrt(pi L); V in units of L^2/G; per unit loading magnitude.
K1_0   = 1.0                          # sqrt(pi L), undrained elastic SIF
K1_inf = 1.0                          # sqrt(pi L), drained elastic SIF
V1_0   = np.pi * (1.0 - nu_u)         # crack volume, undrained (t = 0+)
V1_inf = np.pi * (1.0 - nu)           # crack volume, drained (t -> inf)
A      = 16.0 * (nu_u - nu) / np.sqrt(np.pi)   # V ~ V1_0 + A sqrt(c t / L^2)
p_faces_0 = 2.0 * (nu_u - nu) / (alpha * (1.0 - 2.0 * nu))
K2_0, K2_inf, V2_inf = 0.0, 0.0, 0.0  # loading 2: see ANSWER.md

def westergaard_Z(X, Y):
    """Z = z/sqrt(z^2-1) - 1 for unit internal pressure, crack |x|<=1 (L=1).
    sigma_kk = 2 Re Z (independent of Poisson ratio)."""
    Z = X + 1j * Y
    with np.errstate(divide="ignore", invalid="ignore"):
        # sqrt(z^2-1) on the branch ~ z at infinity (cut on [-1,1] only):
        W = Z / (np.sqrt(Z - 1.0) * np.sqrt(Z + 1.0)) - 1.0
    return np.nan_to_num(W, nan=0.0, posinf=0.0, neginf=0.0)

# ---------------- CHECK 1: int sigma_kk^2 dA = 2 pi L^2 (L = 1) ----------------
def check1():
    R, dx = 4.0, 0.004
    n = int(2 * R / dx)
    x = (np.arange(n) + 0.5) * dx - R
    X, Y = np.meshgrid(x, x)
    sig_kk_sq = (2.0 * westergaard_Z(X, Y).real) ** 2
    inner = float(sig_kk_sq.sum()) * dx * dx
    tail = np.pi / (2.0 * R * R)   # |Z|~1/(2r^2) => int_R^inf 4/(4r^4) 2pi r dr
    total = inner + tail
    print("CHECK 1 (reciprocity consistency, long-time limit)")
    print(f"  int sigma_kk^2 dA = {total:.6f}  (grid {inner:.6f} + tail {tail:.6f})")
    print(f"  expected 2*pi*L^2 = {2*np.pi:.6f},  rel.err = "
          f"{abs(total-2*np.pi)/(2*np.pi):.3e}")
    return total

# ---------------- CHECK 2: short-time coefficient from 2-D slit diffusion ----------------
def check2():
    """Solve dq/dt = Laplacian q (c=1), q=0 initially, q=-1 held on the slit
    y=0, |x|<=1 (p0 = 1). Compute M(t) = int sigma_kk^R q dA.
    Boundary-layer/reciprocity prediction: M ~ 16 p0 L sqrt(t/pi)."""
    Ldom, dx = 1.4, 0.003
    n = int(2 * Ldom / dx)
    dt = 0.20 * dx * dx          # explicit FTCS stability: dt <= dx^2/4
    x = (np.arange(n) + 0.5) * dx - Ldom
    X, Y = np.meshgrid(x, x)
    sig_kk = 2.0 * westergaard_Z(X, Y).real
    q = np.zeros((n, n))
    j0 = n // 2                  # grid row at the crack line y = 0
    slit = (np.abs(x) <= 1.0)
    q[j0, slit] = -1.0
    t = 0.0
    t_end, nprint = 1.2e-2, 40
    marks = np.geomspace(2.6e-4, t_end, nprint)
    im = 0
    nsteps = int(t_end / dt)
    pred = 16.0 / np.sqrt(np.pi)
    ts, Ms = [], []
    print("CHECK 2 (2-D slit diffusion -> short-time volume growth)")
    print(f"  face-layer prediction: M_face = 16 sqrt(t/pi)  "
          f"(M/sqrt(t) -> {pred:.6f});  tips/curvature add O(t^0.75), O(t)")
    for k in range(1, nsteps + 1):
        lap = (np.roll(q, 1, 0) + np.roll(q, -1, 0)
               + np.roll(q, 1, 1) + np.roll(q, -1, 1) - 4.0 * q) / (dx * dx)
        q += dt * lap
        q[0, :] = q[-1, :] = 0.0
        q[:, 0] = q[:, -1] = 0.0
        q[j0, slit] = -1.0
        t = k * dt
        if im < nprint and t >= marks[im]:
            M = float((sig_kk * q).sum() * dx * dx)
            ts.append(t); Ms.append(M)
            if im % 8 == 0:
                Mf = 16.0 * np.sqrt(t / np.pi)
                print(f"  t={t:10.3e}  M={M:9.4f}  M_face={Mf:9.4f}  "
                      f"M/M_face={M/Mf:8.4f}")
            im += 1
    # least squares: M = C0 sqrt(t/pi) + C1 t^0.75 + C2 t ;  C0 should be ~ 16
    T = np.array(ts); Mv = np.array(Ms)
    B = np.column_stack([np.sqrt(T / np.pi), T ** 0.75, T])
    C0, C1, C2 = np.linalg.lstsq(B, Mv, rcond=None)[0]
    rms = np.sqrt(np.mean((Mv - B @ [C0, C1, C2]) ** 2))
    print(f"  fit M = C0 sqrt(t/pi) + C1 t^0.75 + C2 t over {len(T)} samples:")
    print(f"  C0 = {C0:.4f} (predicted 16), C1 = {C1:.4f}, C2 = {C2:.4f}, "
          f"rms resid = {rms:.2e}")
    return

print("=== MP-24 report values (part c: nu=0.2, nu_u=0.4, alpha=0.89) ===")
print(f"K1_0        = {K1_0:.10f}   [sqrt(pi L)]")
print(f"K1_inf      = {K1_inf:.10f}   [sqrt(pi L)]")
print(f"V1_0        = {V1_0:.10f}   [L^2/G]")
print(f"V1_inf      = {V1_inf:.10f}   [L^2/G]")
print(f"V1_short_coeff (A) = {A:.10f}   [L^2/G],  V ~ V1_0 + A sqrt(c t/L^2)")
print(f"p_faces_0   = {p_faces_0:.10f}   [per unit applied traction]")
print(f"K2_0        = {K2_0:.10f}   [sqrt(pi L)]")
print(f"K2_inf      = {K2_inf:.10f}   [sqrt(pi L)]")
print(f"V2_inf      = {V2_inf:.10f}   [L^2/G]")
print()
check1()
print()
check2()

results = {
    "K1_0": K1_0,
    "K1_inf": K1_inf,
    "V1_0": V1_0,
    "V1_inf": V1_inf,
    "V1_short_coeff": A,
    "p_faces_0": p_faces_0,
    "K2_0": K2_0,
    "K2_inf": K2_inf,
    "V2_inf": V2_inf,
    "formulas": {
        "K1_0": "sigma0*sqrt(pi*L)  [= sqrt(pi*L) per unit loading; time-independent for loading 1]",
        "K1_inf": "sigma0*sqrt(pi*L)",
        "V1_0": "sigma0*pi*(1-nu_u)*L^2/G   (undrained crack volume)",
        "V1_inf": "sigma0*pi*(1-nu)*L^2/G   (drained crack volume)",
        "V1_short_coeff": "A = sigma0*16*(nu_u-nu)/sqrt(pi) * L^2/G  in  V ~ V1_0 + A*sqrt(c*t/L^2)",
        "p_faces_0": "sigma0*2*(nu_u-nu)/(alpha*(1-2*nu))  (undrained pore pressure at the faces)",
        "K2_0": "0  (no traction applied; undrained response is identically zero)",
        "K2_inf": "0  (t->inf state: p = 1 uniform, sigma = 0; crack neither opens nor loads tips)",
        "V2_inf": "0  (uniform strain alpha*(1-2*nu)/(2G) everywhere; zero-thickness slit does not open)"
    },
    "reproduce_command": "python3 solve.py"
}
with open("results.json", "w") as f:
    json.dump(results, f, indent=2)
print()
print("results.json written")
