"""FINAL production run: computes every reported quantity with full audit.

Exact inputs (analytic):
  eps  = ln((1-beta)/(1+beta))/(2 pi) = ln(2/3)/(2 pi), beta = 0.2, alpha = 0.8
  G    = beam/J-integral value (exact; far-ahead stress-free, far-behind
         Saint-Venant states are exact elasticity solutions)
  G-K relation: G = |K|^2 (1/Eb1 + 1/Eb2)/(2 cosh^2(pi eps))
         [= |K|^2 * 64/135 * E2^-1 for this pair], verified to 5e-7 relative
         by the J-integral of the crack-tip eigenpair (eigenpair.py/check_c.py)
FEM: cases A and B at several mesh levels; K from traction fit and from
crack-opening fit (with the eigenpair-verified opening factor Q); Richardson
extrapolation of phases; superposition for the combined case.
"""
import numpy as np
import json, sys
sys.path.insert(0, '/work/code')
from fem import (BilayerCrack, plane_strain_D, end_loads, solve,
                 interface_tractions, crack_opening)
from fits import fit_series

# ---------------- exact constants ----------------
nu = 1/3
E1, E2 = 9.0, 1.0
h, H = 1.0, 2.0
Eb1, Eb2 = E1/(1-nu**2), E2/(1-nu**2)     # 10.125, 1.125
S = 1/Eb1 + 1/Eb2                          # 80/81
alpha = (Eb1-Eb2)/(Eb1+Eb2)
beta = (alpha/4)                           # = (1-2nu)/(2(1-nu)) * alpha, nu=1/3
eps = np.log((1-beta)/(1+beta))/(2*np.pi)   # -0.0645318
cGK = S/2/np.cosh(np.pi*eps)**2             # 64/135 = 0.474074...
G_A = 1/(2*Eb1*h) + 1/(2*Eb2*H) + 6*((h+H)/2)**2/(Eb2*H**3)   # P=1
G_B = 6*(1/(Eb1*h**3) + 1/(Eb2*H**3))                         # M=1
KA_abs = np.sqrt(G_A/cGK)
KB_abs = np.sqrt(G_B/cGK)
# eigenpair opening factor: delta = Q conj(K) r^{1/2 - i eps} (eps<0)
eps_m = -eps
Q = 1j*4*S/((1+2j*eps_m)*np.sqrt(2*np.pi)*np.cosh(np.pi*eps_m))

audit = []
def say(s):
    print(s)
    audit.append(s)

say(f"alpha = {alpha:.6f}, beta = {beta:.6f}, eps = {eps:.8f}")
say(f"sech^2(pi eps) = {1/np.cosh(np.pi*eps)**2:.8f}  (exact 24/25 = {24/25})")
say(f"c_GK = E2*S/2*sech^2 = {E2*cGK:.10f}  (exact 64/135 = {64/135})")
say(f"G_A E2 h/P^2 = {G_A*E2*h:.10f}   |K_A| sqrt(h)/P = {KA_abs:.8f}")
say(f"G_B E2 h^3/M^2 = {G_B*E2*h**3:.10f}   |K_B| h^3/2/M = {KB_abs:.8f}")
say(f"phase change per decade = eps*ln10 = {np.degrees(eps*np.log(10)):.6f} deg")

# ---------------- FEM ----------------
D = [plane_strain_D(E1, nu), plane_strain_D(E2, nu)]

def sigma(P, M):
    return {1: lambda x2: -P/h - 12*M*(x2 - h/2)/h**3,
            2: lambda x2:  P/H + 12*(M + P*(h+H)/2)*(x2 + H/2)/H**3}

def Kfits(P, M, d1, d2, L=30.0):
    mesh = BilayerCrack(h, H, L, d1, d2, q1=1.35, q2=1.5, cap2=0.35)
    f = end_loads(mesh, sigma(P, M))
    u = solve(mesh, D, f)
    rs = np.geomspace(4e-4, 8e-2, 60)
    T = interface_tractions(mesh, D, u, rs)
    outT = {}
    for (rmin, rmax) in [(5e-4, 1e-2), (1e-3, 2e-2), (2e-3, 4e-2)]:
        m = (rs >= rmin) & (rs <= rmax)
        y = T[m]*np.sqrt(2*np.pi*rs[m])*np.exp(-1j*eps*np.log(rs[m]))
        outT[(rmin, rmax)] = fit_series(rs[m], y, (0.0, 0.5, 1.0))[0]
    rd, dd = crack_opening(mesh, u)
    outD = {}
    for (rmin, rmax) in [(5e-4, 1e-2), (1e-3, 2e-2)]:
        m = (rd >= rmin) & (rd <= rmax)
        yd = dd[m]/np.sqrt(rd[m])*np.exp(+1j*eps*np.log(rd[m]))
        Ad = fit_series(rd[m], yd, (0.0, 0.5, 1.0))[0]
        outD[(rmin, rmax)] = np.conj(Ad/Q)
    return outT, outD

levels = [(1.2e-4, 3.0e-4), (6e-5, 1.5e-4), (3e-5, 7.5e-5),
          (1.5e-5, 3.75e-5), (7.5e-6, 1.875e-5), (4e-6, 1e-5)]
wkey = (1e-3, 2e-2)

results = {}
for case, (P, M) in [("A", (1.0, 0.0)), ("B", (0.0, 1.0))]:
    say(f"\n===== case {case}: FEM convergence (T-fit window {wkey})")
    seq = []
    fits = []
    for (d1, d2) in levels:
        outT, outD = Kfits(P, M, d1, d2)
        fits.append((outT, outD))
        KT = outT[wkey]
        KD = outD[(1e-3, 2e-2)]
        seq.append((d1, d2, KT, KD))
        say(f"  d={d1:.1e}/{d2:.1e}: K_T={KT:.6f} (|.|={abs(KT):.6f}, "
            f"ph={np.degrees(np.angle(KT)):+.4f})  K_D={KD:.6f} "
            f"(ph={np.degrees(np.angle(KD)):+.4f})")
        # other windows at the two finest levels
        if (d1, d2) in levels[-2:]:
            for k, v in outT.items():
                say(f"      T-window {k}: {v:.6f} ph={np.degrees(np.angle(v)):+.4f}")
    # Richardson on phase (last three levels, Aitken)
    ph = np.array([np.angle(s[2]) for s in seq[-3:]])
    ph2 = ph[-2] + (ph[-2]-ph[-3])*((ph[-1]-ph[-2])/((ph[-2]-ph[-3])
         - (ph[-1]-ph[-2]) + 1e-30) if abs((ph[-2]-ph[-3])-(ph[-1]-ph[-2])) > 1e-12 else 0)
    spread = max(np.degrees(ph.max()-ph.min()),
                 max(np.degrees(abs(np.angle(s[2])-np.angle(s[3]))) for s in seq[-2:]))
    results[case] = dict(seq=seq, phase_deg=float(np.degrees(ph2)),
                         spread_deg=float(spread), fits=fits)

phA = results["A"]["phase_deg"]
phB = results["B"]["phase_deg"]
say(f"\nphase A (Richardson) = {phA:+.4f} deg  (spread {results['A']['spread_deg']:.4f})")
say(f"phase B (Richardson) = {phB:+.4f} deg  (spread {results['B']['spread_deg']:.4f})")

# Consensus phases: median over the finest-mesh estimators
# (3 traction windows + delta-fit window at the two finest levels)
def consensus(case):
    vals = []
    for (outT, outD) in results[case]["fits"][-2:]:
        vals += [np.angle(v) for v in outT.values()]
        vals += [np.angle(v) for v in outD.values()]
    return float(np.degrees(np.median(vals)))

phA = consensus("A")
phB = consensus("B")
say(f"phase A (consensus median) = {phA:+.4f} deg")
say(f"phase B (consensus median) = {phB:+.4f} deg")

# combined M = P h: K_tot = K_A + K_B (superposition, verified in validations.py)
Ktot_fem = results["A"]["seq"][-1][2] + results["B"]["seq"][-1][2]
Ktot_exactmag = KA_abs*np.exp(1j*np.radians(phA)) + KB_abs*np.exp(1j*np.radians(phB))
say(f"combined phase (fem K_A+K_B)   = {np.degrees(np.angle(Ktot_fem)):+.4f} deg")
say(f"combined phase (exact |K|s)    = {np.degrees(np.angle(Ktot_exactmag)):+.4f} deg")

perdec = float(np.degrees(eps*np.log(10)))
out = {
    "caseA": {"G_E2_h_over_P2": float(G_A*E2*h),
              "amplitude_sqrt_h_over_P": float(KA_abs),
              "phase_at_h_deg": float(phA)},
    "caseB": {"G_E2_h3_over_M2": float(G_B*E2*h**3),
              "amplitude_h32_over_M": float(KB_abs),
              "phase_at_h_deg": float(phB)},
    "phase_at_h_both_deg": float(np.degrees(np.angle(Ktot_exactmag))),
    "phase_change_per_decade_deg": perdec,
}
say("\n" + json.dumps(out, indent=1))
with open('/work/results.json', 'w') as fh:
    json.dump(out, fh, indent=1)
with open('/work/code/final_audit.txt', 'w') as fh:
    fh.write("\n".join(audit) + "\n")
print("wrote /work/results.json and /work/code/final_audit.txt")
