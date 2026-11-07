"""Independent checks for submission C. Does not write into C/."""
import json
import os
import sys
import time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from crack_solver2 import Params, solve

nu, nuu, al = 0.2, 0.4, 0.89
P = Params(nu, nuu, al)
eta = al * (1 - 2 * nu) / (2 * (1 - nu))
beta = (nuu - nu) / (al * (1 - 2 * nu))
gam = (nuu - nu) / (1 - nu)
S = al**2 * (1 - 2 * nu) ** 2 * (1 - nuu) / (2 * (nuu - nu) * (1 - nu))  # G=1

ref = {
    "K1_0": 0.807593,
    "K1_inf": 1.0,
    "V1_0": np.pi * (1 - nuu),
    "V1_inf": np.pi * (1 - nu),
    "V1_short_coeff": 2.440220,
    "p_faces_0": 2 * (nuu - nu) / (al * (1 - 2 * nu)),
    "K2_0": 0.0,
    "K2_inf": -eta,
    "V2_inf": -np.pi * al * (1 - 2 * nu) / 2.0,
}
paper = {
    "K1_0": (1 - nuu) / (1 - nu),
    "V1_short_coeff_faces_only": 16 * gam * (1 - nuu) / np.sqrt(np.pi),
}

with open(os.path.join(HERE, "results.json")) as f:
    reported = json.load(f)

print("=== closed-form identities ===")
print(f"eta={eta:.12f} beta={beta:.12f} gamma={gam:.12f} S={S:.12f}")
print(f"2*beta={2*beta:.12f}  gamma==2*beta*eta? {abs(gam - 2*beta*eta):.3e}")
print(f"V1_inf/V1_0={(1-nu)/(1-nuu):.12f}")
A_face = 16 * gam * (1 - nuu) / np.sqrt(np.pi)
A_tip_fo = 4 * np.sqrt(np.pi) * gam * (1 - nuu)
A_fo = 4 * (4 + np.pi) * gam * (1 - nuu) / np.sqrt(np.pi)
print(f"A_face exact={A_face:.12f}")
print(f"A_tip first-order={A_tip_fo:.12f}")
print(f"A first-order closed={A_fo:.12f}")
print(f"paper K1_0 heuristic={(1-nuu)/(1-nu):.12f}")
print(f"J upper bound sqrt((1-nu_u)/(1-nu))={np.sqrt((1-nuu)/(1-nu)):.12f}")
print(f"Yu-like 0.75**0.743={((1-nuu)/(1-nu))**0.743:.12f}")
print(f"log-capacity D=4*ln2-2*gamma_E={4*np.log(2)-2*np.euler_gamma:.12f}")

print("\n=== nine-value comparison (1% nonzero, 0.01 for K2_0) ===")
rows = []
for k, rref in ref.items():
    rval = reported[k]
    if k == "K2_0":
        ok = abs(rval - rref) <= 0.01
        err = abs(rval - rref)
        print(f"{k:16s} reported={rval:.8g}  ref={rref:.8g}  abs={err:.3e}  {'PASS' if ok else 'FAIL'}")
    else:
        err = abs(rval - rref) / abs(rref)
        ok = err <= 0.01
        print(f"{k:16s} reported={rval:.8g}  ref={rref:.8g}  rel={err:.3e}  {'PASS' if ok else 'FAIL'}")
    rows.append((k, rval, rref, err, ok))

print("\n=== principle checks on agent's own numbers ===")
print(f"V1_inf/V1_0={reported['V1_inf']/reported['V1_0']:.12f}  expected={(1-nu)/(1-nuu):.12f}")
print(f"K1_inf={reported['K1_inf']}")

print("\n=== re-extrapolation from saved large-s data ===")
data = {}
for name in ["large_s_1e+02.json", "large_s_1e+06.json", "large_s_1e+07.json"]:
    for r in json.load(open(os.path.join(HERE, name))):
        if r["s"] not in data or r["N"] > data[r["s"]]["N"]:
            data[r["s"]] = r
extra = {
    100.0: (0.8650232778471896, 2.028587164286242, -0.0734811080445404, -0.11336879472227795),
    1e3: (0.83963065, 1.93942880, -0.04170787, -0.04167034),
    1e4: (0.82546865, 1.90398655, -0.02352233, -0.01433188),
    1e5: (0.81759171, 1.89131882, -0.01323970, -0.00475444),
}
for s, v in extra.items():
    if s not in data:
        data[s] = dict(s=s, sK1=v[0], sV1=v[1], sK2=v[2], sV2=v[3])

ss = np.array(sorted(data))
sK1 = np.array([data[s]["sK1"] for s in ss])
sV1 = np.array([data[s]["sV1"] for s in ss])
sK2 = np.array([data[s]["sK2"] for s in ss])
sV2 = np.array([data[s]["sV2"] for s in ss])
b2 = 2 * P.beta
Vu = np.pi * (1 - P.nuu)
W = sK1 + b2 * sK2
VW = sV1 + b2 * sV2
print(" s        sK1        W            (sV1-Vu)rt(s)  (VW-Vu)rt(s)   sV2*rt(s)   s^{1/4}sK2")
for i, s in enumerate(ss):
    print(
        f"{s:8.0e} {sK1[i]:.10f} {W[i]:.10f} {(sV1[i]-Vu)*np.sqrt(s):.8f} "
        f"{(VW[i]-Vu)*np.sqrt(s):.8f} {sV2[i]*np.sqrt(s):.8f} {sK2[i]*s**0.25:.8f}"
    )

def fit(x, y, powers):
    A = np.vstack([x ** (-p) for p in powers]).T
    coef, residuals, rank, sv = np.linalg.lstsq(A, y, rcond=None)
    pred = A @ coef
    return coef, pred, y - pred

sel = ss >= 1e3
a_coef, a_pred, a_res = fit(ss[sel], W[sel], [0, 0.5, 1.0])
pw = [0, 0.25, 0.5, 0.75] if sel.sum() >= 5 else [0, 0.25, 0.5]
cW_coef, _, _ = fit(ss[sel], (VW[sel] - Vu) * np.sqrt(ss[sel]), pw)
c2_coef, _, _ = fit(ss[sel], sV2[sel] * np.sqrt(ss[sel]), pw)
c1_coef, _, _ = fit(ss[sel], (sV1[sel] - Vu) * np.sqrt(ss[sel]), pw)
cK2_coef, _, _ = fit(ss[sel], sK2[sel] * ss[sel] ** 0.25, [0, 0.5, 1.0])
print("W residuals at selected s:", a_res)
print(f"extrapolated a={a_coef[0]:.12f}  claimed 0.807592986")
print(f"max |W residual|={np.max(np.abs(a_res)):.3e}")
print(f"W(s=1e7)-a = {W[ss==1e7][0]-a_coef[0]:.3e}")
cW = cW_coef[0]
A_tip = 2 / np.sqrt(np.pi) * cW
A_direct = 2 / np.sqrt(np.pi) * c1_coef[0]
A_face_num = -b2 * 2 / np.sqrt(np.pi) * c2_coef[0]
print(f"cW={cW:.12f} chat={cW/(1-P.nuu):.12f} A_tip={A_tip:.12f}")
print(f"A_face numeric={A_face_num:.12f}  formula={A_face:.12f}  diff={A_face_num-A_face:.3e}")
print(f"A_direct={A_direct:.12f}  A_exactface+tip={A_face+A_tip:.12f}")
print(f"claimed uncertainty a +/- 5e-6, A +/- 3e-4")
print(f"|a-0.807593|={abs(a_coef[0]-0.807593):.3e}")
print(f"|A_exactface+tip-2.440220|={abs(A_face+A_tip-2.440220):.3e}")

print("\n=== weak-coupling a(gamma) from saved a_of_gamma files ===")
import glob
rows_g = []
for f in sorted(glob.glob(os.path.join(HERE, "a_of_gamma_*.json"))):
    rows_g += json.load(open(f))
rows_g.sort(key=lambda r: r["gamma"])
print(" gamma    a_fit      1-3g/4    (1-a)/g   chat      2*pi*g   chat/(2pi g)")
for r in rows_g:
    g = r["gamma"]
    d = r["data"]
    s = np.array([x["s"] for x in d])
    Wg = np.array([x["W"] for x in d])
    cWg = np.array([x["cW"] for x in d])
    A = np.vstack([np.ones(3), s ** -0.5, s ** -1.0]).T
    a = np.linalg.solve(A, Wg)[0]
    B = np.vstack([np.ones(3), s ** -0.25, s ** -0.5]).T
    ch = np.linalg.solve(B, cWg)[0]
    print(f"{g:6.3f}  {a:.8f}  {1-0.75*g:.8f}  {(1-a)/g:.5f}  {ch:.5f}  {2*np.pi*g:.5f}  {ch/(2*np.pi*g):.5f}")

print("\n=== time-domain short/long checks from saved Talbot table ===")
td = json.load(open(os.path.join(HERE, "time_domain_M12.json")))
A_rep = reported["V1_short_coeff"]
K10 = reported["K1_0"]
C_k1 = 0.1947
C_k2 = 0.2599
for row in td:
    t = row["t"]
    if t <= 0.01:
        V_pred = reported["V1_0"] + A_rep * np.sqrt(t)
        K1_pred = K10 + C_k1 * t ** 0.25
        K2_pred = -C_k2 * t ** 0.25
        print(
            f"t={t:.4g}  V1={row['V1']:.6f} vs V0+A sqrt(t)={V_pred:.6f} dV={row['V1']-V_pred:.4e}  "
            f"K1={row['K1']:.6f} vs {K1_pred:.6f}  K2={row['K2']:.6f} vs {K2_pred:.6f}"
        )
tlast = td[-1]["t"]
row = td[-1]
K2_log = -eta * (1 - 1 / (np.log(16 * tlast) - np.euler_gamma))
V2_log = -np.pi * eta * (1 - nu) * (1 - 1 / (np.log(16 * tlast) - np.euler_gamma))
print(f"t={tlast:.0e} K1={row['K1']:.8f} V1={row['V1']:.8f} K2={row['K2']:.6f} (log {K2_log:.6f}) V2={row['V2']:.6f} (log {V2_log:.6f})")
print(f"K2(t=1e4)/K2_inf={row['K2']/(-eta):.4f}  still far from 1 (2D log approach)")

print("\n=== representative solver reruns (copied crack_solver2.py) ===")
# s=100, N=20 should match large_s_1e+02.json
t0 = time.time()
o100 = solve(100.0, 20, P)
dt = time.time() - t0
saved100 = json.load(open(os.path.join(HERE, "large_s_1e+02.json")))[0]
print(f"s=100 N=20 ({dt:.2f}s)")
for name, keyK, keyV in [("L1", "sK1", "sV1"), ("L2", "sK2", "sV2")]:
    sK = 100.0 * o100[name]["K"]
    sV = 100.0 * o100[name]["V"]
    print(
        f"  {name} sK={sK:.12f} saved={saved100[keyK]:.12f} d={sK-saved100[keyK]:.3e}  "
        f"sV={sV:.12f} saved={saved100[keyV]:.12f} d={sV-saved100[keyV]:.3e}"
    )

# s=1, N=24: drained-ish not yet; energy identity
t0 = time.time()
o1 = solve(1.0, 24, P)
h = 1e-4
op = solve(1.0 + h, 24, P)["L1"]["V"]
om = solve(1.0 - h, 24, P)["L1"]["V"]
dV = (op - om) / (2 * h)
Jen = o1["L1"]["V"] / 1.0 + 0.5 * dV
Jmech = (1 - P.nu) * np.pi * o1["L1"]["K"] ** 2 / 2
Jfl = np.pi * np.sum(o1["L1"]["b"]) ** 2 / (2 * P.kappa * 1.0)
print(f"s=1 N=24 ({time.time()-t0:.2f}s) sK1={o1['L1']['K']:.8f} sV1={o1['L1']['V']:.8f} sK2={o1['L2']['K']:.8f} sV2={o1['L2']['V']:.8f}")
print(f"  energy J_en={Jen:.10f} J_mech+J_fl={Jmech+Jfl:.10f} rel={abs(Jen-Jmech-Jfl)/abs(Jen):.3e}")
print(f"  fluid fraction of J={Jfl/(Jmech+Jfl):.4f}")

# s=1e-4, N=16: long-time loading 1 should be essentially drained; loading 2 not
t0 = time.time()
o4 = solve(1e-4, 16, P)
print(f"s=1e-4 N=16 ({time.time()-t0:.2f}s)")
print(
    f"  sK1={1e-4*o4['L1']['K']:.10f} (expect ~1)  sV1={1e-4*o4['L1']['V']:.10f} (expect {np.pi*(1-nu):.10f})"
)
print(
    f"  sK2={1e-4*o4['L2']['K']:.10f} (inf {-eta:.6f})  sV2={1e-4*o4['L2']['V']:.10f} (inf {-np.pi*eta*(1-nu):.6f})"
)
print("  matches small_s.log line 1: sK1=0.999976784 sV1=2.513213157 sK2=-0.302926321 sV2=-0.761333698")

# s=10, N=20: intermediate, compare with ANSWER table at tau~ corresponding Laplace not equal, skip
print("\n=== face-condition / kernel sanity ===")
# drained large-k: 1-gamma factor
print(f"(1-gamma)/(1-nu_u) should equal 1/(1-nu): {(1-gam)/(1-nuu):.12f} vs {1/(1-nu):.12f}")
print(f"solver uses drained K conversion 1/(1-nu)={1/(1-nu)}")
print(f"loading 1 rhs: sigma=-1/s, p=0; loading 2: sigma=0, p=1/s")
print(f"W loading (sigma=-1, p=2 beta) cancels face layer: 2*beta={b2:.12f}")

print("\nDONE")
