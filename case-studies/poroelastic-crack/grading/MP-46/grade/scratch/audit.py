"""Local grading audit; leaves the original submission untouched."""
from pathlib import Path
import json
import math
import runpy
import numpy as np

root = Path(__file__).resolve().parent.parent
submission = root / "submission"
copy = root / "scratch" / "reproduction"
reported = json.loads((submission / "results.json").read_text())
rerun = json.loads((copy / "results.json").read_text())
print("Reproduction results.json exactly matches original:", reported == rerun)
print("FEM output exactly matches delivered fem_out.txt:",
      (copy / "fem_rerun.txt").read_text() == (submission / "fem_out.txt").read_text())

key = {
    "K1_0": 0.807593,
    "K1_inf": 1.0,
    "V1_0": 0.6 * math.pi,
    "V1_inf": 0.8 * math.pi,
    "V1_short_coeff": 2.440220,
    "p_faces_0": 0.4 / (0.89 * 0.6),
    "K2_0": 0.0,
    "K2_inf": -0.33375,
    "V2_inf": -0.89 * math.pi * 0.6 / 2,
}
for name, expected in key.items():
    actual = reported[name]
    tolerance = 0.01 if expected == 0 else 0.01 * abs(expected)
    print(f"{name}: actual={actual:.12g}, key={expected:.12g}, "
          f"within_band={abs(actual - expected) <= tolerance}")

# Only instrument the scratch copy's return value, retaining the original
# assembly, boundaries, and solve. The ordinary rerun above used no changes.
source = (copy / "fem_poro.py").read_text()
assert source.count("    return log\n") == 1
instrumented = root / "scratch" / "fem_instrumented.py"
instrumented.write_text(source.replace("    return log\n", "    return log, un, pn\n"))
fem = runpy.run_path(str(instrumented), run_name="grading_audit")
nu, nu_u, alpha, G = (fem[x] for x in ("nu", "nu_u", "alpha", "G"))
beta = alpha**2 * (1 - 2*nu) * (1 - 2*nu_u) / (2*G*(nu_u - nu))
print("Constitutive coupling su: implemented", fem["su"], "required", alpha)
print("Constitutive storage bcoef: implemented", fem["bcoef"], "required", beta)

xy, tris, grad = (fem[x] for x in ("xy", "tris", "GRAD"))
for e in (0, 1):
    points = xy[tris[e]]
    # Use differences from the first node to avoid cancellation at x=-3.
    gx = (points[:, 0] - points[0, 0]) @ grad[e]
    gy = (points[:, 1] - points[0, 1]) @ grad[e]
    print(f"Element {e}: computed grad(x)={gx}, grad(y)={gy}; "
          "required [1,0] and [0,1]")

log, u, p = fem["run"](1, elastic=True)
_, V, Kd, Ku, center_opening = log[-1]
print("Elastic calibration dimensional K:", Kd,
      "K normalized by sqrt(pi L):", Kd / math.sqrt(math.pi*fem["L"]))
for i in (fem["iL"], fem["iR"]):
    w_tip = u[2*fem["up"][i]+1] - u[2*fem["dn"][i]+1]
    print("Opening at nominal crack tip x=", fem["xs"][i], ":", w_tip,
          "; continuous crack tips require zero")

Cp, Mp, rowsum = (fem[x] for x in ("Cp", "Mp", "rowsum_p"))
consistent = Cp @ u
implemented = Mp @ (consistent / rowsum)
print("Previous-step displacement contribution relative discrepancy:",
      np.linalg.norm(implemented-consistent) / np.linalg.norm(consistent))
print("Implemented RHS is su*Mp*M_lump^(-1)*Cp*u_n; "
      "consistent backward-Euler RHS is su*Cp*u_n (plus pressure storage).")

print("Reported transient V first/last ratio:", 0.93640 / 1.64975)
print("Loading 1 K last/first:", 1.43016 / 0.87947)
print("Loading 2 |K| last/first:", 0.38227 / 0.15114)
print("Permeability is included in Kp, then multiplied again by kap in A.")
