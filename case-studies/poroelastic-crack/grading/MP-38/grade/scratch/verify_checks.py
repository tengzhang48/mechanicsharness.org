"""Recreate the submission's reported checks using the scratch copy only."""
import importlib.util
import json
import math
import pathlib
import time

ROOT = pathlib.Path(__file__).resolve().parent.parent
spec = importlib.util.spec_from_file_location("scratch_solver", ROOT / "scratch/reproduction/solve.py")
solver = importlib.util.module_from_spec(spec)
spec.loader.exec_module(solver)
import numpy as np
from scipy.integrate import quad

records = {}

def run(name, n, cutoff, spacing, qs, nu=.2, nu_u=.4):
    start = time.monotonic()
    values = solver.spectral(n, cutoff, spacing, qs, nu=nu, nu_u=nu_u)
    record = {"n": n, "cutoff": cutoff, "spacing": spacing, "values": values.tolist(), "elapsed_seconds": time.monotonic()-start}
    if len(qs) >= 3:
        x = 1/np.sqrt(values[:, 0])
        record["initial_tip_first_three_quadratic"] = float(np.polyfit(x[:3], values[:3, 3], 2)[-1])
        record["short_coefficient_first_four_cubic"] = float(np.polyfit(x[:4], values[:4, 2], 3)[-1]) if len(qs) >= 4 else None
        record["short_coefficient_sliding_quadratics"] = [float(np.polyfit(x[i:i+3], values[i:i+3, 2], 2)[-1]) for i in range(max(0, len(qs)-3))]
    records[name] = record
    (ROOT / "scratch/report_checks.json").write_text(json.dumps(records, indent=2)+"\n")
    print(json.dumps({"case": name, "elapsed_seconds": record["elapsed_seconds"], "initial_tip_first_three_quadratic": record.get("initial_tip_first_three_quadratic"), "short_coefficient_first_four_cubic": record.get("short_coefficient_first_four_cubic"), "short_coefficient_sliding_quadratics": record.get("short_coefficient_sliding_quadratics")}), flush=True)
    return values

no_coupling = run("uncoupled", 16, 4000, .125, (1, 100), nu_u=.2)
assert np.allclose(no_coupling[:, 1], .8*math.pi, atol=1e-12)
assert np.allclose(no_coupling[:, 2], 0, atol=1e-12)
assert np.allclose(no_coupling[:, 3], 1, atol=1e-12)
drained = run("drained", 48, 40000, .125, (.001, .01, .1, 1, 4))
assert abs(drained[0, 1]/(.8*math.pi)-1) < .01
assert abs(drained[0, 3]-1) < .01

for n in (48, 96, 160):
    run(f"basis_{n}", n, 40000, .25, (100, 200, 400, 800, 1600, 3200))
run("volume_cutoff_120000", 96, 120000, .25, (400, 800, 1600, 3200, 6400, 12800))
run("volume_spacing_0125", 160, 40000, .125, (400, 800, 1600, 3200, 6400, 12800))
for name, n, cutoff, spacing, qs in (
    ("tip_160_80000", 160, 80000, .125, (20, 80, 320, 1280, 5120)),
    ("tip_96_120000", 96, 120000, .125, (20, 80, 320, 1280, 5120)),
    ("tip_48_120000", 48, 120000, .25, (20, 80, 320, 1280, 5120)),
    ("tip_96_300000", 96, 300000, .5, (20, 80, 320, 1280, 5120)),
    ("tip_64_200000", 64, 200000, .125, (20, 80, 320, 1280)),
    ("production_320000", 80, 320000, .125, (80, 320, 1280, 5120)),
):
    run(name, n, cutoff, spacing, qs)

production = records["production_320000"]
assert abs(production["initial_tip_first_three_quadratic"]-.807587689880131) < 1e-7
assert abs(production["short_coefficient_first_four_cubic"]-2.440069858483082) < 1e-7

eta = .89*.6/(2*.8)
v1 = quad(lambda x: 2*.8*math.sqrt(1-x*x), -1, 1, epsabs=1e-12)[0]
v2 = quad(lambda x: -2*.8*eta*math.sqrt(1-x*x), -1, 1, epsabs=1e-12)[0]
records["opening_integrals"] = {"V1_integrated": v1, "V1_expected": .8*math.pi, "V2_integrated": v2, "V2_expected": -.89*.6*math.pi/2}
assert abs(v1-.8*math.pi) < 1e-12
assert abs(v2+.89*.6*math.pi/2) < 1e-12
(ROOT / "scratch/report_checks.json").write_text(json.dumps(records, indent=2)+"\n")
print("Reported convergence configurations and opening-integral checks completed.", flush=True)
