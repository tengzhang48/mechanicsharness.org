"""Independent local grading checks; never modifies submission/."""
import importlib.util
import json
import math
from pathlib import Path
import re
import subprocess

import sympy as sp

root = Path(__file__).resolve().parents[1]
copy = root / "scratch" / "rerun"
run = subprocess.run(
    ["python3", "solve.py"], cwd=copy, capture_output=True, text=True, check=True
)
(root / "scratch" / "rerun_stdout.txt").write_text(run.stdout)
reported = json.loads((root / "submission" / "results.json").read_text())
rerun = json.loads((copy / "results.json").read_text())
keys = [
    "K1_0", "K1_inf", "V1_0", "V1_inf", "V1_short_coeff",
    "p_faces_0", "K2_0", "K2_inf", "V2_inf",
]
nu_value, nu_u_value, alpha_value = 0.2, 0.4, 0.89
reference = {
    "K1_0": 0.807593,
    "K1_inf": 1.0,
    "V1_0": math.pi * (1 - nu_u_value),
    "V1_inf": math.pi * (1 - nu_value),
    "V1_short_coeff": 2.440220,
    "p_faces_0": 2 * (nu_u_value - nu_value) / (alpha_value * (1 - 2 * nu_value)),
    "K2_0": 0.0,
    "K2_inf": -alpha_value * (1 - 2 * nu_value) / (2 * (1 - nu_value)),
    "V2_inf": -alpha_value * math.pi * (1 - 2 * nu_value) / 2,
}
print("Rerun JSON exactly equals delivered JSON:", rerun == reported)
printed = dict(re.findall(r"^\s+(\w+)\s*=\s*([-+0-9.eE]+)$", run.stdout, re.M))
print("All nine numbers match printed output to its rounding:", all(
    abs(float(printed[k]) - reported[k]) <= 5.1e-11 for k in keys
))
for key in keys:
    tolerance = 0.01 if reference[key] == 0 else 0.01 * abs(reference[key])
    ok = abs(reported[key] - reference[key]) <= tolerance
    print(f"{key}: reported={reported[key]:.10f}; key={reference[key]:.10f}; pass={ok}")
print("Normalized submitted K fields:", {
    k: reported[k] / math.sqrt(math.pi) for k in keys if k.startswith("K")
})

# Derive identities directly from the supplied two-dimensional constitutive laws.
nu, nu_u, alpha, p, trace = sp.symbols("nu nu_u alpha p trace")
sigma11 = sp.symbols("sigma11")
zeta_rhs = alpha * (1 - 2 * nu) * trace + alpha**2 * (1 - 2 * nu)**2 * p / (nu_u - nu)
p_undrained = sp.solve(zeta_rhs, p)[0]
strain11_rhs = sigma11 - nu * trace + alpha * (1 - 2 * nu) * p
print("Undrained compliance residual (correct target):", sp.simplify(
    strain11_rhs.subs(p, p_undrained) - (sigma11 - nu_u * trace)
))
trace_undrained = sp.solve(zeta_rhs, trace)[0]
effective_trace = sp.factor(trace_undrained + 2 * alpha * p)
claimed_effective_trace = alpha * (1 - 2 * nu_u) * p / (nu_u - nu)
print("Undrained effective-stress trace derived from laws:", effective_trace)
print("Residual of solver's claimed positive-sign trace identity:", sp.factor(
    effective_trace - claimed_effective_trace
))
L, G, q = sp.symbols("L G q", positive=True)
x = sp.symbols("x", real=True)
opening = 2 * q * (1 - nu) * sp.sqrt(L**2 - x**2) / G
volume = sp.integrate(opening, (x, -L, L))
print("Ellipse volume integration residual:", sp.simplify(volume - sp.pi * q * (1 - nu) * L**2 / G))

spec = importlib.util.spec_from_file_location("copied_submission", copy / "solve.py")
solver = importlib.util.module_from_spec(spec)
spec.loader.exec_module(solver)
equal_moduli = solver.results_core(nu=0.2, nu_u=0.2, alpha=0.89)
print("Solver equal-moduli short-time growth (must vanish):", equal_moduli["V1_short_coeff"])
print("Solver equal-moduli induced pressure (must vanish):", equal_moduli["p_faces_0"])
print("Full hydraulic loading: solver normalized long-time K =", (
    reported["K1_inf"] + reported["K2_inf"]
) / math.sqrt(math.pi), "; own claimed anchor = 1; key superposition =", 1 + reference["K2_inf"])
