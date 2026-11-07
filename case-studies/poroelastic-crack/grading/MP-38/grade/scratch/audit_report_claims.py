"""Compare each reported spectral table and fit to recomputed scratch data."""
import importlib.util
import json
import pathlib
import re

ROOT = pathlib.Path(__file__).resolve().parent.parent
spec = importlib.util.spec_from_file_location("scratch_solver", ROOT / "scratch/reproduction/solve.py")
solver = importlib.util.module_from_spec(spec)
spec.loader.exec_module(solver)
import numpy as np

records = json.loads((ROOT / "scratch/report_checks.json").read_text())
if "production_160000" not in records:
    production = solver.spectral(80, 160000, .125, (80, 320, 1280, 5120))
    records["production_160000"] = {"n": 80, "cutoff": 160000, "spacing": .125, "values": production.tolist()}
    (ROOT / "scratch/report_checks.json").write_text(json.dumps(records, indent=2)+"\n")
answer = (ROOT / "submission/ANSWER.md").read_text()
names = ["uncoupled", "drained", "basis_48", "basis_96", "basis_160", "volume_cutoff_120000", "volume_spacing_0125", "tip_160_80000", "tip_96_120000", "tip_48_120000", "tip_96_300000", "tip_64_200000", "production_160000", "production_320000"]
tables = re.findall(r"\[\[(.*?)\]\]", answer, re.S)
assert len(tables) == len(names), (len(tables), len(names))
audit = {}
for name, table in zip(names, tables):
    reported = np.array([[float(x) for x in line.replace("[", "").replace("]", "").split()] for line in table.splitlines() if line.strip()])
    actual = np.array(records[name]["values"])[:, :reported.shape[1]]
    matches = bool(np.allclose(actual, reported, rtol=1e-8, atol=1e-8))
    audit[name] = {"printed_table_matches": matches, "max_absolute_difference": float(np.max(np.abs(actual-reported)))}
    assert matches, name

fit_specs = [
    ("volume_cutoff_120000", 2, 4, 2, [2.440061382113711, 2.440124224292605, 2.438870923602375]),
    ("volume_spacing_0125", 2, 4, 2, [2.4403607284375415, 2.440645750314374, 2.43940756974985]),
    ("tip_48_120000", 3, 3, 2, [.8075495387906676, .8064007065492166, .771111773897171]),
    ("tip_96_300000", 3, 3, 2, [.8075496261092631, .8075892226710333, .806460525912287]),
    ("tip_64_200000", 3, 3, 2, [.8075431965629399]),
    ("tip_64_200000", 2, 4, 3, [2.4398607431830546]),
    ("production_160000", 3, 3, 2, [.8075897671067087]),
    ("production_160000", 2, 4, 3, [2.4400664995666825]),
    ("production_320000", 3, 3, 2, [.807587689880131]),
    ("production_320000", 2, 4, 3, [2.440069858483082]),
]
for name, column, window, degree, expected in fit_specs:
    values = np.array(records[name]["values"])
    x = 1/np.sqrt(values[:, 0])
    start = 1 if name == "tip_64_200000" and column == 3 else 0
    actual = [float(np.polyfit(x[i+start:i+start+window], values[i+start:i+start+window, column], degree)[-1]) for i in range(len(expected))]
    matches = bool(np.allclose(actual, expected, rtol=0, atol=1e-10))
    audit[f"{name}_column{column}_window{window}_degree{degree}"] = {"reported": expected, "recomputed": actual, "matches": matches}
    assert matches, (name, actual, expected)
(ROOT / "scratch/report_checks.json").write_text(json.dumps(records, indent=2)+"\n")
(ROOT / "scratch/report_claims_audit.json").write_text(json.dumps(audit, indent=2)+"\n")
print(f"All {len(tables)} reported spectral tables and {len(fit_specs)} groups of fitted values reproduce.")
print("Production fits:", audit["production_160000_column3_window3_degree2"]["recomputed"], audit["production_160000_column2_window4_degree3"]["recomputed"])
