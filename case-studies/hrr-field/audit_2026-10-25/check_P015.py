"""Automatic numeric checker for P015 T1 and T2 (audit draft; the source problems have no checks.py).

Bands copied from the rubrics:
  T1 (T1/RUBRIC.md:7,9,10): p full 1e-5 / partial 1e-3; H full 0.005 / partial 0.05; C full 0.002 / partial 0.02.
  T2 (T2/RUBRIC.md:9,10):   p full 3e-4 / partial 3e-3; R full 0.002 / partial 0.02.
Key values: builder's reference_results.json, reproduced independently in recompute_key_results.json
(T1 H, C to 2e-9; T2 p, R to 1e-12).
Usage: python -I check_P015.py results.json   (results: {"T1": {"p","H","C"}, "T2": {"p","R"}}; either part optional)
       python -I check_P015.py --selftest      (key, saved runs, mutants -> check_P015_selftest.json)
"""
import json
import math
import sys
from pathlib import Path

KEY = {"T1": {"p": 2 / 15, "H": 1.893569725213, "C": 0.457305783974},
       "T2": {"p": 0.156175260286, "R": 0.373262313312}}
BANDS = {"T1": {"p": (1e-5, 1e-3), "H": (0.005, 0.05), "C": (0.002, 0.02)},
         "T2": {"p": (3e-4, 3e-3), "R": (0.002, 0.02)}}


def check(res):
    out = {}
    for task, items in BANDS.items():
        if task not in res:
            continue
        out[task] = {}
        for q, (full, part) in items.items():
            v = res[task].get(q)
            if v is None or (isinstance(v, float) and math.isnan(v)):
                out[task][q] = "missing"
                continue
            e = abs(v - KEY[task][q])
            out[task][q] = "full" if e <= full else ("partial" if e <= part else "fail")
        out[task]["all_full"] = all(x == "full" for k, x in out[task].items() if k != "all_full")
    return out


def selftest(here):
    mut = json.loads((here / "mutant_values.json").read_text())
    cases = []
    cases.append(("key", {"T1": KEY["T1"], "T2": KEY["T2"]}, "all full", None))
    # saved runs (numbers copied from the ANSWER.md that was graded; file:line in AUDIT.md)
    runs = [
        ("run T1 with-glm (R11)", {"T1": {"p": 0.133333, "H": 1.8936, "C": 0.4573}}, "pass", "p 10/10, H 20/20, C 15/15"),
        ("run T1 with-dsaz (R35)", {"T1": {"p": 0.1333333333, "H": 1.8935697252, "C": 0.4573057839}}, "pass", "10, 20, 15"),
        ("run T1 without-dsaz graded copy (R31)", {"T1": {"p": 0.1333333, "H": 2.15092, "C": 0.29178}}, "fail H, C", "p 10/10, H 0/20, C 0/15"),
        ("run T1 without-dsaz post-deadline overwrite (runs/ folder, not graded)", {"T1": {"p": 0.1333333, "H": 1.4706, "C": 1.6978}}, "fail H, C", None),
        ("run T2 with-dsaz (R32)", {"T2": {"p": 0.15618, "R": 0.3733}}, "pass", "p 15/15, R 15/15"),
        ("run T2 without-glm (R09)", {"T2": {"p": 0.1562, "R": 0.3733}}, "pass", "p 15/15, R 15/15"),
        ("run T2 without-dsaz (R23)", {"T2": {"p": 0.08764, "R": 0.3330}}, "fail both", "p 4/15 (procedure), R 4/15 (procedure)"),
    ]
    cases += runs
    for name, v in mut["T1"].items():
        cases.append((f"mutant T1 {name}", {"T1": {k: v[k] for k in ("p", "H", "C")}}, "fail", None))
    for name, v in mut["T2"].items():
        cases.append((f"mutant T2 {name}", {"T2": {k: v[k] for k in ("p", "R")}}, "fail", None))
    rows = []
    for name, res, expected, grader in cases:
        r = check(res)
        rows.append(dict(name=name, values=res, expected=expected, checker_result=r, grader_numeric_rows=grader))
    (here / "check_P015_selftest.json").write_text(json.dumps(rows, indent=1))
    for row in rows:
        print(f"{row['name']:<72s} expected={row['expected']:<10s} -> {json.dumps(row['checker_result'])}")


if __name__ == "__main__":
    here = Path(__file__).resolve().parent
    if len(sys.argv) > 1 and sys.argv[1] == "--selftest":
        selftest(here)
    else:
        print(json.dumps(check(json.loads(Path(sys.argv[1]).read_text())), indent=1))
