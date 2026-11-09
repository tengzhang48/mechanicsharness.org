"""Minimal automatic checker for the numerical rows of P038 T1 and T2 (audit draft, 2026-10-25).

It implements only what the rubrics make decidable from numbers: the tolerance bands of
T1/RUBRIC.md lines 14-15 and T2/RUBRIC.md lines 26-36. Method, verification and scope rows,
and all of T3, need a human or blind grader.

Input JSON schema (results file):
  {"T1": {"Sigma_MPa": float, "L_mm_per_end": float},
   "T2": {"M_Nm_each_strip": float, "L_mm": float, "delta_mid_um": float}}
Usage: python -I check_P038.py results.json   -> prints a JSON verdict
       python -I check_P038.py --selftest      -> key, runs and mutants (writes check_P038_selftest.json)
"""
import json, os, sys

KEY = {"T1": {"Sigma_MPa": 133.730286053811, "L_mm_per_end": 4.08560913824553},
       "T2": {"M_Nm_each_strip": 0.212625787376721, "L_mm": 10.8001023708973,
              "delta_mid_um": 38.2188979260933}}
# T1/RUBRIC.md:14-15  (pass = 25 points band)
T1_BANDS = [(0.003, 25), (0.01, 20), (0.05, 12), (0.15, 5)]
# T2/RUBRIC.md:26-35  weights 15/15/10; fractions of weight
T2_WEIGHTS = {"M_Nm_each_strip": 15, "L_mm": 15, "delta_mid_um": 10}
T2_BANDS = [(0.005, 1.0), (0.02, 0.8), (0.05, 0.5), (0.15, 0.2)]


def band(rel, bands):
    for tol, pts in bands:
        if rel <= tol:
            return pts
    return 0


def check(res):
    out = {}
    if "T1" in res:
        r = {}
        for k in KEY["T1"]:
            v = res["T1"].get(k)
            if v is None:
                r[k] = {"value": None, "pass": False, "points": 0}
                continue
            rel = abs(v - KEY["T1"][k]) / KEY["T1"][k]
            r[k] = {"value": v, "rel_err": rel, "pass": rel <= 0.003, "points": band(rel, T1_BANDS)}
        r["numeric_points_of_50"] = sum(x["points"] for x in r.values())
        out["T1"] = r
    if "T2" in res:
        r = {}
        for k, w in T2_WEIGHTS.items():
            v = res["T2"].get(k)
            if v is None:
                r[k] = {"value": None, "pass": False, "points": 0}
                continue
            rel = abs(v - KEY["T2"][k]) / KEY["T2"][k]
            r[k] = {"value": v, "rel_err": rel, "pass": rel <= 0.005, "points": w * band(rel, T2_BANDS)}
        r["numeric_points_of_40"] = sum(x["points"] for x in r.values())
        out["T2"] = r
    return out


def selftest():
    here = os.path.dirname(os.path.abspath(__file__))
    rec = json.load(open(os.path.join(here, "recompute_P038_results.json")))
    cases = []
    # (a) key-built results
    cases.append(("key", KEY, "all pass"))
    # (a') our independent recomputation
    cases.append(("independent_recompute",
                  {"T1": {"Sigma_MPa": rec["T1"]["Sigma_MPa_abel_route"], "L_mm_per_end": rec["T1"]["L_mm_abel_route"]},
                   "T2": {"M_Nm_each_strip": rec["T2"]["M_Nm_shoot"], "L_mm": rec["T2"]["L_mm"],
                          "delta_mid_um": rec["T2"]["delta_mid_um"]}}, "all pass"))
    # (b) saved runs: headline (bold) numbers as written in ANSWER.md, with the file:line they come from;
    #     grader numeric points from grading/Rxx/GRADE.md
    runs = json.load(open(os.path.join(here, "runs_extracted.json")))
    for name, d in runs.items():
        cases.append(("run:" + name, d["values"], "grader numeric points %s" % d["grader_numeric_points"]))
    # (c) mutants
    for k, (S, L) in rec["T1_mutants"].items():
        cases.append(("T1 mutant:" + k, {"T1": {"Sigma_MPa": S, "L_mm_per_end": L}}, "fail"))
    for k, (M, L, dm) in rec["T2_mutants"].items():
        cases.append(("T2 mutant:" + k, {"T2": {"M_Nm_each_strip": M, "L_mm": L, "delta_mid_um": dm}}, "fail"))
    out = []
    for name, res, expected in cases:
        c = check(res)
        flags = {t: {k: v["pass"] for k, v in c[t].items() if isinstance(v, dict)} for t in c}
        pts = {t: c[t].get("numeric_points_of_50", c[t].get("numeric_points_of_40")) for t in c}
        rel = {t: {k: (round(v["rel_err"], 6) if v.get("rel_err") is not None else None)
                   for k, v in c[t].items() if isinstance(v, dict)} for t in c}
        out.append({"name": name, "expected": expected, "pass_flags": flags, "rel_err": rel, "points": pts})
    json.dump(out, open(os.path.join(here, "check_P038_selftest.json"), "w"), indent=1)
    for o in out:
        print(o["name"], "| expected:", o["expected"], "| pass:", o["pass_flags"], "| pts:", o["points"])


if __name__ == "__main__":
    if sys.argv[1] == "--selftest":
        selftest()
    else:
        print(json.dumps(check(json.load(open(sys.argv[1]))), indent=1))
