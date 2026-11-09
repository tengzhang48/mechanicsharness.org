"""Numeric checker for the P015 matched-test problems V1 and V2 (audit 2026-11-07).

Applies the tolerances stated in the problems (V1/TASK_PROBLEM.md:34: all seven numbers 2% relative;
V2/TASK_PROBLEM.md:48: T 1% relative, the other six 2% relative) to signed values:
    pass  <=>  |value - key| <= tol * |key|.
Keys: reference_output.json route A of each problem (recomputed independently in this folder,
recompute_V1.py / recompute_V2.py, agreement <= 4e-11). This file holds key values: HIDDEN.

Input: JSON {"V1": {slip_0.02, slip_0.08, tau_rz, tau_tz, r_0, r_90, r_180},
             "V2": {T, C_rr, C_tt, C_g, C_y, r_min, r_max}}  (either block may be absent).
Usage:  python -I check_V1_V2.py results.json         -> per-item report, exit 0 iff all present pass
        python -I check_V1_V2.py --audit               -> key self-test, graded answers, mutants;
                                                          writes check_V1_V2_audit.json / .out
Not covered by any numeric check: the justification, verification and initiation/scope sentences
(rubric-only), the 'nonempty interval' statement of V2 item 4.
"""
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
KEY = {
    "V1": {"slip_0.02": 0.001318126777760562, "slip_0.08": 0.0017936937381386464,
           "tau_rz": 38.586596862313044, "tau_tz": 138.01061663235828,
           "r_0": 0.046894423054945324, "r_90": 0.025066123488862637, "r_180": 0.013398406587127234},
    "V2": {"T": 4816.472187302739, "C_rr": 0.003942048237733179, "C_tt": 0.0006331781350710297,
           "C_g": 0.005731131133456488, "C_y": 0.009139018941480233, "r_min": 0.009999999999999998,
           "r_max": 0.06132522854990914},
}
TOL = {"V1": {k: 0.02 for k in KEY["V1"]},
       "V2": {**{k: 0.02 for k in KEY["V2"]}, "T": 0.01}}


def check(res):
    """Return {problem: {item: (value, rel_err, pass)}} for the blocks present."""
    out = {}
    for prob in ("V1", "V2"):
        if prob not in res:
            continue
        out[prob] = {}
        for k, kv in KEY[prob].items():
            v = res[prob].get(k)
            if v is None:
                out[prob][k] = (None, None, False)
                continue
            rel = abs(v - kv) / abs(kv)
            out[prob][k] = (v, rel, rel <= TOL[prob][k] * (1 + 1e-12))
    return out


def summary(chk, prob):
    items = chk[prob]
    npass = sum(1 for v in items.values() if v[2])
    fails = [k for k, v in items.items() if not v[2]]
    return npass, len(items), fails


def audit():
    lines = []
    P = lambda s="": lines.append(s)
    rec = {"key_selftest": {}, "answers": {}, "mutants": {}}
    chk = check(KEY)
    for prob in ("V1", "V2"):
        n, m, f = summary(chk, prob)
        rec["key_selftest"][prob] = f"{n}/{m}"
        P(f"key self-test {prob}: {n}/{m} pass")
    ans = json.load(open(os.path.join(HERE, "answers_extracted.json")))
    P("\ngraded answers (numbers from ANSWER.md; grader verdict from GRADE.md)")
    for prob in ("V1", "V2"):
        for lab, a in ans[prob].items():
            c = check({prob: a["values"]})
            n, m, f = summary(c, prob)
            worst = max(v[1] for v in c[prob].values())
            rec["answers"][lab] = {"run": a["run"], "checker": f"{n}/{m}", "fails": f,
                                   "max_rel_err": worst, "grader": a["grade_verdict"],
                                   "grade_total": a["grade_total"]}
            P(f"  {prob} {lab:9s} {a['run']:24s} checker {n}/{m}  grader {a['grade_verdict']:22s}"
              f" max rel err {worst:.3e}  fails {f}")
    mut = json.load(open(os.path.join(HERE, "mutants.json")))
    P("\nmutants (expected: fail)")
    for name, mm in mut.items():
        prob = mm["problem"]
        vals = dict(KEY[prob]); vals.update(mm["values"])      # untouched items keep key values
        c = check({prob: vals})
        n, m, f = summary(c, prob)
        caught = len(f) > 0
        rec["mutants"][name] = {"problem": prob, "changed": mm["values"], "checker": f"{n}/{m}",
                                "fails": f, "caught": caught,
                                "passes_changed": [k for k in mm["values"] if c[prob][k][2]]}
        P(f"  {prob} {name:44s} checker {n}/{m}  {'CAUGHT' if caught else 'PASSES (not caught)'}"
          f"  fails {f}  changed-but-passing {rec['mutants'][name]['passes_changed']}")
    json.dump(rec, open(os.path.join(HERE, "check_V1_V2_audit.json"), "w"), indent=1)
    open(os.path.join(HERE, "check_V1_V2_audit.out"), "w").write("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "--audit":
        audit()
        sys.exit(0)
    res = json.load(open(sys.argv[1]))
    chk = check(res)
    ok = True
    for prob, items in chk.items():
        for k, (v, rel, p) in items.items():
            ok &= p
            print(f"{prob} {k:10s} value {v!s:>14s}  key {KEY[prob][k]:.10g}  rel err "
                  f"{('%.3e' % rel) if rel is not None else 'missing':>10s}  tol {TOL[prob][k]:.0%}  "
                  f"{'PASS' if p else 'FAIL'}")
    sys.exit(0 if ok else 1)
