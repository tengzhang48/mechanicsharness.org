"""Automatic checks for the Xia-Hutchinson film problem (added 26 Oct 2026, after the audit of 25 Oct).
Usage: python3 checks.py <run_folder>

Reads the rows of <run_folder>/results.json (or numeric_rows.json) and compares the sixteen requested values
(gamma_A, gamma_B, kI_A, kI_B at ell/H = 0.5, 1, 1.5, 2) with the frozen reference, reference_frozen/reference.json.
The statement gives no accuracy and leaves the method open, so a value passes within 0.5 % relative: every physically
wrong answer the audit tried (no 1 - nu^2, no foundation energy, Case B per period, swapped spacings, ell/H read as H/ell,
decay length 2 ell) fails at 0.5 %. The count at 1e-6 absolute (the record's frozen tolerance) is reported beside it."""
import json, os, sys

REF = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "reference_frozen", "reference.json")))
KEYS = ("gamma_A", "gamma_B", "kI_A", "kI_B")
run = sys.argv[1]
rows = None
for name in ("results.json", "numeric_rows.json"):
    try:
        d = json.load(open(os.path.join(run, name)))
    except Exception:
        continue
    rows = d if isinstance(d, list) else d.get("rows") if isinstance(d, dict) else None
    if isinstance(rows, list):
        break
if not isinstance(rows, list):
    print(json.dumps({"results_json": "no rows in results.json or numeric_rows.json"})); sys.exit()
num = lambda x: isinstance(x, (int, float)) and not isinstance(x, bool)
out, n_rel, n_abs, missing = {}, 0, 0, []
for ref in REF["rows"]:
    lh = ref["ell_over_H"]
    match = [r for r in rows if isinstance(r, dict) and num(r.get("ell_over_H")) and abs(r["ell_over_H"] - lh) < 1e-9]
    for k in KEYS:
        x = match[0].get(k) if len(match) == 1 else None
        if not num(x):
            out[f"{k}_{lh}_ok"] = False; missing.append(f"{k} at ell/H = {lh}"); continue
        out[f"{k}_{lh}_ok"] = abs(x - ref[k]) <= 5e-3 * abs(ref[k])
        n_rel += out[f"{k}_{lh}_ok"]; n_abs += abs(x - ref[k]) <= 1e-6
out["values_within_0.5pct"] = n_rel
out["values_within_1e-6_abs"] = n_abs
out["missing"] = missing
print(json.dumps(out, indent=1))
