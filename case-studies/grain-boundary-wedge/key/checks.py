"""Automatic checks for problem E (grading v2, log P344). Usage: python3 checks.py <run_folder>"""
import json, sys

L1, L2, SHAPE = 11.74265, 264.6348, 0.8020
try:
    r = json.load(open(f"{sys.argv[1]}/results.json"))
except Exception as e:
    print(json.dumps({"results_json": f"unreadable: {e}"})); sys.exit()
num = lambda k: r.get(k) if isinstance(r.get(k), (int, float)) else None
l1, l2, sh, kin, tr = num("lambda1"), num("lambda2"), num("mode1_shape_half"), num("K_inf_norm"), num("relaxation_time_t0")
out = {}
out["lambda1_class"] = None if l1 is None else ("exact" if abs(l1 - L1) / L1 < 0.01 else
                                                "approximate" if abs(l1 - L1) / L1 < 0.15 else "wrong")
out["lambda1_ok"] = None if l1 is None else abs(l1 - L1) / L1 < 0.01   # 26 Oct fix: counted with the other *_ok flags
out["lambda2_ok"] = None if l2 is None else abs(l2 - L2) / L2 < 0.02
out["mode1_shape_ok"] = None if sh is None else abs(sh - SHAPE) / SHAPE < 0.005   # 26 Oct fix: 2% admitted the Eq. B6 estimate 0.811
out["K_inf_ok"] = None if kin is None else 1.09 <= kin <= 1.15
out["relaxation_time_ok"] = None if tr is None else 0.05 <= tr <= 0.15
# principle checks on the agent's own numbers
out["own_lambda2_gt_lambda1"] = None if None in (l1, l2) else l2 > l1
out["own_time_consistent_with_lambda1"] = None if None in (l1, tr) or l1 <= 0 else 0.5 <= tr * l1 <= 2.0
print(json.dumps(out, indent=1))
