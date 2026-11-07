"""Data for the round's HTML summary (log P390): one record per run from the label map, the run
logs (start, exit), the automatic checks and the blind grades (profiles.json). Writes summary_data.json.
Usage: python3 grade_v2.py profile; python3 summary_data.py"""
import json, os, datetime as dt
H = os.path.dirname(os.path.abspath(__file__))
B = "<round-runs>"
LAB = json.load(open(f"{H}/LABELS_runs_coordinator_only.json"))
prof = {r["label"]: r for r in json.load(open(f"{H}/profiles.json"))}
iso = lambda s: dt.datetime.fromisoformat(s.strip().replace("Z", "+00:00"))
out = []
for l, v in sorted(LAB.items()):
    lg = f"{B}/logs/{l}"
    ex = open(lg + ".exit").read().split() if os.path.exists(lg + ".exit") else []
    mins = None
    if ex and os.path.exists(lg + ".start"):
        try:
            mins = round((iso(ex[-1]) - iso(open(lg + ".start").read())).total_seconds() / 60)
        except Exception:
            pass
    p = prof.get(l, {})
    a = p.get("auto_checks") or {}
    keyok = [x for k, x in a.items() if k.endswith("_ok") and x is not None]
    g = p.get("grade") if isinstance(p.get("grade"), dict) else None
    out.append(dict(label=l, agent=v["agent"], problem=v["problem"], repeat=v.get("repeat"), minutes=mins,
                    exit=" ".join(ex[:2]) if ex else None, answered=bool(p.get("answer")),
                    key_checks=[sum(1 for x in keyok if x is True), len(keyok)] if keyok else None,
                    claims=p.get("claims"),
                    parts=g["answers"] if g else None, method=g["method_class"] if g else None,
                    errors=g["errors"] if g else None,
                    checks=[g["checks_could_fail"], g["checks_reported"]] if g else None,
                    source=g.get("source_named") if g else None,
                    challenges=[c.get("claim") for c in (g.get("challenges") or [])] if g else None))
json.dump(out, open(f"{H}/summary_data.json", "w"), indent=1)
print(len(out), "runs;", sum(r["answered"] for r in out), "answered;", sum(1 for r in out if r["parts"]), "graded")
