"""Record review completion without assigning an organizer scientific grade."""
from pathlib import Path
import datetime
import hashlib
import json
import math
import sys
import time

review = Path(__file__).resolve().parent
base = review.parents[1]
run_id = sys.argv[1]
run_folder = base / "runs" / run_id

def read(path):
    return json.loads(path.read_text())

def save(path, value):
    tmp = path.with_suffix(".watcher.tmp")
    tmp.write_text(json.dumps(value, indent=2) + "\n")
    tmp.replace(path)

while True:
    try:
        run = read(run_folder / "RUN.json")
    except (OSError, ValueError):
        time.sleep(5)
        continue
    if run.get("ended_utc"):
        break
    time.sleep(15)

workspace = run_folder / "workspace"
try:
    grade = read(workspace / "GRADING.json")
except (OSError, ValueError):
    grade = {}
score = grade.get("total_score")
score_is_number = isinstance(score, (int, float)) and not isinstance(score, bool) and math.isfinite(score)
oneshot = run_id == "R05_grok46_oneshot_azure_deepseek"
required = ["REVIEW.md", "GRADING.json"] if oneshot else ["REVIEW.md", "GRADING.json", "audit.py", "audit_results.json"]
artifacts = {name: (workspace / name).is_file() for name in required}
files_delivered = all(artifacts.values()) and score_is_number
status = "review files delivered; process audit pending" if files_delivered else "review incomplete; process exited"
freeze = read(review / ("ONESHOT_INPUT_FREEZE.json" if oneshot else "CONTINUATION_INPUT_FREEZE.json"))
completion_path = review / ("ONESHOT_COMPLETION_OBSERVED.json" if oneshot else "CONTINUATION_COMPLETION_OBSERVED.json")
input_mismatches = [name for name, digest in freeze["input_hashes"].items()
                    if not (workspace / name).is_file() or hashlib.sha256((workspace / name).read_bytes()).hexdigest() != digest]
result = {"checked_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
          "review_run": run_id, "status": status, "exit_code": run.get("exit_code"),
          "time_limit_reached": run.get("time_limit_reached", False), "ended_utc": run["ended_utc"],
          "artifacts": artifacts, "reported_score": score if score_is_number else None,
          "reported_verdict": grade.get("overall_verdict"), "input_mismatches": input_mismatches,
          "scientific_status": "Recorded external model output only; Codex process audit and consistency check still required."}
save(completion_path, result)
try:
    state = read(review / "STATE.json")
    if state.get("review_run") == run_id:
        state.update(status=status, ended_utc=run["ended_utc"], exit_code=run.get("exit_code"), reported_score=result["reported_score"], process_audit="pending")
        save(review / "STATE.json", state)
    path = base / "runs/P01_deepseek/EVALUATION_STATUS.json"
    evaluation = read(path)
    if evaluation.get("external_review", {}).get("run_id") == run_id:
        evaluation["scientific_review"] = status
        evaluation["external_review"].update(result)
        save(path, evaluation)
except (OSError, ValueError) as error:
    result["status_update_error"] = type(error).__name__
    save(completion_path, result)
print(json.dumps(result))
