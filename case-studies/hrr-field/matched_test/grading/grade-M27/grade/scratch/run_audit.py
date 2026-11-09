"""Run the supplied programs unchanged, retaining outputs and child PIDs."""
import json
from pathlib import Path
import subprocess
import sys
import time

root = Path(__file__).resolve().parents[1]
jobs = []
for label, program in (
    ("submission", root / "scratch" / "solve.py"),
    ("reference", root / "solve_reference.py"),
):
    output = (root / "scratch" / f"{label}_run.txt").open("w")
    started = time.monotonic()
    child = subprocess.Popen([sys.executable, "-u", str(program)], cwd=root,
                             stdout=output, stderr=subprocess.STDOUT)
    (root / "scratch" / f"{label}.pid").write_text(f"{child.pid}\n")
    print(f"Started {label}: PID {child.pid}", flush=True)
    jobs.append((label, child, output, started))

status = {}
for label, child, output, started in jobs:
    code = child.wait()
    output.close()
    status[label] = {"pid": child.pid, "exit_code": code,
                     "elapsed_seconds": time.monotonic() - started}
    print(f"Finished {label}: exit code {code}", flush=True)
(root / "scratch" / "run_status.json").write_text(json.dumps(status, indent=2) + "\n")
