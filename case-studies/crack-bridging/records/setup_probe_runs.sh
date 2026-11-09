#!/bin/bash
# Build probe-run folders: setup_probe_runs.sh <PID> <unit_dir_name> <solver> [<solver> ...]   (log P211)
# One folder per probe (T1..T3) x condition (with, without) x solver.
set -euo pipefail
B=<corpus-pilot-folder>; BUILD=<probe-runs>/$1-build-codex
RUNS=<probe-runs>/runs; pid=$1; unit=$B/knowledge_v1/$2; shift 2
[ -f "$unit/FREEZE.md" ] || { echo "unit not frozen: $unit"; exit 1; }
mkdir -p $RUNS
head_common='# Task

Solve the problem below. Write your complete answer to `ANSWER.md` in this folder: the requested
numbers or decisions, your method, and your checks. You may write and run code in this folder
(`python3` with numpy and scipy). Budget: 45 minutes; write `ANSWER.md` early, update it, and make
sure it is complete by 40 minutes. Keep `PROCESS.md`: one line per step, each starting with a time
pasted from `date -u +%H:%M:%S` (never type a time).
'
rules_with='Rules: use only this folder (including `./knowledge/`) and your own knowledge. Do not open files
outside this folder, and do not use web search or web fetch. The session log will be checked.
Never stop or kill a process you did not start. Kill only by the PID you recorded when you started it; never use pkill, killall, or pgrep-based loops (this machine is shared with other runs).

Provided material: `./knowledge/` contains a knowledge unit about a mechanics paper (an entry with a
"Using this unit" section, an extraction, a report, Python code and its tests). Use it as you see fit.
'
rules_without='Rules: use only this folder and your own knowledge. Do not open files outside this folder, and do
not use web search or web fetch. The session log will be checked.
Never stop or kill a process you did not start. Kill only by the PID you recorded when you started it; never use pkill, killall, or pgrep-based loops (this machine is shared with other runs).
'
for solver in "$@"; do for T in T1 T2 T3; do for cond in with without; do
  d=$RUNS/$pid-$T-$cond-$solver; [ ! -e "$d" ] || { echo "$d exists"; exit 1; }; mkdir -p "$d"
  if [ $cond = with ]; then rules=$rules_with; mkdir "$d/knowledge"
    (cd "$unit" && cp ENTRY.md EXTRACT.md REPORT.md *.py "$d/knowledge/")
  else rules=$rules_without; fi
  { printf '%s\n%s\n---\n\n' "$head_common" "$rules"; cat "$BUILD/$T/TASK_PROBLEM.md"; } > "$d/TASK.md"
done; done; done
ls $RUNS | grep "^$pid-" | wc -l
