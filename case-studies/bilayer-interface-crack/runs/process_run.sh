#!/bin/bash
# After one bilayer run has ended: audit its log, score it, apply the rubric, grade it with Grok, and copy its record
# into this folder (runs/NAME/). Organizer only.
# usage: bash process_run.sh NAME PACKET        e.g. process_run.sh kimi_1 B2   (PACKET: a neutral name for the grader)
# For a Codex run the copy of the login is deleted first. The mapping of packet names to runs is kept in
# RUNS/grading/MAPPING.txt and is never inside a packet.
set -u
NAME="${1:?name}"; PACKET="${2:?packet name}"
HERE="$(cd "$(dirname "$0")" && pwd)"; BENCH="$(cd "$HERE/.." && pwd)"
SHARED="$(cd "$BENCH/../trapped_front_mixed/runs" && pwd)"
RUNS=<run-folders>; RUN=$RUNS/$NAME
[ -d "$RUN/workspace" ] || { echo "no run $RUN"; exit 2; }
docker ps --format '{{.Names}}' | grep -qx "bilayer_$NAME" && { echo "bilayer_$NAME is still running"; exit 2; }
rm -f "$RUN/codexhome/auth.json"
if [ -f "$RUN/codex.log" ]; then python3 -I "$SHARED/audit_codex_log.py" "$RUN" > "$RUN/AUDIT.md"
else python3 -I "$SHARED/audit_log.py" "$RUN" > "$RUN/AUDIT.md"; fi
echo "audit exit $?"
VARIANT=$(cat "$RUN/task_variant.txt" 2>/dev/null || echo none)
python3 -I "$BENCH/scorer/score.py" "$RUN/workspace" "$RUN/score" "$VARIANT"; echo "score exit $? (variant $VARIANT)"
python3 -I "$BENCH/scorer/rubric.py" "$RUN/score"; echo "rubric exit $?"
mkdir -p "$RUNS/grading"; echo "$PACKET $NAME" >> "$RUNS/grading/MAPPING.txt"
bash "$BENCH/grader/grade.sh" "$RUN/workspace" "$RUN/score" "$RUNS/grading/$PACKET" "$VARIANT" > "$RUNS/grading/$PACKET.run.txt" 2>&1
echo "grader exit $? ($(grep -c start "$RUNS/grading/$PACKET.times.txt" 2>/dev/null) attempts)"
RECORD_DIR="$HERE" python3 -I "$SHARED/collect_record.py" "$RUNS" "$NAME" "$PACKET"; echo "record exit $?"
