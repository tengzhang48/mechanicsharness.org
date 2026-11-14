#!/bin/bash
# Grade one packet with Grok 4.7 (Azure deployment grok-4.7) through opencode (written for opencode 2.2.6).
# usage: run_grader.sh PACKET_DIR        (a folder made by make_packet.py)
#
# What the grader can do is set by opencode.grok.json: every action is denied except read, glob and grep, and reading
# outside the packet is denied. So the grader has no tool to edit, to run a shell or code, to fetch or search the web
# or to start a sub-agent (tested 13 Nov 2026: see the README). Reasoning effort high (the deployment refuses "max").
# --standalone is needed: only a private server reads OPENCODE_CONFIG; the background service keeps its own settings.
# No --auto. extract_grade.py then accepts a grade only if every completed tool call stayed inside the packet.
#
# opencode sometimes ends a Grok run with "Duplicate tool input start" (its parallel tool calls), about one run in
# five in our tests: up to five attempts are made, and every attempt's stream is kept. The streams, the error output
# and the grade are written beside the packet (PACKET_DIR.out.N.jsonl, PACKET_DIR.err.N.txt, PACKET_DIR.times.txt,
# PACKET_DIR.grade/), never inside it. PACKET_DIR.grade/ holds grade.json only if an attempt was accepted.
# GRADER_EXTRA_PROMPT (environment, optional): one more sentence for the prompt. The bilayer benchmark uses it to ask
# for one tool call at a time, after five of five attempts on a folder of twenty files ended with that error.
set -u
HERE="$(cd "$(dirname "$0")" && pwd)"
# The form of the grade: this benchmark's form.json unless another benchmark gives its own. Since 14 Nov 2026 (Teng,
# "yes to all") it derives the summary word of question 2 from the listed checks (see extract_grade.py).
export GRADE_FORM="${GRADE_FORM:-$HERE/form.json}"
PACKET="$(cd "$1" && pwd)" || exit 2
cd "$PACKET" || exit 2
rm -rf "$PACKET.grade"
: > "$PACKET.times.txt"
for attempt in 1 2 3 4 5; do
  date -u +"attempt $attempt start %Y-%m-%dT%H:%M:%SZ" >> "$PACKET.times.txt"
  OPENCODE_CONFIG="$HERE/opencode.grok.json" timeout 2400 "$HOME/.opencode/bin/opencode" run --standalone \
    -m azure/grok-4.7 --format json \
    "Read GRADING.md in this folder and grade the submission exactly as it says. Read every file under submission/ before you answer.${GRADER_EXTRA_PROMPT:+ $GRADER_EXTRA_PROMPT}" \
    < /dev/null > "$PACKET.out.$attempt.jsonl" 2> "$PACKET.err.$attempt.txt"
  echo "attempt $attempt opencode exit $?" >> "$PACKET.times.txt"
  date -u +"attempt $attempt end %Y-%m-%dT%H:%M:%SZ" >> "$PACKET.times.txt"
  if python3 -I "$HERE/extract_grade.py" "$PACKET.out.$attempt.jsonl" "$PACKET" "$PACKET.grade"; then
    echo "attempt $attempt accepted" >> "$PACKET.times.txt"; exit 0
  fi
  echo "attempt $attempt not accepted" >> "$PACKET.times.txt"
done
exit 1
