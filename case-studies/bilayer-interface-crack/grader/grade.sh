#!/bin/bash
# Grade one bilayer submission with Grok 4.7: build the packet, then run the grader. Organizer only.
# usage: bash grade.sh SUBMISSION_DIR SCORE_DIR PACKET_DIR [VARIANT]  (PACKET_DIR must not exist; give it a neutral name)
# VARIANT no_definitions: the instructions and form of variant_no_definitions/, and the task text the run received
# (SUBMISSION_DIR/PUBLIC_TASK.md) as TASK.md of the packet.
# The tools are the trapped front's (../../trapped_front_mixed/grader: make_packet.py, run_grader.sh, extract_grade.py,
# opencode.grok.json: every action denied except reading and searching inside the packet). This benchmark gives them
# its own instructions (GRADING.md, KEY_NOTES.md) and its own form (form.json). One difference from the trapped front:
# the summary word of question 2 is derived by the tool from the checks the grader lists (form.json, derive_summary),
# because on the trapped front Grok's own summary word disagreed with its list in two gradings and the grades were lost.
set -u
HERE="$(cd "$(dirname "$0")" && pwd)"
SHARED="$(cd "$HERE/../../trapped_front_mixed/grader" && pwd)"
VARIANT="${4:-none}"
if [ "$VARIANT" = none ]; then GD="$HERE"; TASK=""
else GD="$HERE/variant_$VARIANT"; TASK="$1/PUBLIC_TASK.md"; [ -f "$GD/GRADING.md" ] || { echo "no instructions in $GD"; exit 2; }; fi
GRADER_DIR="$GD" PACKET_BUNDLE=1 python3 -I "$SHARED/make_packet.py" "$1" "$2" "$3" $TASK || exit 2
# opencode often ends a run of this model with "Duplicate tool input start"; the more tool calls, the more often. On the
# first real folder (twenty files) ten of ten attempts ended with it, five of them with the request for one call at a
# time. So the packet also holds SUBMISSION_ALL.txt, every file of the submission in one, and the prompt points to it.
export GRADER_EXTRA_PROMPT="SUBMISSION_ALL.txt holds every file of submission/ with line numbers: read it instead of opening those files one by one, and cite the file names and line numbers it shows. Make one tool call at a time."
GRADE_FORM="$GD/form.json" bash "$SHARED/run_grader.sh" "$3"
