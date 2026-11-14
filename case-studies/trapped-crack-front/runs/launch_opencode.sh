#!/bin/bash
# Launch one opencode run of the trapped-front benchmark in a Docker container that sees only its workspace (the task
# text) and a fresh opencode home (the model's thinking at its maximum, web tools denied, the key of one provider
# through the environment). No /home and no /media of this machine. The network stays on, because the model is reached
# through it; the run's log is audited afterwards for any use of the network other than the model.
# Teng, 12 Nov 2026: "call codex gpt 6.1-sol, kimi k3.1 and deepseek v4.1 to do the tests, like what we did before for
# the magnet case". 13 Nov 2026: "Apply, commit, push, launch"; "Web tools off, log audited"; "We probably do not need
# to test codex"; "we can see the other 3 first" (the three opencode agents of the magnet case: Kimi, DeepSeek, GLM).
# Image: the magnet case's (Python 3.12, NumPy 1.26, SciPy 1.11; nothing of any benchmark in it). Limit 3 hours, 2 CPUs.
# The launch is refused unless this benchmark's folder is committed and unchanged: the commit is the frozen task.
# Organizer only. usage: bash launch_opencode.sh NAME PROVIDER MODEL [VARIANT]      e.g. kimi_1 azure 'azure/Kimi-K3.1-Pro#max'
# VARIANT: none (default; the task text as it is) or reminder (the task text with the one-paragraph reminder of
# task/REMINDER.md, see ../make_workspace.py). The variant is recorded in the run folder (task_variant.txt).
set -e
NAME="${1:?name}"; PROVIDER="${2:?provider}"; MODEL="${3:?model}"; VARIANT="${4:-none}"
HERE=$(cd "$(dirname "$0")" && pwd)
BENCH=$(cd "$HERE/.." && pwd)
RUNS=<run-folders>
RUN=$RUNS/$NAME; WS=$RUN/workspace; OH=$RUN/ochome
IMAGE=magnet-bench-agent:20261111
OC_DIR=$(dirname "$(readlink -f "$(which opencode)")")
PROMPT="Read PUBLIC_TASK.md in this folder, then solve the problem it states. Work only inside this folder and do not use the web. End with a short summary of what you reached and what remains unverified."
[ -e "$RUN" ] && { echo "$RUN exists; refusing to overwrite"; exit 1; }
[ -z "$(git -C "$BENCH" status --porcelain -- .)" ] || { echo "uncommitted changes in $BENCH; commit before a run"; exit 1; }
docker image inspect "$IMAGE" > /dev/null
mkdir -p "$RUN"
python3 -I "$BENCH/make_workspace.py" "$WS" "$VARIANT"
echo "$VARIANT" > "$RUN/task_variant.txt"
mkdir -p "$WS/tmp"
( cd "$WS" && find . -type f -exec sha256sum {} \; | sort -k2 ) > "$RUN/workspace_hashes.txt"
git -C "$BENCH" rev-parse HEAD > "$RUN/benchmark_commit.txt"
python3 -I "$HERE/make_opencode_home.py" "$OH" "$PROVIDER"
opencode --version > "$RUN/opencode_version.txt"; echo "$MODEL" > "$RUN/model.txt"; echo "$PROMPT" > "$RUN/prompt.txt"
nohup timeout 3h docker run --rm --name "trapped_$NAME" --cpus 2 --user "$(id -u):$(id -g)" --env-file "$OH.env" \
  -v "$WS":/work -v "$OH":/ochome -v "$OC_DIR":/opt/opencode:ro \
  -e HOME=/ochome -e TMPDIR=/work/tmp -e LANG=C.UTF-8 -e TERM=dumb -w /work "$IMAGE" \
  /opt/opencode/opencode run --standalone --auto -m "$MODEL" --format json "$PROMPT" \
  < /dev/null > "$RUN/opencode.jsonl" 2> "$RUN/opencode.err" &
echo "$NAME pid $!"
sleep 20; rm -f "$OH.env"     # the container has the key in its environment by now
date -u +"launched $NAME ($MODEL) %Y-%m-%dT%H:%M:%SZ" | tee -a "$RUNS/launch_times.txt"
