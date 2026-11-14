#!/bin/bash
# Launch one Codex run of the trapped-front benchmark in a Docker container that sees only its workspace (the task
# text) and a fresh Codex home. No /home and no /media of this machine. The network stays on, because the model is
# reached through it; the run's log is audited afterwards for any use of the network other than the model.
# Teng, 13 Nov 2026, after the first two results: "If it is good, we can call codex to do the benchmark."
# Same image, workspace, prompt, time limit and CPUs as launch_opencode.sh. Web tools off, as for the other agents:
#   web_search = "disabled"   Codex's own web search;
#   [features] apps = false   the account's connectors. They include a web search: in a smoke test with web_search
#                             disabled and the connectors on, a search through them worked. With both settings Codex
#                             reported no search or web tool (13 Nov 2026).
# Codex runs with --dangerously-bypass-approvals-and-sandbox INSIDE the container (its own sandbox cannot start here).
# A copy of the login (auth.json) is in the fresh Codex home during the run: never commit it, delete it after the run.
# The launch is refused unless this benchmark's folder is committed and unchanged.
# Organizer only. usage: bash launch_codex.sh NAME      (for example codex_1)
set -e
NAME="${1:?usage: launch_codex.sh NAME}"
HERE=$(cd "$(dirname "$0")" && pwd)
BENCH=$(cd "$HERE/.." && pwd)
RUNS=<run-folders>
RUN=$RUNS/$NAME; WS=$RUN/workspace; CH=$RUN/codexhome
IMAGE=magnet-bench-agent:20261111
CODEX_DIR=$(dirname "$(readlink -f <home>/.local/bin/codex)")
PROMPT="Read PUBLIC_TASK.md in this folder, then solve the problem it states. Work only inside this folder and do not use the web. End with a short summary of what you reached and what remains unverified."
[ -e "$RUN" ] && { echo "$RUN exists; refusing to overwrite"; exit 1; }
[ -z "$(git -C "$BENCH" status --porcelain -- .)" ] || { echo "uncommitted changes in $BENCH; commit before a run"; exit 1; }
docker image inspect "$IMAGE" > /dev/null
mkdir -p "$RUN" "$CH/home"
python3 -I "$BENCH/make_workspace.py" "$WS"
mkdir -p "$WS/tmp"
( cd "$WS" && find . -type f -exec sha256sum {} \; | sort -k2 ) > "$RUN/workspace_hashes.txt"
git -C "$BENCH" rev-parse HEAD > "$RUN/benchmark_commit.txt"
cp <home>/.codex/auth.json "$CH/" && chmod 600 "$CH/auth.json"
printf 'model = "gpt-6.1-sol"\nmodel_reasoning_effort = "max"\nmodel_reasoning_summary = "detailed"\nweb_search = "disabled"\n\n[features]\napps = false\n' > "$CH/config.toml"
<home>/.local/bin/codex --version > "$RUN/codex_version.txt"; echo "gpt-6.1-sol (reasoning effort max)" > "$RUN/model.txt"; echo "$PROMPT" > "$RUN/prompt.txt"
nohup timeout 3h docker run --rm --name "trapped_$NAME" --cpus 2 --user "$(id -u):$(id -g)" \
  -v "$WS":/work -v "$CH":/codexhome -v "$CODEX_DIR":/opt/codex:ro \
  -e HOME=/codexhome/home -e CODEX_HOME=/codexhome -e TMPDIR=/work/tmp -e LANG=C.UTF-8 -e TERM=dumb -w /work "$IMAGE" \
  /opt/codex/codex exec --strict-config --skip-git-repo-check --dangerously-bypass-approvals-and-sandbox -C /work -o /work/codex_final.md "$PROMPT" \
  < /dev/null > "$RUN/codex.log" 2>&1 &
echo "$NAME pid $!"
date -u +"launched $NAME (gpt-6.1-sol, effort max) %Y-%m-%dT%H:%M:%SZ" | tee -a "$RUNS/launch_times.txt"
