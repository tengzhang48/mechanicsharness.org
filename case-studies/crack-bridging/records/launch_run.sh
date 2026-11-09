#!/bin/bash
# Launch one batch-1 run detached from the harness (lesson 39). Usage: launch_run.sh <name> <timeout_s>
set -euo pipefail
ROOT=<probe-runs>/runs; LOGS=<probe-runs>-logs
name=$1; tmo=$2; d=$ROOT/$name; agent=${name##*-}; mkdir -p $LOGS
P='Read TASK.md in this folder and complete it.'
case $agent in
  deepseek) cmd="opencode run --auto -m deepseek/deepseek-flash --dir $d --format json \"$P\"";;
  dsaz)     cmd="cline -P openai-compatible -m DeepSeek-V4.1-Flash --json -c $d -t $tmo \"$P\"";;
  k3az)     cmd="cline -P openai-compatible -m FW-Kimi-K3 --json -c $d -t $tmo \"$P\"";;
  gptaz)    cmd="opencode run --auto -m azure/gpt-6.1-sol --dir $d --format json \"$P\"";;
  grok)     cmd="opencode run --auto -m azure/grok-4.6 --dir $d --format json \"$P\"";;
  glm)      cmd="opencode run --auto -m zai-coding-plan/glm-5.3 --dir $d --format json \"$P\"";;
  kimi)     cmd="kimi -m kimi-code/kimi-for-coding -p \"$P\" --output-format stream-json";;
  codex)    cmd="codex exec -s danger-full-access --skip-git-repo-check --json -C $d -o $LOGS/$name.last.md -- \"$P\"";;
esac
cd "$d"
setsid nohup bash -c "date -u +%FT%TZ > $LOGS/$name.start; OPENCODE_EXPERIMENTAL_OUTPUT_TOKEN_MAX=131072 timeout $tmo $cmd > $LOGS/$name.jsonl 2> $LOGS/$name.err; echo \"exit \$? \$(date -u +%FT%TZ)\" > $LOGS/$name.exit" > /dev/null 2>&1 < /dev/null &
echo "launched $name (timeout $tmo s)"
