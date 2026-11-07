#!/usr/bin/env python3
"""Proxy OpenCode with metadata-only stdout and discarded stderr.

Solver files are the scientific record. Private reasoning, tool arguments,
tool outputs and response text are not copied into coordinator logs.
"""
import json
import subprocess
import sys

OPENCODE = "<home>/.opencode/bin/opencode"


def public_metadata(event):
    if not isinstance(event, dict):
        return None
    part = event.get("part")
    part = part if isinstance(part, dict) else {}
    if event.get("type") in {"reasoning", "thinking"} or part.get("type") in {"reasoning", "thinking"}:
        return None
    result = {key: event[key] for key in ("type", "timestamp", "sessionID")
              if isinstance(event.get(key), (str, int, float))}
    details = {key: part[key] for key in ("id", "type", "sessionID", "tool", "reason")
               if isinstance(part.get(key), (str, int, float))}
    state = part.get("state")
    if isinstance(state, dict) and isinstance(state.get("status"), str):
        details["state"] = {"status": state["status"]}
    tokens = part.get("tokens")
    if isinstance(tokens, dict):
        details["tokens"] = {key: tokens[key] for key in ("input", "output", "reasoning", "total")
                             if isinstance(tokens.get(key), (int, float))}
        cache = tokens.get("cache")
        if isinstance(cache, dict):
            details["tokens"]["cache"] = {key: value for key, value in cache.items()
                                           if key in {"read", "write"} and isinstance(value, (int, float))}
    if details:
        result["part"] = details
    error = event.get("error")
    if isinstance(error, dict):
        safe_error = {key: error[key] for key in ("name", "type", "statusCode")
                      if isinstance(error.get(key), (str, int))}
        if isinstance(error.get("data"), dict):
            data = error["data"]
            safe_error["data"] = {key: data[key] for key in ("statusCode", "isRetryable")
                                   if isinstance(data.get(key), (int, bool))}
        result["error"] = safe_error
    return result or None


def main():
    process = subprocess.Popen([OPENCODE, *sys.argv[1:]], stdin=subprocess.DEVNULL,
                               stdout=subprocess.PIPE, stderr=subprocess.DEVNULL,
                               text=True, encoding="utf-8", errors="replace")
    for line in process.stdout:
        try:
            event = public_metadata(json.loads(line))
        except ValueError:
            continue
        if event is not None:
            print(json.dumps(event, sort_keys=True), flush=True)
    return process.wait()


if __name__ == "__main__":
    raise SystemExit(main())
