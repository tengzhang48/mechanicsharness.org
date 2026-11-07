"""Retain local tool actions only, without private reasoning or provider state."""
from __future__ import annotations

import argparse
import datetime as dt
import json
from pathlib import Path
import subprocess
import time


REMOTE = r'''
import json
import sqlite3

connection = sqlite3.connect(
    "file:/home/solver/.local/share/opencode/opencode.db?mode=ro", uri=True
)
sql = """
SELECT m.id, m.seq,
       json_extract(p.value, '$.name'),
       json_extract(p.value, '$.state.status'),
       json_extract(p.value, '$.state.input'),
       json_extract(p.value, '$.state.content'),
       json_extract(p.value, '$.time')
FROM session_message m, json_each(m.data, '$.content') p
WHERE json_extract(p.value, '$.type') = 'tool'
ORDER BY m.seq, p.key
"""
actions = []
for mid, seq, name, status, inputs, content, times in connection.execute(sql):
    result = json.loads(content) if content else []
    action = {
        "message_id": mid, "message_sequence": seq,
        "tool": name, "status": status,
        "input": json.loads(inputs) if inputs else {},
        "time": json.loads(times) if times else None,
        "result_text_characters": sum(
            len(block.get("text", ""))
            for block in result if isinstance(block, dict)
        ),
    }
    if name == "shell":
        action["tool_result"] = result
    actions.append(action)
print(json.dumps(actions))
'''


def collect(container: str) -> list:
    result = subprocess.run(
        ["docker", "exec", "-i", container, "python", "-"],
        input=REMOTE, capture_output=True, text=True, timeout=20, check=True,
    )
    return json.loads(result.stdout)


def save(path: Path, record: dict) -> None:
    temporary = path.with_suffix(".tmp")
    temporary.write_text(json.dumps(record, indent=2) + "\n")
    temporary.replace(path)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("container")
    parser.add_argument("run_json", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--follow", action="store_true")
    args = parser.parse_args()
    failures = 0
    record = {}
    while True:
        try:
            actions = collect(args.container)
            record = {
                "observed_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
                "purpose": "Actual local tool actions for grading-process audit",
                "collection": (
                    "Read-only SQL selects tool inputs/status, result sizes and "
                    "shell tool results only. No assistant text, private reasoning, "
                    "provider state, credentials or complete database copied."
                ),
                "container": args.container,
                "collector_status": "observing" if args.follow else "snapshot",
                "actions": actions,
            }
            save(args.output, record)
            failures = 0
        except (subprocess.SubprocessError, json.JSONDecodeError) as error:
            failures += 1
            if record:
                record["collector_last_error_type"] = type(error).__name__
                save(args.output, record)
            if not args.follow:
                raise
        if not args.follow:
            print(json.dumps({"actions_retained": len(record["actions"])}))
            return
        run = json.loads(args.run_json.read_text())
        if run.get("status") != "running" or failures >= 3:
            if record:
                record["collector_status"] = "stopped; latest successful snapshot retained"
                record["collector_stopped_utc"] = dt.datetime.now(dt.timezone.utc).isoformat()
                record["snapshot_limit"] = (
                    "The container is ephemeral. Actions after the last successful "
                    "snapshot may be absent; retained event logs and output files "
                    "provide additional completion evidence."
                )
                save(args.output, record)
            return
        time.sleep(5)


if __name__ == "__main__":
    main()
