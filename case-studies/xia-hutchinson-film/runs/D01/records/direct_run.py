"""A fresh DeepSeek task attempt using the existing owned phase controller."""
import json
from pathlib import Path
import signal
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
import runner

runner.PROMPTS.update({
    "initial_analysis": (
        "Phase: initial_analysis. Your FIRST tool action must write STATUS.md "
        "with UTC time and incomplete work. Read TASK.md, PROBLEM.md and "
        "RUN_PROTOCOL.md. Solve the task from those authorized files only. "
        "Save derivations, ANSWER.md, results.json and executed checking outputs "
        "incrementally, including clearly marked partial work. The reporting "
        "schedule in RUN_PROTOCOL.md replaces only the older timing paragraph "
        "in TASK.md: provisional reporting starts at minute 5 and ends at "
        "minute 6; final reporting runs from minute 45 to minute 50. Finish "
        "and exit earlier if done. Do not delegate, fetch outside sources, "
        "inspect other runs/references or use pkill/killall. Work only in "
        "this folder and bound/record computations as TASK.md requires."
    ),
    "checkpoint_report": (
        "Phase: checkpoint_report, the minute-5-to-6 reporting reserve. FIRST "
        "write PROVISIONAL_REPORT.md and update STATUS.md: quantities already "
        "derived, actual checks and saved outputs, incomplete work and "
        "uncertainty. Save ANSWER.md and results.json even if explicitly "
        "partial. Use this same session and existing files; do not start a "
        "new calculation. No mathematical hint or desired result is supplied."
    ),
})


class DirectController(runner.Controller):
    def plan(self):
        result = super().plan()
        result.update(
            framework="OpenCode 1.17.18 through metadata-only proxy",
            provider="Direct DeepSeek API",
            thinking_configuration={"type": "enabled", "reasoningEffort": "high"},
            privacy="Only session/tool-status/usage metadata persisted; stderr and private content discarded",
            original_attempt="S02",
            comparison_limits="Fresh follow-up; provider, explicit reasoning and reporting schedule differ from S02. Not a controlled thinking off/on pair.",
            reference_or_peer_artifacts_provided=False,
        )
        return result


if __name__ == "__main__":
    args = runner.parser().parse_args()
    controller = DirectController(args)
    controller.validate()
    if args.dry_run:
        print(json.dumps(controller.plan(), indent=2, sort_keys=True))
    else:
        for signum in (signal.SIGTERM, signal.SIGINT):
            signal.signal(signum, lambda *_: setattr(controller, "stop_requested", True))
        record = controller.run()
        print(json.dumps({key: record[key] for key in ("outcome", "ended_at_utc", "elapsed_seconds", "error", "delivery_audit")}))
        raise SystemExit(1 if record["error"] else 0)
