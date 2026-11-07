"""Independent blinded review with the existing owned phase controller."""
import json
from pathlib import Path
import signal
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
import runner

runner.PROMPTS.update({
    "initial_analysis": (
        "Independent scientific review: initial_analysis. FIRST write STATUS.md "
        "with UTC time and incomplete obligations. Read TASK.md, PROBLEM.md "
        "and RUBRIC.md. Review the neutral submission and frozen reference "
        "under TASK.md. Derive deciding quantities independently; preserve "
        "your executed checking code and outputs. Save GRADE.md and grade.json "
        "incrementally, with partial judgments and incomplete checks explicit, "
        "plus the delivery aliases ANSWER.md and results.json. Do not edit "
        "submission/ or reference/. Provisional reporting starts at minute 8 "
        "and ends at minute 10; total review budget is 25 minutes including "
        "reporting. Finish and exit earlier when done. Work only here; no web, "
        "peer reviews, identity lookup, delegation or unowned process signals. "
        "Bound any calculation with a timeout or record its owned PID. "
        "No desired grade or coordinator interpretation is supplied."
    ),
    "checkpoint_report": (
        "Independent review checkpoint, minute 8 to 10. FIRST save "
        "PROVISIONAL_REPORT.md and update STATUS.md: completed judgments, "
        "actual independent checks, evidence locations, uncertainty and "
        "unfinished obligations. Use this session and saved work only. "
        "Preserve partial GRADE.md/grade.json and delivery aliases. Do not "
        "start new investigations or infer a desired verdict."
    ),
    "continued_analysis": (
        "Continue the independent review under TASK.md from your own saved "
        "work. FIRST update STATUS.md with UTC time and completed/incomplete "
        "obligations. Keep original submission/reference files unchanged. "
        "Save all seven rubric judgments, actual deciding checks, uncertainty "
        "and missing work in GRADE.md/grade.json, with ANSWER.md/results.json "
        "aliases for delivery. Final reporting begins at minute 20 and ends "
        "at minute 25. No scientific feedback or desired verdict is supplied."
    ),
    "final_report": (
        "Final independent-review reporting reserve, minute 20 to 25. FIRST "
        "update STATUS.md. Finalize GRADE.md and grade.json from existing "
        "derived evidence, with the seven rubric judgments, deciding checks, "
        "numerical agreement, uncertainty and incomplete obligations. Save "
        "identical ANSWER.md/results.json aliases and REPRODUCE.md. Do not "
        "start new investigations, repair the solver's answer, or invent "
        "checks or confidence."
    ),
})


class ReviewController(runner.Controller):
    def plan(self):
        result = super().plan()
        result.update(
            framework="OpenCode 1.17.18 through metadata-only proxy",
            role="Independent scientific grader",
            neutral_submission="submission/; identity and solver wall time withheld",
            coordinator_findings_supplied=False,
            desired_grades_supplied=False,
            privacy="Session/tool-status/usage metadata only; no private reasoning or stderr exported",
        )
        return result


if __name__ == "__main__":
    args = runner.parser().parse_args()
    controller = ReviewController(args)
    controller.validate()
    if args.dry_run:
        print(json.dumps(controller.plan(), indent=2, sort_keys=True))
    else:
        for signum in (signal.SIGTERM, signal.SIGINT):
            signal.signal(signum, lambda *_: setattr(controller, "stop_requested", True))
        result = controller.run()
        print(json.dumps({key: result[key] for key in ("outcome", "ended_at_utc", "elapsed_seconds", "error", "delivery_audit")}))
        raise SystemExit(1 if result["error"] else 0)
