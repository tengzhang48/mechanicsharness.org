from pathlib import Path
import sys,json,signal
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
import runner
from subscription_runner import SubscriptionController

runner.PROMPTS.update({
    "initial_analysis": "Continuation phase: initial_analysis. Read TASK.md and REVIEW_GAPS.md. Your FIRST tool action must save STATUS.md and CONTINUATION_CHECKPOINT.md. Finish the final write-up from previous_attempt/ before expensive checks, save pending verification honestly, and continue only within this working folder. First handoff at minute 5; overall 18-minute budget includes reporting. No summary editing, peer/reference access or delegation.",
    "checkpoint_report": "Continuation phase: checkpoint_report. FIRST save PROVISIONAL_REPORT.md and update STATUS.md with current UTC time, results already saved, checks actually executed, and remaining gaps. This is the reserved minute-5-to-6 reporting phase. Use saved work only; do not start a new calculation. Save ANSWER.md/results.json/REPRODUCE.md even if explicitly partial.",
    "continued_analysis": "Continuation phase: continued_analysis. Continue TASK.md from your own saved files and REVIEW_GAPS.md. FIRST update STATUS.md. Audit the finite-difference check and physical scope honestly. Keep individual computations below 120 seconds and save outputs. Final reporting begins at minute 15 inside this new 18-minute continuation budget. Write only in this continuation folder.",
    "final_report": "Continuation phase: final_report. Reporting reserve from minute 15 to minute 18. FIRST update STATUS.md. Finish ANSWER.md, results.json, REPRODUCE.md and VERIFICATION.md from existing work, listing every incomplete obligation. Do not start new investigations or invent checks/grades. State this is a feedback-assisted continuation and the original attempt remains unchanged."
})
class ContinuationController(SubscriptionController):
    def plan(self):
        record=super().plan()
        record.update(scope="User-authorized feedback-assisted continuation; delivery audit only, no scientific grade",original_attempt="S04",new_scored_attempt=False,context="Previous public artifacts and reviewer gaps; prior private conversation not imported")
        return record

parser=runner.parser();parser.add_argument("--skills",type=Path,required=True)
args=parser.parse_args();controller=ContinuationController(args)
for signum in (signal.SIGTERM,signal.SIGINT):
    signal.signal(signum,lambda *_:setattr(controller,"stop_requested",True))
record=controller.run()
print(json.dumps({key:record[key] for key in ("outcome","ended_at_utc","elapsed_seconds","error","delivery_audit")}))
