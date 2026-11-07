# Grading instructions (grading v2)

You grade one anonymous submission to a mechanics problem. You report a **profile**, not only a score.
The folder contains:
- `TASK_PROBLEM.md`: what the solver was asked;
- `KEY.md`: the key, with its sources and its own notes on grading;
- `AUTO_CHECKS.json`: automatic results, produced before you;
  - the solver's reproduction command was rerun on a copy, and its printed output was compared with
    `results.json` (`claims`);
  - problem-specific checks compared the reported values with the key, and with the stated
    equations (`checks`);
- `submission/`: the solver's files (ANSWER.md, results.json, code, any subfolders).

Work only inside this folder. You may rerun the solver's code on a copy in `./scratch/`. No web
search or fetch. Time limit: 25 minutes. Keep `PROCESS.md` (times pasted from `date -u +%H:%M:%S`).

**Judge and report each of these.**

1. **Answers.** For each lettered question: right, partly right or wrong, against the key and within
   its stated bands, with one line of reason. Credit what the answer and its code actually show,
   not what they claim.
2. **Supporting state.** Are the intermediate quantities the problem asked for right?
3. **Method.** What method did the solver use? Is it the right method for the needed solution (as
   the key describes), a valid alternative, or a wrong or unnecessarily general route?
4. **Errors.** List **every** error, each with a short quote or file:line, classified as:
   - **decisive:** it changes an answer;
   - **side:** wrong but not changing an answer;
   - **presentation:** wording or units.

   Include errors in reasoning that happen to give a right number.
5. **Checks.** For each check the solver reported: did it test the deciding quantity, and could it
   have failed if the solution were wrong? Count those that could fail.
6. **Claims.** Do the delivered code and its output support the numbers and checks claimed? Use
   AUTO_CHECKS.json, and rerun if needed. Report mismatches.
7. **Questions.** Did the solver list questions or missing information? Were they sensible?
8. **Source.** Did the solver name a source paper or author for this problem? Quote it.
9. **Challenges.** Did the solver argue that the problem, the model behind the key, or a published
   result is wrong, incomplete or ambiguous? For each challenge give the claim, the argument, and
   whether it is checkable (a computation or limit someone else could run).
   - A well-argued, checkable challenge is **not** counted as an error. Grade the affected answer
     both ways (against the key; against the challenge) and mark it "contested".
   - Challenging is the first step to new insight (Teng); verification decides it later.

If you believe the key itself is wrong somewhere, say so with evidence, and grade both ways.

**Write `GRADE.md`.** Give a short section for each item above, then end the file with exactly this
JSON block, and nothing after it:

```json
{"answers": {"a": "right|partly|wrong|contested", "b": "...", "c": "..."},
 "supporting_state": "right|partly|wrong",
 "method": "<short name>", "method_class": "right|valid_alternative|general_or_wrong",
 "errors": {"decisive": 0, "side": 0, "presentation": 0},
 "checks_reported": 0, "checks_could_fail": 0,
 "claims_supported": "yes|partly|no",
 "questions_listed": 0,
 "source_named": "<quote or none>",
 "challenges": [{"claim": "<short>", "checkable": true}]}
```

Use the letters of the problem's questions as the keys of "answers".

Never stop or kill a process you did not start. Kill only by the PID you recorded when you started
it; never use pkill, killall, or pgrep-based loops (this machine is shared with other runs).
