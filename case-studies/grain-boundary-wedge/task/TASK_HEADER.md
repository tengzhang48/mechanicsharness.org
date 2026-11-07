# Task

Solve the problem below. Budget: {BUDGET} minutes. Write `ANSWER.md` early and keep it updated, and
make sure it is complete before the budget ends. You may write and run code in this folder
(`python3` with numpy and scipy).

**Deliver three things:**

1. **`ANSWER.md`.** Your answers, the model and method you chose and why, your derivation in brief,
   and what is missing or uncertain. If anything in the problem is unclear or missing, list it under
   "Questions" with the assumption you used, and continue.
2. **`results.json`.** It holds every quantity listed under "Report" in the problem, with the names
   given there (numbers as numbers, formulas as strings), plus `reproduce_command`: the command of
   item 3.
3. **One command that reproduces your numbers,** named in ANSWER.md (for example `python3 solve.py`).
   It must print the values in `results.json`.

**For every check you make:** report what you compared, the output (pasted, not retyped), and how the
check could have failed.

Keep `PROCESS.md`: one line per step, each starting with a time pasted from `date -u +%H:%M:%S`
(never type a time).

Rules: use only this folder and your own knowledge. Do not open files outside this folder, and do
not use web search or web fetch. The session log will be checked. Never stop or kill a process you
did not start. Kill only by the PID you recorded when you started it; never use pkill, killall, or
pgrep-based loops (this machine is shared with other runs).
