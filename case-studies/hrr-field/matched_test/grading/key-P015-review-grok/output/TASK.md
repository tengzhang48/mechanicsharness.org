# Task: review a mechanics knowledge topic, its checks and one knowledge entry

You are reviewing a prototype in `./library/crack_tip_asymptotics/`:
- `TOPIC.md`: what crack-tip and notch-tip asymptotic fields share (equations, assumptions,
  how the exponent is fixed, J and conventions, what changes outside the scope).
- `checks/crack_tip_checks.py`: constraint checks meant to accept any solver's angular stress
  field and report pass/fail with stated conditions and tolerances; `checks/test_checks.py`
  and `checks/mutant_detection.txt` show how they were validated.
- `entries/hutchinson_1968/`: the knowledge entry for Hutchinson J.W. (1968), J. Mech. Phys.
  Solids 16, 13–31 (`ENTRY.md`, `reference_data.json`, optional `reference_code/hrr.py`).
- `KNOWLEDGE_ENTRY_DESIGN.md`: the design these follow. The paper is here as `paper.txt` and
  page images `pages/page-NN.png` (300 dpi).

The purpose: an agent facing a *related new problem* should get the governing equations,
assumptions and constraints from the topic, be free to choose its own method, and use the
checks to test its own solution. The material must not act as a recipe that constrains it.

Report in `REVIEW.md`, **created within your first 10 minutes** with a checklist of the items
below, then updated as you go; it must be complete by minute 50. For each finding give
severity (high/medium/low), file:line, what is wrong, and your evidence (a derivation, a page of
the paper, or a script you ran in `./scratch/`).

1. **Mechanics correctness** of `TOPIC.md` and `ENTRY.md`: equations, assumptions and where
   they enter, the statements about how the exponent is fixed (crack vs notch), J and its
   conditions, conventions, the scope table. Check the entry against the paper.
2. **Try to break the checks.** (a) Construct fields that are wrong but pass (false negatives);
   (b) find correct fields that fail (false positives), e.g. other n, plane stress, other notch
   angles, coarse grids. Run your attempts with scripts in `./scratch/` and report each,
   including those the checks handled correctly. Are the stated conditions and tolerances right?
3. **Does anything constrain a solver** (recipe-like wording, hidden assumptions, missing
   parameters) or **mislead** one facing a variant (another geometry, mode, law)?
4. **What is missing** that a capable agent would need for a related problem?

End `REVIEW.md` with a verdict: usable as is / usable after the listed fixes / not usable.

Rules: work only inside this folder; do not edit `./library/` (copy to `./scratch/` to
experiment); no web search or fetch. Keep `PROCESS.md`: one line per step, each starting with a
time pasted from `date -u +%H:%M:%S` (never type a time). Time limit: 60 minutes.
Never stop or kill a process you did not start. Kill only by the PID you recorded when you
started it; never use pkill, killall, or pgrep-based loops (this machine is shared with other runs).
