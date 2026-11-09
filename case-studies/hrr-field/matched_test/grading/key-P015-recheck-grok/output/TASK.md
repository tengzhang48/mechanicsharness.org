# Task: re-check a revised mechanics topic, its checks and one knowledge entry (version 0.2)

You (or another reviewer) reviewed version 0.1; that review is `PREVIOUS_REVIEW.md`. The revised
material is in `./library/crack_tip_asymptotics/` (`TOPIC.md` with a change log at the end,
`checks/`, `entries/hutchinson_1968/`). The paper is `paper.txt` with page images `pages/`.

Report in `RECHECK.md`, **created within your first 10 minutes** with a checklist, updated as
you go, complete by minute 35. For each item give status and evidence (a script in `./scratch/`,
a derivation, or a page of the paper).

1. **Each finding of `PREVIOUS_REVIEW.md`**: resolved / partly / not resolved / rejected with a
   correct reason. Rerun your earlier breaking attempts against the revised checks.
2. **New section 6 of `TOPIC.md`** ("Where the separable solution holds") and the new existence
   table in section 3: is each row's physics right? Is each basis label honest (derived /
   standard / from memory)? Is anything important missing or overstated?
3. **New problems introduced by the revision**, if any (try to break the new checks: crack
   opening, non-singular report, mode II Williams exponent, the `leading` switch, the
   field-scale normalization).

End with a verdict: freeze / freeze after the listed small fixes / not ready.

Rules: work only inside this folder; do not edit `./library/` (copy to `./scratch/`); no web
search or fetch. Keep `PROCESS.md`: one line per step, each starting with a time pasted from
`date -u +%H:%M:%S` (never type a time). Time limit: 40 minutes.
Never stop or kill a process you did not start. Kill only by the PID you recorded when you
started it; never use pkill, killall, or pgrep-based loops (this machine is shared with other runs).
