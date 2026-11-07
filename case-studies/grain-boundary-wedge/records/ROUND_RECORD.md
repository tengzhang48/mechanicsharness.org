# Multiphysics crack-tip round: record (logs P336–P353; interim, results pending)

> **Excerpt for the grain-boundary wedge case, made on 7 November 2026.** Kept, unchanged: the title, incidents 76-79 and the challenges raised on problem E, with their tables' headers. Left out: the other incidents, the interim results, the other problems' challenges and the other sections.

| # | What | Handling |
|---|---|---|
| 76 | **K3 rate limit (Opus's error, Teng: "You probably launch too many k3 at on time").** Six at once exceeded the Azure deployment's token rate; 11 of 12 K3 runs failed | Attempts kept; rerun one at a time, A, B, D first and E last (from 15:14) |
| 77, 79 | **Cline ends a run on a text-only compaction-summary turn.** DeepSeek on E: MP-07 and MP-20 at about 6.5 min; MP-07's fresh rerun again at 19 min | Resuming the session with `--id` does not work in JSON mode (P348). MP-20 rerun running. MP-07 on hold: count as a harness failure, or rerun with `--compaction basic` (labelled) — Teng's call |
| 78 | **GLM output cap.** MP-32 (E): one reasoning step of 130,716 tokens hit the 131,072 cap; ended with nothing written | Rerun, as in U7 |

| Problem, part | Challenge (short) | Runs | Status |
|---|---|---|---|
| E (c) | The relaxation criterion is unspecified (e-folding against 95%) | GLM MP-03 | **A statement ambiguity** to fix |
| E | A uniform initial traction conflicts with zero surface traction at t = 0 | K3 MP-02 | open |
| E | The stated kernel has an unphysical sign; a π normalization is ambiguous | DeepSeek MP-18 | open |
