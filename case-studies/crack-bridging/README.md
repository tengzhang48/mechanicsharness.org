# Crack bridging: can AI agents couple a bridging law to a finite crack?

A MechanicsHarness case. The case page explains the mechanics, the solution and what each agent did:
https://mechanicsharness.org/case-studies/crack-bridging.html

The problems come from G. Bao and Z. Suo (1992), "Remarks on crack-bridging concepts", *Applied Mechanics Reviews* 45,
355–366: when a bridging law can be condensed into a toughness and when the load has to be computed from the body and
the law together. The agents received the problems without the paper. The paper is not included here.

This folder lets you attempt the two problems, check a new answer, and check every saved answer and grade against the
keys.

## The problems in brief

Probe round, 6 October 2026:

- **T1**, the paper's own finite-flaw problem at new inputs: a slit of half-length a = 2.40 mm in an infinite
  plane-strain plate under remote tension Σ, bridged as its fronts advance by connectors with a constant traction of
  160 MPa up to a full opening of 9.00 μm, plus an intrinsic front work of 22.0 J/m². Asked: Σ and the bridged length L
  per end when the first connector fails. Key: Σ = 133.73 MPa, L = 4.086 mm.
- **T2**, a variant: two unequal strips joined along a weak plane, opened by equal and opposite end couples, with a
  linear softening bridging law and a front work of 1.40 J/m². Asked: the couple M on each strip, the bridged length L
  and the opening at its midpoint when the first connector fails. Key: M = 0.2126 N·m, L = 10.80 mm, δ(L/2) = 38.22 μm
  (our computation; the paper solves only the equal-arm splitting beam). The audit found that δ(L/2) cannot detect a
  wrong arm stiffness (`FIXES.md`).

A third problem of the round, T3, is not here: it rested on invented observations (the case page, Section 10).

## What is here

| Path | Content |
|---|---|
| `task/T1`, `task/T2` | The problems exactly as given. Each run's `workspace/TASK.md` is the whole task it received, with the shared header; runs with the knowledge unit also have `workspace/knowledge/`. |
| `key/T1`, `key/T2` | The references, rubrics and reference scripts, unchanged since the runs. |
| `key/check_P038.py` | The numeric checker (added on 26 October 2026). |
| `FIXES.md` | What the audit of 25 October found and what was added. Read it before using the keys. |
| `audit_2026-10-25/` | That audit: its recomputation of the keys (`recompute_P038.py`, output `recompute_P038_results.json`), its checker run on the key, the graded answers and mutants, and its checks of T2's third output. Its passages on a restatement of T3 are left out (`RELEASE_EDITS.md`), because a restated T3 has not been run. |
| `regrade_2026-10-26/` | The checker run on the saved answers. |
| `runs/<run>/workspace/` | Every saved file of the eight runs, as delivered; `runs/<run>/session/`: their start, exit and error stamps. |
| `grading/` | The blind grade of each run, with the graders' scratch work, and the key builder's notes. |
| `records/` | Records cited by the case page: the probe pilot summary, the blind-label map and the harness notes, cut to this case (each says inside what was kept), the scores, the grader instructions, the run hashes, the setup and launch scripts, and the knowledge unit's freeze record. |
| `reference_checks/` | Our check of the reference values (`bridging_check.py`, output `bridging_check.txt`). |
| `check_package.py` | Checks this folder (below). |
| `RELEASE_EDITS.md` | Every change made to a file for this release. |
| `SHA256SUMS` | Hashes of every file. |
| `LICENSE`, `LICENSE-CONTENT` | MIT for the code; CC BY 4.0 for the text and data. |

Not included: session transcripts of the agents and graders, the source paper, and T3.

## Check an answer

Write `{"T1": {"Sigma_MPa": ..., "L_mm_per_end": ...}, "T2": {"M_Nm_each_strip": ..., "L_mm": ..., "delta_mid_um": ...}}`
(either part may be left out) to a file and run

    python3 key/check_P038.py results.json

It prints, for each number, its relative error, whether it is within the full-credit band (T1 0.3 %, T2 0.5 %) and the
rubric's numeric points (out of 50 for T1, 40 for T2). A value that is not a number is reported as such; an input
without a T1 or T2 block exits with status 1. Python 3, standard library only.

## Check this package

    python3 check_package.py

It copies this folder to a temporary directory, checks every file against `SHA256SUMS` (and that no other file is
present), reruns the checker's self-test (the key, the audit's independent recomputation, the eight graded answers and
the audit's mutants) and compares it with the saved one, and runs the checker on the key, on values just inside the
bands, on deliberately wrong answers (R03's numbers, Σ 1 % high, the moment summed over both strips), on a value given
as text and on inputs with no problem block. Its last line is `OK: all checks passed`. Python 3, standard library only.

## What needs human judgment

The checker compares numbers. Whether a run kept the front work separate from the connector work, coupled the bridged
zone to the finite geometry and verified its result independently needs a reader, and the rubrics' method and
verification rows are judgements. The blind grades in `grading/` are by Codex, as written at the time.

## The runs

| Agent | Model | Interface | Runs |
|---|---|---|---|
| GLM | GLM-5.3 | OpenCode | 4 (T1, T2; with and without the knowledge unit) |
| DeepSeek | DeepSeek-V4.1-Flash | Cline, Azure, reasoning off | 4 (the same cells) |

One run per cell, each with a 45-minute budget. Minutes, grades and numbers within tolerance for every run are in the
case page's Appendix A.

## Credit and status

Paper: G. Bao and Z. Suo. Answers in `runs/` were written by the AI agents named above. The problems, keys and rubrics
were prepared by Codex from the paper; the blind grades are by Codex; the audit and checker by AI assistants, as each
states. The case is maintained by Teng Zhang (Syracuse University); its checks and records were prepared with AI
assistants, as each record states. No expert has reviewed the keys or the grades. Corrections are welcome as issues in
the MechanicsHarness repositories.

Licence: the code (the `.py` files) is under the MIT License (`LICENSE`); the text and data are under CC BY 4.0
(`LICENSE-CONTENT`). The source paper is not included and is not covered.
