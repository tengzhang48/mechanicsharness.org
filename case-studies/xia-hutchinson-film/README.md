# Xia–Hutchinson film: parallel channel cracks in a prestressed film

A MechanicsHarness case. The case page explains the mechanics, the solution and what each agent did:
https://mechanicsharness.org/case-studies/xia-hutchinson-film.html

The problem comes from Section 4.2 of Z.C. Xia and J.W. Hutchinson (2000), "Crack patterns in thin films", *J. Mech.
Phys. Solids* 48, 1107–1131 (the authors' copy: https://groups.seas.harvard.edu/hutchinson/papers/crackpatterns.pdf). The
agents received the problem without the paper. The paper is not included here. The statement did tell every agent that it
was a source-faithful reproduction task (`task/PROBLEM.md`, lines 3–4), so the agents knew that a published source
existed; `FIXES.md` explains how to read the runs with that cue.

This folder lets you attempt the task, check a new answer, and check every saved answer and grade against the key.

## The task in one paragraph

A prestressed elastic film (equibiaxial tension σ₀, plane stress, ν = 0.3) sits on a bed of springs and cracks through its
thickness along parallel lines. The task asks for the decay length ℓ of the crack wake, the energy release rate per
advancing tip for simultaneous cracks at spacing H (case A) and for new cracks midway between existing ones (case B), the
stress intensity, the critical stress and its limits, and sixteen normalized values: γ_A, γ_B and both K_I at
ℓ/H = 0.5, 1, 1.5 and 2. The decisive point is the crack face: it is free of total traction, σ₀ + Ē u₂′ = 0, so the
prestress enters the face condition and the energy. The key's closed forms are γ_A = tanh(H/2ℓ) and
γ_B = 2 tanh(H/2ℓ) − tanh(H/ℓ), with K_I/(σ₀√ℓ) = √((1 − ν²) γ).

## What is here

| Path | Content |
|---|---|
| `task/PROBLEM.md`, `task/TASK.md` | The problem and the task brief, exactly as given to every run. `task/D01/SETUP.json`: the set-up of the later DeepSeek run, which also received `runs/D01/workspace/RUN_PROTOCOL.md`, a one-page reporting schedule. |
| `key/reference_frozen/` | The key: the reference derivation (`REFERENCE.md`), its code (`reproduce.py`), the values (`reference.json`), one correction of a transcribed sentence (`CORRECTION_RECORD.md`), and Grok's review of the reference (`REVIEW_*`). |
| `key/RUBRIC.md`, `key/checks.py`, `key/KEY_FREEZE.json` | The grading rubric, the automatic checker (added 26 October 2026) and the record of the frozen key. |
| `FIXES.md` | What the audit of 25 October found, what was fixed, and Teng Zhang's decisions of 6 November 2026. Read it before using the key. `DECISION_PROBLEM_TEXT_2026-10-29.md`: the draft of the release annotation it quotes. |
| `audit_2026-10-25/` | The audit before the fixes, as written, with the statement edits it proposed (not applied). |
| `runs/<run>/workspace/` | Every saved file of the five released runs, as delivered: answers, results, code and outputs. `runs/D01/records/` and `runs/K01-continuation/records/`: their launch and completion records; K01's task is `CONTINUATION_TASK.md`. |
| `grading/R01`, `R03`, `R04`, `R05` | Grok 4.6's blind reviews of S02, S01, S04 and D01 (`output/GRADE.md`, `grade.json`; R05 in three phases). |
| `grading/<run>-replay/` | Unchanged reruns of each run's own reproduction command. |
| `grading/key-reference-review/` | Grok's review of the reference, as delivered (the reviewed copy is in `key/reference_frozen/`). |
| `regrade_2026-10-26/` | The checker run on every released answer: the expected checker outputs, and a summary. |
| `records/` | The study's records, as written at the time: graded profiles, the numerical comparison, model capabilities, the review mapping, the solver launches. |
| `records/Case_Studies_Review_2026-10-23.html`, `..._Response_...` | Excerpts of a review of six case pages and of its response: the parts about this case (marked inside). |
| `reference_checks/xia_hutchinson_check.py` | That review's own check of the sixteen values from the strip equation. |
| `check_package.py` | Checks this folder (below). |
| `RELEASE_EDITS.md` | Every change made to a file for this release. |
| `SHA256SUMS` | Hashes of every file. |
| `LICENSE`, `LICENSE-CONTENT` | MIT for the code; CC BY 4.0 for the text and data. |

Not included: session transcripts of the agents and reviewers, the source paper, and Kimi's two attempts through Azure
(S03), which stopped before any work was saved, with their review.

## Check an answer

Write `results.json` as an object with a list `rows`, one row per sampling location, each with `ell_over_H`, `gamma_A`,
`gamma_B`, `kI_A` and `kI_B` (numbers as numbers; `task/TASK.md` describes the format), then run

    python3 key/checks.py <folder containing results.json>

It prints a flag `<quantity>_<ell/H>_ok` for each of the sixteen values, the number within 0.5 % of the key, the number
within 10⁻⁶ absolute (the tolerance of the original record), and the values missing or not numbers. The statement gives
no accuracy and leaves the method open, so the checker passes a value within 0.5 %; every physically wrong answer the
audit tried fails at that band. Example: `python3 key/checks.py runs/S02/workspace` prints `regrade_2026-10-26/S02.json`
(7 of 16 within 0.5 %). It also reads `numeric_rows.json`, the file of run S04. Python 3, standard library only.

## Check this package

    python3 check_package.py

It copies this folder to a temporary directory, checks every file against `SHA256SUMS` (and that no other file is
present), recomputes the sixteen key values from the closed forms and compares them with `reference.json`, reruns the
checker on the five saved answers and compares each output with the saved one, and runs the checker on deliberately wrong
answers (the audit's: K_I without 1 − ν², the foundation energy left out, case B per period, the case B spacings swapped,
ℓ/H read as H/ℓ, a decay length 2ℓ), on values just outside and inside the band, on values that are not numbers and on
malformed files. Its last line is `OK: all checks passed`. Python 3, standard library only.

## Reproduce the key

The reference code and Grok's independent check need NumPy and SciPy and write output files next to themselves, so run
them in a copy outside this folder (otherwise `check_package.py` reports the new files):

    mkdir -p ../xh-check/reference_submission
    cp key/reference_frozen/reproduce.py key/reference_frozen/REVIEW_check_independent.py ../xh-check/
    cp key/reference_frozen/reference.json ../xh-check/reference_submission/
    python3 ../xh-check/reproduce.py
    python3 ../xh-check/REVIEW_check_independent.py

`reproduce.py` solves the strip problem numerically by two methods, integrates the stated energy, checks the closed forms,
runs three deliberately wrong controls, and writes `reference.json`: on 7 November 2026 its output was identical, byte for
byte, to `key/reference_frozen/reference.json`. Grok's check expects the reference in `reference_submission/`, the folder
it was written for; its output matched `REVIEW_check_independent_output.txt` apart from the lines naming the files it
wrote. The review's check, `python3 reference_checks/xia_hutchinson_check.py`, prints the sixteen values, the split between the
prestress work and the stored energy, and the other numbers of the review's section on this case, in about a second; no
output of it was saved at the time.

## What needs human judgment

The checker compares numbers, which is one of seven items of `key/RUBRIC.md`. A reviewer judges the rest: the wake
formulation and the face condition, the energy accounting of both geometries, the fracture modes, the thresholds and
their limits, the verification and its shared assumptions, and the physical scope of the foundation model. Close values
do not show a right argument: run S02's finite-element values are within 0.2–1.6 % of the key, while its answer states
that no relaxed state exists behind the tips. The statement's flaws (the source cue, no stated accuracy, a recognizable
source) are named in `FIXES.md`.

## The runs

| Run | Agent | Model | Interface | Reasoning setting | Condition | Start (UTC) | Minutes |
|---|---|---|---|---|---|---|---:|
| S01 | GPT | azure/gpt-6.1-sol | OpenCode | high (requested, listed in model profile) | task only | 2026-10-15 13:02 | 6.2 |
| S02 | DeepSeek | azure/DeepSeek-V4.1-Flash | OpenCode via Azure | provider default (requested high unsupported) | task only | 2026-10-15 13:02 | 50.0 |
| S04 | Kimi | kimi-code/k3 | Kimi Code CLI (subscription) | max thinking | task only | 2026-10-15 14:37 | 30.4 |
| D01 | DeepSeek | deepseek/deepseek-flash | OpenCode, direct API | thinking enabled; effort high | task + RUN_PROTOCOL.md reporting schedule | 2026-10-15 21:18 | 10.2 |
| K01-continuation | Kimi | kimi-code/k3 | Kimi Code CLI (subscription) | not recorded | continuation of S04 with the reviewer's comments; not part of the comparison | 2026-10-15 19:34 | 9.4 |

S01, S02, S04 and D01 had a fifty-minute budget, which included reporting (`task/TASK.md`; D01's `RUN_PROTOCOL.md`
changed only the reporting schedule); K01 had 18 minutes (`runs/K01-continuation/records/CONTINUATION_TASK.md`). S04 stopped at 30 minutes without a final
answer file; its saved files hold the derivation and all sixteen values (`numeric_rows.json`). D01 was a later run with
thinking explicitly enabled; provider, reasoning setting and reporting schedule all changed from S02, so their separate
effects are not measured. The reviews of S01, S02 and S04 were blind to the solver (`records/REVIEW_MAPPING.json`). The
run settings come from each run's launch records.

## Credit and status

Problem and paper: Z.C. Xia and J.W. Hutchinson. Answers in `runs/` were written by the AI agents named above. The problem
statement and the rubric were drafted by the coordinating GPT agent, from the same model family as run S01; the reference
by GLM-5.3, reviewed by Grok 4.6; the reviews by Grok 4.6. The case is maintained by Teng Zhang (Syracuse University); its
checker, checks and records were prepared with AI assistants, as each record states. No expert has reviewed the key or the
grades. Corrections are welcome as issues in the MechanicsHarness repositories.

Licence: the code (the `.py` files) is under the MIT License (`LICENSE`); the text and data are under CC BY 4.0
(`LICENSE-CONTENT`). The source paper is not included and is not covered.
