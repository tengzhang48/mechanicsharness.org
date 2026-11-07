# Poroelastic crack: a suddenly loaded crack in a fluid-saturated solid

A MechanicsHarness case, level 1 (classical problems with known answers). The case page explains the mechanics, the
solution and what each agent did: https://mechanicsharness.org/case-studies/poroelastic-crack.html

The problem comes from E. Detournay and A. H.-D. Cheng (1991), "Plane strain analysis of a stationary hydraulic fracture
in a poroelastic medium", *Int. J. Solids Struct.* 27, 1645–1662, doi:10.1016/0020-7683(91)90067-P. The agents received
the problem without the paper. The paper is not included here.

This folder lets you attempt the task, check a new answer, and check every saved answer and grade against the key.

## The task in one paragraph

A crack of length 2L in an infinite Biot solid is loaded on its faces by a sudden unit pressure (loading 1) or a sudden
unit pore pressure (loading 2). The task asks for the stress intensity K and the crack volume V just after loading and at
long times, the leading growth of V at small times, and the pore pressure induced next to the faces, for ν = 0.2,
ν_u = 0.4 and α = 0.89. The decisive point is the tip: just after loading the crack as a whole responds undrained, but
the tip drains at once, so the key's K₁(0⁺) is 0.807593√(πL), not √(πL).

## What is here

| Path | Content |
|---|---|
| `task/PROBLEM.md` | The problem, exactly as given to every run. |
| `task/TASK_HEADER.md` | The instructions of the original round; each run received them with the budget set to 50 minutes, followed by the problem (`runs/MP-37/workspace/TASK.md`). The trials of 24 October received a shorter instruction (`runs/<trial>/records/prompt.txt`). |
| `key/KEY.md` | The key: the solution, its sources in the paper, the correction of 8 October and the grading choices. |
| `key/checks.py` | The automatic checker, as corrected on 26 October 2026; on 6 November 2026 it also stopped accepting `true`/`false` as numbers (no saved result changed). |
| `key/KEY_2026-10-08.md` | The key of 8 October, which the Grok reviewers received on 24 October; superseded by `KEY.md`. |
| `key/FE_CHECK.md`, `key/fe_check_*` | The finite-element check of the key (FEniCSx scripts, logs and outputs). |
| `FIXES.md` | What the audit of 25 October found, what was fixed in the key and the checker, and Teng Zhang's decisions of 6 November 2026. Read it before using the key. |
| `runs/<run>/workspace/` | Every saved answer of the eighteen graded runs, as delivered: `ANSWER.md`, `results.json` and code; for the original round also `PROCESS.md` and `TASK.md` (the prompt as received); for the trials `PROBLEM.md`. |
| `runs/<trial>/records/` | Prompt and run records of the three trials of 24 October. |
| `grading/MP-*/` | The blind grade of each original-round run (`grade/GRADE.md`; instructions in `grade/TASK.md`) and the automatic checks as first run (`auto/`). |
| `grading/grok-4.6-*/` | Grok 4.6's reviews of the trials of 24 October. |
| `grading/review_B_K10-*/` | The two factorization reviews (Codex with gpt-6.1-sol, and Grok) that led to the corrected key. |
| `regrade_2026-10-26/` | The corrected checker rerun on every saved answer: the expected checker outputs. |
| `independent_check_2026-10-29/` | An independent check of the tip term of the volume coefficient. |
| `audit_2026-10-25/` | The audit before the fixes, as written. |
| `reference_checks/` | Closed-form checks, the factorization rerun of 25 October, the checks of the case page's tip derivation (6 November) and the tip-factor curve of the case page. |
| `records/` | The records cited by the case page, as written at the time; four are excerpts that leave out other cases (marked inside). |
| `check_package.py` | Checks this folder (below). |
| `RELEASE_EDITS.md` | Every change made to a file for this release. |
| `SHA256SUMS` | Hashes of every file. |
| `LICENSE`, `LICENSE-CONTENT` | MIT for the code; CC BY 4.0 for the text and data. |

Not included: session transcripts of the agents and graders, the source paper, and three attempts that a provider rate
limit stopped before they were graded.

## Check an answer

Write the nine requested values to `results.json` with the names given under "Report" in `task/PROBLEM.md` (numbers as
numbers), then run

    python3 key/checks.py <folder containing results.json>

It prints JSON flags: `<name>_ok` for each value (true, false, or null when the value is missing or not a number;
`true` and `false` are not numbers), which
value K₁(0⁺) and the volume coefficient A match (`K1_0_matches`, `V1_short_matches`), and two checks of the answer
against itself. A value passes within 1% of the key, except K₁(0⁺) within 2%, A within 5% and K₂(0⁺) within 0.01
absolute. Example: `python3 key/checks.py runs/MP-37/workspace` prints the `checks` of `regrade_2026-10-26/MP-37.json`.
Python 3, standard library only.

## Check this package

    python3 check_package.py

It copies this folder to a temporary directory, checks every file against `SHA256SUMS` (and that no other file is
present), reruns the checker on all eighteen saved answers and compares each output with the saved one, runs the checker
on fourteen deliberately wrong answers (for example K₁(0⁺) = 1, the face term alone for A, −α for K₂(∞), `true` for a
number) to confirm that each is caught, on two changes inside the bands to confirm that nothing changes, and on two
malformed inputs. Its last line is `OK: all checks passed`.

## Reproduce the key's tip values

    python3 grading/review_B_K10-codex/output/scratch/wilson_factor.py --counts 16384 65536 262144 --scales 0.01

Codex's factorization prints K = 0.8075930 (K₁(0⁺) per √(πL)) and A = 2.4402196 on its last line, in about ten seconds
(NumPy). The same run on 25 October is in `reference_checks/wilson_rerun_2026-10-25.txt`. The tip term of A (1.086165)
was checked independently on 29 October to 1.09 ± 0.01 (`independent_check_2026-10-29/REPORT.md`).

The case page shows the tip calculation step by step. The numbers it quotes (the symbol's properties, n₁₁, n₂₁, K₁(0⁺),
the energy identity, the two terms of A and the weak-coupling slope) come from

    python3 reference_checks/tip_derivation_check.py grading/review_B_K10-codex/output/scratch/wilson_factor.py

which runs the same factorization, unchanged, in about 30 s; its output of 6 November 2026 is
`reference_checks/tip_derivation_check.txt`.

## What needs human judgment

The checker compares numbers. A grader decides the rest, following `grading/MP-37/grade/TASK.md` and the grading notes
in `key/KEY.md`:

- whether K₁(0⁺) was reasoned with a drained tip or taken as the crack-scale value;
- whether the small-time volume growth was reasoned (diffusion through the faces, the tip regions, reciprocity);
- the sign and the confinement in loading 2;
- every error, classified as decisive, side or presentation.

The statement says only "just after loading (t = 0⁺)", without the order of the tip and time limits. The key takes the
tip limit first (0.807593); the crack-scale value 1 is the other order and gets partial credit from a grader. The
checker marks it `outer_limit_1` and fails it. Two further flaws of the statement are named in `FIXES.md`: the face
shear stress is not prescribed, and the definition of `p_faces_0` conflicts with the face condition p = 0.

## The runs

| Run | Agent | Model | Interface | Reasoning setting | Start (UTC) | Minutes |
|---|---|---|---|---|---|---:|
| MP-06 | GLM | glm-5.3 | opencode run --auto (JSON) | not set in launch command | 2026-10-08 14:28 | 54.7 |
| MP-10 | GLM | glm-5.3 | opencode run --auto (JSON) | not set in launch command | 2026-10-08 14:28 | 28.2 |
| MP-13 | DeepSeek | DeepSeek-V4.1-Flash | Cline (JSON mode) | no --thinking flag (reasoning off by mistake) | 2026-10-08 14:34 | 16.5 |
| MP-21 | GLM | glm-5.3 | opencode run --auto (JSON) | not set in launch command | 2026-10-08 15:02 | 39.3 |
| MP-22 | Kimi | FW-Kimi-K3 | Cline (JSON mode) | not set in launch command | 2026-10-08 15:58 | 11.3 |
| MP-24 | Kimi | FW-Kimi-K3 | Cline (JSON mode) | not set in launch command | 2026-10-08 16:09 | 21.5 |
| MP-27 | DeepSeek | DeepSeek-V4.1-Flash | Cline (JSON mode) | no --thinking flag (reasoning off by mistake) | 2026-10-08 14:35 | 12.3 |
| MP-30 | DeepSeek | DeepSeek-V4.1-Flash | Cline (JSON mode) | no --thinking flag (reasoning off by mistake) | 2026-10-08 14:39 | 5.6 |
| MP-33 | Kimi | FW-Kimi-K3 | Cline (JSON mode) | not set in launch command | 2026-10-08 16:31 | 16.8 |
| MP-37 | GPT | gpt-6.1-sol | opencode run --auto (JSON) | not set in launch command | 2026-10-08 15:11 | 32.5 |
| MP-38 | GPT | gpt-6.1-sol | opencode run --auto (JSON) | not set in launch command | 2026-10-08 15:11 | 47.1 |
| MP-39 | GPT | gpt-6.1-sol | opencode run --auto (JSON) | not set in launch command | 2026-10-08 15:11 | 48.7 |
| MP-41 | DeepSeek | DeepSeek-V4.1-Flash | Cline (JSON mode) | --thinking high | 2026-10-08 20:54 | 15.6 |
| MP-46 | DeepSeek | DeepSeek-V4.1-Flash | Cline (JSON mode) | --thinking high | 2026-10-08 21:12 | 29.7 |
| MP-50 | DeepSeek | DeepSeek-V4.1-Flash | Cline (JSON mode) | --thinking high | 2026-10-08 21:32 | 14.6 |
| kimi-k31-high-2026-10-24 | Kimi | kimi-code/k3.1-pro | Kimi Code CLI 2.3.1 | high | 2026-10-24 12:18 | 7.9 |
| kimi-k31-max-2026-10-24 | Kimi | kimi-code/k3.1-pro | Kimi Code CLI 2.3.1 | max | 2026-10-24 12:52 | 28.8 |
| opus-5.5-subagent-2026-10-24 | Claude (Opus) | claude-opus-5-5 | Claude Code sub-agent | effort max | 2026-10-24 13:48 | 93.6 |

The original round (8 October) had a 50-minute budget; the three trials of 24 October had none. The blind grader of the
original round was Codex with gpt-6.1-sol, the same model family as the GPT runs. The run settings come from each run's
launch record; "not set in launch command" means the provider default.

## Credit and status

Problem and paper: E. Detournay and A. H.-D. Cheng. Answers in `runs/` were written by the AI agents named above, and
the reviews by the AI reviewers named in each file. The case is maintained by Teng Zhang (Syracuse University); its key,
checker, checks and records were prepared with AI assistants, as each record states. No expert has reviewed the key or
the grades. Corrections are welcome as issues in the MechanicsHarness repositories.

Licence: the code (the `.py` files) is under the MIT License (`LICENSE`); the text and data are under CC BY 4.0
(`LICENSE-CONTENT`). The source paper is not included and is not covered.
