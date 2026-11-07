# Grain-boundary diffusion wedge: how a stressed film relaxes into a crack

A MechanicsHarness case. The case page explains the mechanics, the solution and what each agent did:
https://mechanicsharness.org/case-studies/grain-boundary-wedge.html

The problem comes from H. Gao, L. Zhang, W. D. Nix, C. V. Thompson and E. Arzt (1999), "Crack-like grain-boundary
diffusion wedges in thin metal films", *Acta Mater.* 47, 2865–2878, doi:10.1016/S1359-6454(99)00178-0 (the single-wedge
problem of Sections 2.1 and 3.1). The agents received the problem without the paper. The paper is not included here.
It was problem E of a round of four coupled-field crack problems run on 8 October 2026.

This folder lets you attempt the task, check a new answer, and check every saved answer and grade against the key.

## The task in one paragraph

A polycrystalline film of thickness h lies on a substrate with the same elastic constants; a grain boundary runs from
the free surface to the interface and starts under a uniform tension σ₀. Atoms diffuse from the surface into the
boundary and plate onto it as a wedge, which lowers the traction. The traction is zero at the surface, and at the
interface there is no flux and no opening rate. The task asks (a) how the boundary traction σ_gb(z, t) evolves, with the
late-time decay rate λ₁ and the next rate λ₂; (b) what the boundary becomes as t → ∞, and the stress intensity factor at
the interface in units of σ₀√(πh); and (c) the time scale of relaxation, in units of t₀. In the key, λ₁ = 11.74265 and
λ₂ = 264.6348 (the paper's Table 1; recomputed as 11.742662 and 264.635), the boundary becomes a traction-free edge crack
with K∞ = 1.1215 σ₀√(πh) (the model's own limit, 1.121522), and the time scale is 1/λ₁ ≈ 0.085 t₀, graded in a band of
0.05–0.15 t₀ because the statement defines no measure (`FIXES.md`).

## What is here

| Path | Content |
|---|---|
| `task/PROBLEM.md`, `task/TASK_HEADER.md` | The problem and the round's task header (the answer contract), exactly as given to every run. |
| `key/KEY.md`, `key/checks.py` | The key and the automatic checker as they stand, fixed on 26 October 2026 (`FIXES.md`). |
| `key/KEY_2026-10-08.md`, `key/checks_2026-10-08.py` | The key and checker of 8 October, with which the runs were graded and counted; superseded (a notice at the top says so). |
| `key/reproduce_keys.py` | The round's reproduction of the key's own numbers, the part for this problem. |
| `FIXES.md` | What the audit of 25 October found and what was fixed. Read it before using the key. |
| `audit_2026-10-25/` | The audit, its proposed statement and checker edits (not applied), and its recomputation of the key (`key_recompute.py`, output `key_recompute.out`, collocation up to N = 480). |
| `runs/<run>/workspace/` | Every saved file of the twelve graded runs, as delivered. |
| `grading/<run>/grade/` | The blind grade of each answered run (`GRADE.md`), with the grader's scratch work; `grading/<run>/auto/`: the automatic checks at the time. |
| `regrade_2026-10-26/` | The fixed checker run on every saved answer, and a summary. |
| `records/` | Excerpts of the round's records cited by the case page (each says inside what was kept): the round summary, round notes, harness notes, the milestone summary, the process log and the round's run data (`summary_data.json`, problem E only), and the two scripts behind the summary, whole. `MP-18_web_requests_excerpt.txt`: the web requests of run MP-18, which broke the round's rule against web search; the content it fetched is left out. |
| `reference_checks/` | A check of the reference values cited by the case page (`grain-boundary_check.py`, output `grain-boundary_check.txt`). Its sections 1–5 (eigenvalues, traction-free limit, Appendix B, decay ratios, late-time amplitude) rerun from this folder and print the saved numbers; sections 6 and 7 read the runs in the case library's folders, so they stop here, and their output is kept as written. |
| `check_package.py` | Checks this folder (below). |
| `RELEASE_EDITS.md` | Every change made to a file for this release. |
| `SHA256SUMS` | Hashes of every file. |
| `LICENSE`, `LICENSE-CONTENT` | MIT for the code; CC BY 4.0 for the text and data. |

Not included: session transcripts of the agents and graders, the source paper, and the attempts that the
infrastructure stopped and that were rerun (rate limits, a compaction turn, an output cap, a text-only turn; the round
notes record each).

## Check an answer

Write `results.json` with `lambda1`, `lambda2`, `mode1_shape_half` (the slowest mode at z/h = 0.5 with f(1) = 1),
`K_inf_norm`, `relaxation_time_t0` and the method, then run

    python3 key/checks.py <folder containing results.json>

It prints the class of λ₁ (exact within 1 %, approximate within 15 %, otherwise wrong), pass or fail for λ₁ (1 %), λ₂
(2 %), the mode shape (0.5 %), K∞ (1.09 to 1.15) and the relaxation time (0.05 to 0.15 t₀), and two checks of the
answer's own numbers (λ₂ > λ₁; the relaxation time consistent with λ₁). A value that is not a number is not counted.
Example: `python3 key/checks.py runs/MP-40/workspace` prints `regrade_2026-10-26/MP-40.json`'s flags (λ₁ approximate).
Python 3, standard library only. Passing these flags does not show the eigenproblem was solved: K∞ = 1.1215 can be
recalled, and several runs typed it in (the case page, Appendix A).

## Check this package

    python3 check_package.py

It copies this folder to a temporary directory, checks every file against `SHA256SUMS` (and that no other file is
present), reruns the checker on the eleven saved answers it was run on and compares the flags with the saved ones,
and runs it on the key's own values and on deliberately wrong answers: the paper's Appendix B estimates (λ₁ = 13.36,
shape 0.811), K = 1, λ₁ just outside and inside its band, λ₂ below λ₁, a relaxation time out of its band, a value given
as text, and a missing results file. Its last line is `OK: all checks passed`. Python 3, standard library only.

## Reproduce the key

`python3 key/reproduce_keys.py` (NumPy, SciPy) solves the eigenproblem by the paper's Gauss–Chebyshev collocation and
prints, at N = 120, λ₁ = 11.742662, λ₂ = 264.6346, f₁(0.5)/f₁(1) = 0.8020, the traction-free limit 1.1197 (the audit
extrapolates 1.121522) and the crack-shape estimate of Eq. B8, 13.3595, in a few seconds (run on 7 November 2026). The audit's `key_recompute.py` goes to
N = 480 with Richardson extrapolation.

## What needs human judgment

The checker compares numbers. Whether a run solved the eigenproblem, kept the boundary conditions and computed K∞ from
the traction-free equation, rather than recalling it, needs a reader; so does the relaxation time, which the statement
does not define. The blind grades in `grading/` are by Codex, as written at the time; MP-47 was not graded. With the fixed checker, λ₁ is now one of the counted flags; the round's summary counted four (`FIXES.md`).

## The runs

| Run | Agent | Model | Interface | Reasoning | Minutes |
|---|---|---|---|---|---:|
| MP-02 | Kimi | Kimi K3 | Cline, Azure | not set | 35 |
| MP-09 | Kimi | Kimi K3 | Cline, Azure | not set | 85 |
| MP-15 | Kimi | Kimi K3 | Cline, Azure | not set | 43 (no answer) |
| MP-03 | GLM | GLM-5.3 | OpenCode | not set | 46 |
| MP-05 | GLM | GLM-5.3 | OpenCode | not set | 67 |
| MP-32 | GLM | GLM-5.3 | OpenCode | not set | no answer |
| MP-07 | DeepSeek | DeepSeek-V4.1-Flash | Cline, Azure | off (configuration error) | no answer |
| MP-18 | DeepSeek | DeepSeek-V4.1-Flash | Cline, Azure | off (configuration error) | 107 |
| MP-20 | DeepSeek | DeepSeek-V4.1-Flash | Cline, Azure | off (configuration error) | 36 |
| MP-40 | DeepSeek | DeepSeek-V4.1-Flash | Cline, Azure | high | 96 |
| MP-44 | DeepSeek | DeepSeek-V4.1-Flash | Cline, Azure | high | 46 |
| MP-47 | DeepSeek | DeepSeek-V4.1-Flash | Cline, Azure | high | 28 |

The problem had a 170-minute budget and a 3-hour timeout. The three runs with reasoning off were followed by three
reruns with reasoning on (the case page, Section 4). Minutes and grades as on the case page (Appendix A).

## Credit and status

Problem and paper: H. Gao, L. Zhang, W. D. Nix, C. V. Thompson and E. Arzt. Answers in `runs/` were written by the AI
agents named above. The problem statement, the key and the checker were prepared by the round's coordinating agent, as
the round's records state; the blind grades are by Codex. The case is maintained by Teng Zhang (Syracuse University); its
checks and records were prepared with AI assistants, as each record states. No expert has reviewed the key or the
grades. Corrections are welcome as issues in the MechanicsHarness repositories.

Licence: the code (the `.py` files) is under the MIT License (`LICENSE`); the text and data are under CC BY 4.0
(`LICENSE-CONTENT`). The source paper is not included and is not covered.
