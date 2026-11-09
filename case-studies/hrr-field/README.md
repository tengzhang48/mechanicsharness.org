# HRR field and notch: can AI agents rebuild and adapt it?

A MechanicsHarness case. The case page explains the mechanics, the solution and what each agent did:
https://mechanicsharness.org/case-studies/hrr-field.html

The problems come from J. W. Hutchinson (1968), "Singular behaviour at the end of a tensile crack in a hardening
material", *J. Mech. Phys. Solids* 16, 13–31: the singular field at a tensile crack in a power-law hardening material
(the HRR field), its amplitude from Rice's J integral, and the paper's V-notch appendix. The agents received the problems
without the paper. The paper is not included here.

This folder lets you attempt the four problems, check a new answer, and check every saved answer and grade against the
keys.

## The problems in brief

Probe round, 6 October 2026 (this folder's top level):

- **T1**, the paper's own plane-strain crack problem at a hardening exponent its table does not list, n = 6.5, with the
  fully plastic law alone: the stress exponent p, the amplitude coefficient H in σθθ(r, 0) = H σ₀ [J/(α ε₀ σ₀ r)]^p and
  the energy coefficient C = U(R)/(J R). Key: p = 2/15, H = 1.893570, C = 0.457306.
- **T2**, a variant: elastic plus power-law material at n = 4 with a 120° notch (the material fills −2π/3 < θ < 2π/3):
  the exponent p of the leading symmetric field and R = σrr/σθθ ahead of the notch. Key: p = 0.156175, R = 0.373262
  (our computation; the paper solves only a 90° notch at n = 3).

Matched test, 7 October 2026 (`matched_test/`):

- **V1**, a variant: an antiplane crack in a pure power-law material in shear, n = 3.5, with a given J: the face slip
  at two radii, the signed stresses at one point and three radii of the |γ| = 0.01 contour.
- **V2**, adapted from the paper's piecewise-linear material: a plane-stress membrane whose tangent slope T is to be
  inferred from a given transverse-strain coefficient, then three strain coefficients, a displacement coefficient and a
  screened interval. Its statement has three flaws, named in `matched_test/FIXES.md`.

The keys of V1 and V2 are our computations, recomputed independently by the audit of 7 November 2026. Two further
problems from the same paper are not here: T3 of the probe round and A1 of the matched test rested on invented
measurements (the case page, Section 12).

## What is here

| Path | Content |
|---|---|
| `task/T1`, `task/T2`, `matched_test/task/V1`, `matched_test/task/V2` | The problems exactly as given. Each run's `workspace/TASK.md` is the whole task it received, with the shared header; runs with knowledge also have `workspace/knowledge/`. |
| `key/T1`, `key/T2` | The references, rubrics and reference scripts of T1 and T2 as they stand (rubrics fixed on 26 October 2026); `RUBRIC_2026-10-06.md`: each rubric as the graders had it, superseded (a notice at the top says so). |
| `key/check_P015.py` | The numeric checker of T1 and T2 (added on 26 October 2026). |
| `FIXES.md` | What the audit of 25 October found in T1 and T2 and what was fixed. Read it before using the keys. |
| `audit_2026-10-25/` | That audit: its recomputation of the keys, its mutants, its scan of T2's branches and its checker run on the graded answers. Its paragraphs on a restatement of T3 are left out (`RELEASE_EDITS.md`), because a restated T3 has not been run yet. |
| `regrade_2026-10-26/` | The checker run on the saved T1 and T2 answers, with a summary. |
| `runs/<run>/workspace/` | Every saved file of the eight probe runs, as delivered; `runs/<run>/session/`: their start, exit and error stamps. |
| `grading/` | The blind grades of the six answered probe runs, with the graders' scratch work, and the key builder's notes. |
| `matched_test/key/V1`, `matched_test/key/V2`, `matched_test/key/check_V1_V2.py` | References, rubrics and reference scripts of V1 and V2, unchanged, and their numeric checker (added on 8 November 2026). |
| `matched_test/FIXES.md`, `matched_test/audit_2026-11-07/` | The audit of V1 and V2 (both keys right; V2's statement flaws named), its independent recomputations, its mutants and its checker run on the 24 graded answers (`check_V1_V2_audit.out`). |
| `matched_test/runs/<run>/` | The 24 graded runs of GLM and DeepSeek and the four Kimi K3 runs on V1 that saved an answer or a draft. |
| `matched_test/grading/` | The blind grades (`grade-M*`) and codings (`code-M*`) of the 24 graded answers; the label map is `matched_test/records/matched_results.json`. Also the key builder's notes and two reviews of the knowledge topic. |
| `records/`, `matched_test/records/` | Records cited by the case page: the probe pilot summary and the process log, cut to this case (each says inside what was kept), the scores, the grader instructions, the launch and queue scripts, the matched-test summary, and the start, exit and error records of the twelve Kimi launches that a token-rate limit stopped (`kimi_launches/`). |
| `reference_checks/` | Our check of the reference values (`hrr_check.py`, output `hrr_check.txt`, which the page's figure reads). |
| `check_package.py` | Checks this folder (below). |
| `RELEASE_EDITS.md` | Every change made to a file for this release. |
| `SHA256SUMS` | Hashes of every file. |
| `LICENSE`, `LICENSE-CONTENT` | MIT for the code; CC BY 4.0 for the text and data. |

Not included: session transcripts of the agents and graders, the source paper, the four probe attempts that were cut
off when a daemon was stopped during a disk-full incident and were rerun, and T3 and A1.

## Check an answer

T1 and T2: write `{"T1": {"p": ..., "H": ..., "C": ...}, "T2": {"p": ..., "R": ...}}` (either part may be left out) to a
file and run

    python3 key/check_P015.py results.json

It prints, for each number, `full`, `partial` or `fail` against the rubric's bands (full credit within T1: p 10⁻⁵,
H 0.005, C 0.002; T2: p 3 × 10⁻⁴, R 0.002; absolute) and whether all are at full credit. A value that is not a number
is reported as such. An input without a T1 or T2 block exits with status 1.

V1 and V2: write `{"V1": {"slip_0.02", "slip_0.08", "tau_rz", "tau_tz", "r_0", "r_90", "r_180"}, "V2": {"T", "C_rr",
"C_tt", "C_g", "C_y", "r_min", "r_max"}}` with numbers in the problems' units and run

    python3 matched_test/key/check_V1_V2.py results.json

It applies the stated tolerances (V1: all 2 %; V2: T 1 %, the others 2 %; relative, signed) and exits with status 0
only if all seven numbers of each block given pass (a missing number fails; an input without a V1 or V2 block fails). Count an answer as right only if all seven numbers of a problem pass: some single
numbers pass for wrong models (`matched_test/FIXES.md`). Python 3, standard library only, for both checkers.

## Check this package

    python3 check_package.py

It copies this folder to a temporary directory, checks every file against `SHA256SUMS` (and that no other file is
present), reruns both checkers' self-tests (the saved T1 and T2 answers and the audit's mutants; the 24 graded V1 and V2
answers and the audit's mutants) and compares their output with the saved one, and runs the checkers on the keys, on
values just inside the bands, on deliberately wrong answers (the crack and elastic-wedge exponents at the notch, plane
stress and a missing n/(n+1) for T1, the tensor-shear and dropped-traction mistakes for V1, the incompressible and
plane-strain mistakes for V2) on a value given as text and on inputs with no problem block. Its last line is `OK: all checks passed`. Python 3,
standard library only.

## What needs human judgment

The checkers compare numbers. Whether a run derived the field, chose the right branch and verified its result, rather
than recalling the HRR result, needs a reader, and the rubrics' method, verification and scope rows are judgements. The
blind grades in `grading/` and `matched_test/grading/` are by Codex, as written at the time; the Kimi answers were not
graded. The rubrics' consequential credit gave several answers high totals with their decisive numbers wrong, so report
whether all numbers are within tolerance beside any total (both `FIXES.md`).

## The runs

| Round | Agent | Model | Interface | Runs |
|---|---|---|---|---|
| Probe round, 6 Oct | GLM | GLM-5.3 | OpenCode | 4 (T1, T2; with and without the knowledge unit); 2 reached the 50-minute limit without an answer |
| Probe round, 6 Oct | DeepSeek | DeepSeek-V4.1-Flash | Cline, Azure, reasoning off | 4 (the same cells) |
| Matched test, 7 Oct | GLM | GLM-5.3 | OpenCode | 12 (V1, V2; conditions B0, B1, C; two repeats) |
| Matched test, 7 Oct | DeepSeek | DeepSeek-V4.1-Flash | OpenCode, reasoning on | 12 (the same cells) |
| Matched test, 7 Oct | Kimi | Kimi K3 | Cline, Azure | 16 launches; 3 complete V1 answers and 1 draft, not graded; 13 launches stopped by a token-rate limit |

Conditions in the matched test: B0 the crack-tip topic, its checks and the paper's knowledge entry, without solver code;
B1 the same plus the entry's reference code; C no knowledge (`matched_test/records/MATCHED_TEST_SUMMARY.md`). Each run had
a 45-minute budget and a 50-minute limit. Minutes, grades and numbers within tolerance for every run are in the case
page's Appendix A.

## Credit and status

Paper: J. W. Hutchinson. Answers in `runs/` and `matched_test/runs/` were written by the AI agents named above. The
problems, keys and rubrics were prepared by Codex from the paper; the blind grades are by Codex; the audits and checkers
by AI assistants, as each states. The case is maintained by Teng Zhang (Syracuse University); its checks and records
were prepared with AI assistants, as each record states. No expert has reviewed the keys or the grades. Corrections are
welcome as issues in the MechanicsHarness repositories.

Licence: the code (the `.py` files) is under the MIT License (`LICENSE`); the text and data are under CC BY 4.0
(`LICENSE-CONTENT`). The source paper is not included and is not covered.
