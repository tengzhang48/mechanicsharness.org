# V — Steady crack growth in a viscoelastic strip (Mode III, standard linear solid, Dugdale zone)

Verdict: READY AFTER FIXES

Source paper: C.-Y. Hui and R. Long (2026), "Steady state crack growth in viscoelastic solids: Full-field analysis for the
strip geometry with finite damage zone", J. Mech. Phys. Solids 214, 106662, and its supplement
(`hui_long_2026-10-22/private/sources/Hui_Long_2026_{main,supplement}.pdf`).

Files (folder `<pilot-folder>/`):
- Statement (all a solver gets): `agent/PROBLEM.md`, `agent/CONTRACT.md`, `agent/case.json`. No header file; no knowledge
  entry.
- Key: `private/REFERENCE.md`, `private/reference_results.json` (from `private/reference_solver.py`),
  `private/published_fig9.json`, `private/verification_results.json`.
- Checker: `private/score_submission.py` (numeric gates), `private/SCORING.md` (100-point rubric),
  `private/SCORING_ANALYTICAL_V2.md` (70/20/10, added after launch, README.md:13).
- Runs: `runs/P01_deepseek`, `P02_kimi` (no submission), `P03_gpt61_sol`, `P04_deepseek_direct_1017`, `P05_kimi_recheck`
  (no submission); reviews `runs/R01..R05_*`.
- This audit: `audit_visco.py` (own solver, not the reference code), `audit_visco.out`, `audit_visco_output.json`,
  `our_recomputation.json`, `score_submission_copy.py`. The paper was read through a pdftotext dump and a page-19
  render made for this audit and deleted afterwards (publisher's text is not kept here).

## 1 Statement

Checked line by line against the paper (pages 4-6, 12-14, 18-19 of the PDF):
- Strip, loading and moving origin — PROBLEM.md:7-17 = Eqs. 1-4 (p.4-5) and Fig. 1. Bounded stress at the zone front —
  PROBLEM.md:25 = §3.2 (p.6). Relaxed far-ahead state — PROBLEM.md:25 = §2.3 (p.5).
- Creep compliance — PROBLEM.md:31-35 = Eqs. 5, 7 (p.5). Free energy — PROBLEM.md:39-43, supplied, which removes the
  paper's nonuniqueness of local dissipation (REFERENCE.md:23, p.25). Good.
- Separation criterion 2w(−l_c, 0⁺) = δ_c — PROBLEM.md:45 = Eq. 48 (p.12). Γ — PROBLEM.md:57 = Eq. 9 (Γ = G∞γ∞²h, the
  stored energy of the full height 2h). Normalizations V = πvt_c/h, ε = πl_c/h, δ̄ = G∞δ_c/(τ_c h) = Eqs. 17, 23, 52.
- λ = 0.1, δ̄ = 0.2, V = 0.001, 0.1, 10 and the streamlines y/h = 0.1, 0.5, 1 are the paper's Fig. 9 case (p.19). Sample X
  values are organizer choices (REFERENCE.md:25).

Every parameter is given in dimensionless form; every requested quantity is defined (T_i, E_i, R, Γ_D, admissible,
critical_speed: CONTRACT.md:13-17). No errata found, no invented narrative. Minor points:
- **minor — cue.** "report any finite endpoint" (PROBLEM.md:56) hints that one may exist; Q3's ban on Γ − Γ₀
  (PROBLEM.md:58) states the energy balance. Both are needed to define the task; acceptable.
- **minor — recognizable.** The title, parameters and streamlines identify Hui & Long 2026; CONTRACT.md:34 asks solvers to
  record recognition. Acceptable for a closed-book reproduction test.
- No solver or reviewer raised a statement ambiguity (GPT ANSWER.md has no Questions section; Grok R02 REVIEW.md:100 finds
  the model "source-consistent").

## 2 Key

| Item | Key | Origin | Recomputed here |
|---|---|---|---|
| ε, a, R at V = 0.001, 0.1, 10 | REFERENCE.md:64-68 | Eq. 51 solved for ε; a = (2/π) acos(e^{−ε/2}) (Eq. 18); R = a²/δ̄ (Eq. 22) | **yes**, own code (mpmath, Eq. 62 for I₁): 0.544705555093 / 0.448842489170 / 1.00729790042; 0.74078093167 / 0.51481218674 / 1.32515793808; 6.33330568444 / 0.973162076321 / 4.73522213395 — equal to the key to ~1e-11 |
| Published check | published_fig9.json | Fig. 9 column headings, p.19 — a **figure annotation**, 3 decimals | **yes**: I rendered p.19 and read 0.545/0.449/1.007, 0.741/0.515/1.325, 6.333/0.973/4.735; transcription correct |
| V* = 13.4235715733, R → 1/δ̄ = 5, a → 1, ε → ∞ | REFERENCE.md:54-62 | Eq. 59 (p.13; "V* = 13.4", p.18) | **yes**: (1−λ)k(V*) = 1 − δ̄ with the series of Eq. 59 → 13.42357157330; the series agrees with term-by-term integration of Eq. 58 (by hand) |
| 180 field samples T_i, E_i | reference_results.json | T: Eq. 12; E: Boltzmann convolution of Eq. 38 (REFERENCE.md:74-81) | **yes**: own evaluation, max difference 4.1e-14 |
| Γ_D/Γ₀ | reference_results.json | spatial integral, Eqs. 83-86 | reference integral equals R − 1 to ≤ 2e-7 (Eq. 11 energy balance); I did not redo the 2D integral |

Slow limit (key, REFERENCE.md:52): R → 1, ε → ε₀ = −2 ln cos(π√δ̄/2) = 0.5406 (Eq. 65). All values are computed from
closed forms and quadrature; none is a figure reading.

## 3 Checker

(a) Key file scored against itself: all gates pass. Our independent recomputation also passes.
(b) Saved runs rescored with a copy of `score()`: P01 fails 8 checks, P04 fails 31 — identical lists to
`SUBMITTED_NUMERIC_REPORT.json`. P03 (GPT) passes all; its automatic score was withheld (EVALUATION_STATUS.json), so there
is nothing to compare. P02, P05: no results.json.
(c) Mutants (`audit_visco.out`):

| Mutant | Result |
|---|---|
| elastic (R = 1, ε = ε₀, E = T, unbounded) | fails 34 |
| V without π (V* = 4.27) | fails critical speed |
| unbounded speed (null) | fails |
| V = 10 declared inadmissible (P04's error) | fails |
| criterion w(−l_c) = δ_c (factor 2, δ̄ → 0.1 = λ) | fails 12 |
| τ₁ sign flipped | fails 46 |
| E = λT (instantaneous only) | fails 65 |
| E = T (no history) | fails 23 — but **passes every V = 0.001 field** (slow limit, physics) |
| Γ_D = R − 1 by subtraction | **passes** (scorer cannot see provenance; SCORING.md:14 leaves it to review) |
| Γ_D = 0 at V = 0.001 | **passes** both the bulk gate and the energy-balance gate |
| ε 0.4 % low everywhere | passes (inside 0.5 %) |

(d) Bands: 0.5 % on ε, a, R, V* (CONTRACT.md:25) is far above the key's precision and below every plausible model error
tested. Field band 0.005 absolute: field values range 0.01-1, so it is 0.5-50 % relative; at y/h = 1 T₁ = E₁ = 0 by the
clamped edge (trivial items). Bulk band 0.01·max(1, Γ_D/Γ₀): at V = 0.001 the key is 0.0073, so any value in
[−0.0027, 0.0173] passes — zero included, and P01/P04 passed with +22 % and +49 % errors (reviews/QUALITATIVE_FINDINGS.md,
"An absolute tolerance hides…").
(e) Not covered automatically: the derivation; the endpoint limits (R → 5, a → 1, ε → ∞) and the slow limit; integration
bounds, refinement and tail evidence; the Q4 checks; whether Γ_D came from an integral.

Issues:
- **grade — bulk-dissipation gate is vacuous at V = 0.001.** `score_submission.py:61-67` with CONTRACT.md:27: Γ_D = 0
  passes there.
- **minor — no internal-consistency checks.** a = (2/π) acos(e^{−ε/2}) and R = a²/δ̄ are necessary relations (Eqs. 18, 22);
  the P01 inconsistency (2.45 % in a, QUALITATIVE_FINDINGS.md) was found by hand.
- **statement/grade — rubric not frozen.** SCORING.md:40 calls itself "a proposed rubric … then freeze its version";
  SCORING_ANALYTICAL_V2.md was added after launch (README.md:13). A release needs one frozen rubric.

## 4 For outside solvers

- Received: `agent/PROBLEM.md`, `agent/CONTRACT.md`, `agent/case.json` only (README.md:3; PILOT_STATUS.md:13: separate
  containers, no internet, only the model API). Budget 6 h (CONTRACT.md:36). Python + NumPy + SciPy, "symbolic tools already
  present" (CONTRACT.md:9) — the environment is defined by `private/Dockerfile*`, which a release must publish or describe.
- Output: `python solve.py case.json results.json`; schema CONTRACT.md:13-19.
- Grading: `python score_submission.py reference_results.json results.json report.json` (paths as arguments, portable),
  then a reviewer with SCORING.md (gates + 100 points).
- Machine dependencies: `private/isolated_run.py:125,250-251` (reads `<home>/.local/share/opencode/auth.json`,
  mounts local CLI binaries) and `monitor_pilot.py` — organizer harness only; do not publish as is.
- Public: `agent/PROBLEM.md`, `agent/CONTRACT.md`, `agent/case.json`, `private/score_submission.py` (holds no key values),
  SCORING.md once frozen, the Dockerfile.
- Hidden: `private/reference_results.json`, `REFERENCE.md`, `reference_solver.py`, `verify_reference.py`,
  `verification_results.json`, `published_fig9.json`, `test_scoring.py` (reads the key), the runs and reviews. Never
  distribute `private/sources/*.pdf` (publisher's PDFs). Note the key is published in the paper itself (Fig. 9 headings,
  V* = 13.4), so hiding it only protects the digits.

## 5 Fixes

1. CONTRACT.md:27 "`bulk_dissipation_ratio` must have absolute error at most 0.01 max(1, Γ_D/Γ₀). Also require
   |Γ_D/Γ₀+1−R|/max(1,R) ≤ 0.01" → "… must be within 5 % of the reference or 2e-4 absolute, whichever is larger. Also
   require |Γ_D/Γ₀ + 1 − R| ≤ 0.05 (R − 1) + 2e-4". Same change in `score_submission.py:61-67`. Reason: §3 (zero passes
   at V = 0.001; with the new band zero fails by 0.0073 > 3.6e-4). Report the old result alongside for P01-P05 (both would
   then fail the V = 0.001 bulk item).
2. `score_submission.py`: add two principle checks on the submitter's own numbers, |a − (2/π)acos(e^{−ε/2})| ≤ 0.5 % and
   |R − a²/δ̄| ≤ 0.5 % (source Eqs. 18, 22).
3. Freeze one rubric (SCORING.md, with the anchors of SCORING.md:31-36) and record V2 as a post-hoc reading.
4. Release: strip `<home>` paths from the harness or keep it private; publish the container recipe.

## 6 Not verified

- The 2D spatial dissipation integral itself (I relied on its agreement with R − 1, an identity, and did not rerun
  `reference_solver.py --energy-ny 64`).
- The supplement's derivations (S2-S7) — I used the main-paper equations only.
- The Grok reviews' scores (R02 100/100 for P03): not re-graded.
