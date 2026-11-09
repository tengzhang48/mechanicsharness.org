# P015 — HRR field and notch (probe tasks T1–T3)

Verdict: **T1 READY** (small fixes optional) · **T2 READY AFTER FIXES** (label as a variant; fix the rubric) ·
**T3 RESTATE** (invented lab narrative, missing load, cues; draft restatement from the paper in this folder).
Overall: RESTATE, because of T3 only.

Source paper: J. W. Hutchinson (1968), "Singular behaviour at the end of a tensile crack in a hardening material",
J. Mech. Phys. Solids 16, 13–31. Text layer and page images: `<probe-runs>/P015-build-codex/paper.txt`,
`pages/page-NN.png` (page-NN = journal p. 12+NN). Used here: Eq. 1 and the no-unloading restriction (p. 14), Eqs. 21–24
(p. 20), Eqs. 26–27 and Table 1 (p. 21), Eqs. 29–31 (p. 22), Eqs. 32–34 (p. 23), Sec. 7 (pp. 26–28), Sec. 8 (p. 29),
Appendix and Fig. 9 (pp. 30–31).

Files:
- Problems: `<probe-runs>/P015-build-codex/T{1,2,3}/{TASK_PROBLEM.md, REFERENCE.md, RUBRIC.md, reference.py, reference_results.json}`, `BUILDER_NOTES.md`, `TASK.md` (builder brief).
- Solver header: `<corpus-pilot-folder>/probes/setup_probe_runs.sh` (head_common + rules_with/rules_without). Every P015 run `TASK.md` is byte-identical to that script's output (rebuilt and compared here).
- Knowledge unit ("with" condition): `<corpus-pilot-folder>/knowledge_v1/P015_hutchinson_1968/{ENTRY.md, EXTRACT.md, REPORT.md, hrr.py, test_hrr.py}`; the copies in the six "with" runs are byte-identical.
- Runs: `<probe-runs>/runs/P015-T{1,2,3}-{with,without}-{glm,dsaz}`; logs `<probe-runs>-logs/`; blind grades `<probe-runs>/grading/R{11,35,05,31,22,32,09,23,19,27,06,10}/GRADE.md` (map: `LABELS_coordinator_only.json`); scores `scores_P015.json`; grader brief `GRADER_INSTRUCTIONS.md`.
- Earlier check (lead): `<library-pilot-folder>/problems_open_book/check_sections/hutchinson1968.md`.
- This audit: `recompute_key.py` (+`_results.json`), `mutant_values.py` (+`.json`), `check_P015.py` (+`_selftest.json/.out`), `scan_T2_branches.py` (+`.json/.out`), `fig9_consistency.out`, `lib/hrr_unit_copy.py` (byte copy of the unit's hrr.py, used only for the plane-stress mutant), `T3_TASK_PROBLEM.md.proposed`, `T3_KEY.md.proposed`.

Not in scope here: the later matched-test problems V1, V2, A1 built from the same paper
(`<runs>/library/P015-build-codex/{V1,V2,A1}/`, runs `<matched-runs>/runs/P015-*`). The GPT draft `hrr-notch.html` mixes them with T1–T3.

## 1 Statement

**T1 (`T1/TASK_PROBLEM.md`).** The paper's own plane-strain problem (Sec. 5, Eq. 31, p. 22) at n = 6.5, inside the
paper's tabulated range n = 3–13 (Table 1, p. 21).
- Complete: crack geometry and faces (line 3), fully plastic J₂ law with all factors and the instruction not to add
  elasticity (lines 5–12), W, the J contour and U per unit thickness (lines 14–20), requested p, H, C with
  normalization and tolerance (lines 24–26). Only the method is open.
- Parameters, not narrative: n, α = 0.02 (the paper's Fig. 1 value, p. 14), σ₀ = 400 MPa, E = 200 GPa are parameters of
  a stated problem and cancel from all three answers.
- Cue (minor): line 30 "Cartesian differentiation of displacement must include rotation of the polar basis" warns
  about the pitfall that decides H and C (the knowledge unit needed a correction round for exactly this, FREEZE.md
  "Coordinator checks"). It helps; it does not give the answer.
- Source: not named in the statement (builder brief `TASK.md` forbade it). The problem is recognizable anyway: all
  three graded without-unit answers on T1/T2 and both without-unit T3 answers name "HRR" (grep of `grading/R{31,09,23,06,10}/submission/ANSWER.md`).
  For release, the source must be cited in the problem's metadata.
- No grader contested any T1 item (R11, R35, R31: "No reference error was found").

**T2 (`T2/TASK_PROBLEM.md`).** A variant: the paper applies its method to one V-notch only (90° opening, n = 3, plane
strain, σₑ ~ r^−0.225; Appendix p. 31, Fig. 9 on p. 30). T2 asks 120° opening, n = 4.
- Complete: material wedge −2π/3 < θ < 2π/3 stated twice (line 3), full elastic + power-law law with all
  parameters (lines 5–16), finite-energy condition (line 16), p and R with sign convention (lines 18–21).
- The solver must choose the inner (hardening-dominated) limit; the statement does not say which; this is the test.
- "the leading local asymptote" (line 18) would be ambiguous if several symmetric eigenvalues existed; our scan
  (section 2) finds only the key's branch over the whole finite-energy range, so the wording is adequate.
- Variant status: the answer is not in the paper; the key is our computation (PROBLEM_CHECK.md "Decisions for Teng" 1).

**T3 (`T3/TASK_PROBLEM.md`).** Not the paper's problem; three invented engineering scenarios.
- Invented narrative (lines 3, 13, 15, 17, 19): "An engineer wants ...", specimens with crack half-length, ligament,
  thickness, tip radius and microstructural length; a "measured" plastic-zone extent; "Local measurements ... show
  plastic strains several times the elastic strains, total strain below 0.03 ... less than 2° change in normalized
  deviatoric stress direction"; 100 load cycles with "measured" reversed plastic strain and closure and a 0.1 mm advance;
  "contour measurements give J = 11.8, 12.2, and 12.0 kJ/m²"; "Measured strain ... reaches 0.06". None of this is in
  the paper, which has no specimens or data. This is the open-book (rule 0) failure the author has flagged.
- Missing data: A asks to "Estimate the stress at r = 0.08" (line 15) and B the "peak-load stress" (line 17), but no
  load, K_I or geometry factor is given. Three of four solvers invented an amplitude by inverting a plastic-zone
  formula, with different formulas (R19 ANSWER.md:33 and R06 ANSWER.md:9: K_I ≈ 26 MPa√m; R10 ANSWER.md:18: 582 MPa√mm
  = 18.4 MPa√m); R27 ANSWER.md:15 took the strain bound 0.03 as a measured value. Graders flagged each as unsupported
  (R19, R06, R10 GRADE.md "Scientific errors"; R27 GRADE.md "(a)").
- Cues: line 3 lists the three outcomes ("can be used, should be replaced, or needs conditional use") and names the
  distinction that decides C ("distinguish the choice of a local field from the choice of a global load-to-amplitude
  relation"); line 13 "This curve supplies no cyclic constitutive data" decides B; line 21 "an unsupported refusal is
  not [welcome]". All four runs returned the reference labels (A use, B replace, C conditional).
- [Left out for this release (8 November 2026): the audit's summary of what the paper supports for a restated T3, the basis of a new problem that has not been run.]


## 2 Key

All values recomputed here with our own code (`recompute_key.py`), written from the paper's Eq. 31 (p. 22, page-10.png),
the plane-strain displacements of p. 23 and the J integral of Eq. 24 (p. 20), not from the builder's `reference.py`.
Formulation: mixed first-order system with constitutive inversion; shooting on (s, f″(0)).

Validation against the paper (Table 1, p. 21, plane strain; our value / paper):

| n | I (max σ̃ₑ = 1) | Eₙ |
|---|---|---|
| 3 | 5.5073 / 5.51 | 1.2272 / 1.23 |
| 5 | 5.0235 / **5.01** | 1.3816 / **1.37** |
| 9 | 4.6018 / 4.60 | 1.4850 / 1.48 |
| 13 | 4.4010 / 4.40 | 1.5184 / 1.52 |

The n = 5 entries of the paper are 0.27% (I) and 0.85% (Eₙ) off; the frozen unit records the same gap (FREEZE.md
caveat 3). Paper erratum to state, not a key error: T1 is solved at n = 6.5 directly, not interpolated.
With s left free, the solver returns s = (2n+1)/(n+1) (Eq. 26) to machine precision at every n.

| Task | Quantity | Key (`reference_results.json`) | Ours | Origin |
|---|---|---|---|---|
| T1 | p | 0.133333333333 | 2/15 | exact, Eq. 26 |
| T1 | H | 1.893569725213 | 1.893569726 (Δ 1.1e−9) | computation (eigenproblem + Eq. 24) |
| T1 | C | 0.457305783974 | 0.457305786 (Δ 2e−9) | computation; C = Eₙ/π (Eq. 33–34) |
| T2 | p | 0.156175260286 | 0.156175260286 (Δ <1e−12) | computation (variant) |
| T2 | R | 0.373262313312 | 0.373262313312 (Δ <1e−12) | computation (variant) |

- T1 diagnostics: I(6.5) = 4.818379 (between the paper's 5.01 at n = 5 and 4.60 at n = 9), σ̃θθ(0) = 2.33526 (max σ̃ₑ = 1),
  Eₙ = 1.43667. n = 1 limit: H = √(2/(3π)) and C = 1/4 reproduced to 1e−9.
- T2 checks: the paper's own notch, 90° opening at n = 3, gives p = 0.22502 against the paper's 0.225 (p. 31). Fig. 9
  (p. 30), read off the figure (approximate, ±0.02): σ̃θ(0) ≈ 1.53, σ̃r(0) ≈ 0.83, σ̃ₑ(0) ≈ 0.65, max σ̃rθ ≈ 0.56; ours
  1.532, 0.811, 0.624, 0.560 (`fig9_consistency.out`). n = 1 gives the elastic wedge root (p = 0.384269, R = 0.444272;
  sin 2βλ + λ sin 2β = 0); β = π at n = 4 gives the crack value p = 0.2.
- T2 branch: `scan_T2_branches.json` (80 × 121 grid, s ∈ [1.600, 1.995] i.e. 0.005 ≤ p ≤ 0.40, the whole finite-energy
  range; f″(0) ∈ [−8, 4]). For every s, f(β) = 0 has exactly one root in f″(0) (80 roots), and f′(β) changes sign once
  along that curve: one symmetric eigenpair, the key's (p = 0.1561753, R = 0.3732623, f > 0 inside). "Leading" is
  therefore unambiguous in this window.
- T3 key (`T3/REFERENCE.md:5-9`): judgements. Parts go beyond the paper: ratios to tip radius and microstructural length,
  r/(J/σ₀) = 5, the 6% strain as a finite-strain marker, constraint from plastic zone/thickness, and "contour agreement is
  not proof" (REFERENCE.md:9, 19). The paper says only that the dominant singularity represents the solution in "a very
  small portion of the plastic zone" (p. 26) and that 3D effects are prominent (p. 28). Its numerical "checks"
  (`T3/reference_results.json`) are length ratios of the invented data and a uniaxial unloading example; nothing to
  recompute against the paper.

## 3 Checker

There is no `checks.py`; all three tasks were graded by a Codex grader against `RUBRIC.md` (`GRADER_INSTRUCTIONS.md`).
`check_P015.py` (this folder) applies the rubric's numeric bands (T1/RUBRIC.md:7,9,10; T2/RUBRIC.md:9,10) to a
`{"T1": {p,H,C}, "T2": {p,R}}` file.

(a) Key-built results: all five numbers "full" (`check_P015_selftest.out`).

(b) Saved runs (numbers from the graded `ANSWER.md`): the checker's bands agree with the graders' numeric rows in all
seven graded cases — R11, R35 (T1 all full); R31 (T1 p full, H and C fail: 2.15092, 0.29178); R32, R09 (T2 both full);
R23 (T2 both fail: 0.08764, 0.3330). Two runs have no answer (R05 T1-without-glm, R22 T2-with-glm: timeouts, scored 0).

(c) Mutants (`mutant_values.json`, all computed here). With "all items full" as the pass rule every mutant fails:

| Mutant | Values | Result |
|---|---|---|
| T1 neighbouring n = 6 / n = 7 | p 0.1429 / 0.125, H 1.8347 / 1.9464, C 0.4525 / 0.4614 | p, H fail; C partial |
| T1 plane-stress field (mean of n = 6, 7) | H 0.976, C 0.815 | H, C fail |
| T1 C reported as Eₙ = πC | C 1.4367 | C fail |
| T1 n/(n+1) dropped from W in U | C 0.5277 | C fail |
| T1 I from linear interpolation of Table 1 (n = 5, 9) | H 1.8916, C 0.4492 | **H full**, C partial |
| T1 I and Eₙ from the n = 5 table entry | H 1.8837, C 0.4361 | H partial, C fail |
| T1 mixed normalization (σ̃ max-scaled, I not) | H 2.7334 | H fail |
| T1 H = σ̃θθ(0)/I (power 1) | H 0.4847 | H fail |
| T1 linear-elastic K field | p 0.5, H 0.4607, C 0.25 | all fail |
| T1 runs: R31 graded copy; its post-deadline overwrite | (2.15092, 0.29178); (1.4706, 1.6978) | H, C fail |
| T2 crack exponent (n = 4) | p 0.2, R 0.7931 | both fail |
| T2 elastic wedge | p 0.3843, R 0.4443 | both fail |
| T2 paper's 90°, n = 3 | p 0.2250, R 0.5294 | both fail |
| T2 60° opening (faces at ±150°) | p 0.1931, R 0.6208 | both fail |
| T2 90° opening at n = 4 | p 0.1810, R 0.5109 | both fail |
| T2 n = 3 at the right angle | p 0.1924, R 0.3801 | p fail, R partial |
| T2 R inverted | R 2.679 | R fail |

(d) Tolerances. T1: p 1e−5 against an exact 2/15; H 0.005 (0.26%) and C 0.002 (0.44%) against a key reproduced to
2e−9. H depends on I only through I^(−2/15), so an error in J of up to about 2% passes H (the Table-1 interpolation
mutant does); C is ∝ 1/I and catches it. T2: p 3e−4 and R 0.002 against 1e−12 agreement; no plausible wrong branch
lands inside. Partial bands (T1 H 0.05, C 0.02; T2 R 0.02) admit neighbouring-n and wrong-n answers (rows above).

(e) Not covered by any numeric check: T1 "brief method, and one verification" (line 28); T2 "short physical
justification ... and one quantitative verification" (line 23); all of T3.

Rubric audit:
- T1, T2 numeric rows are decidable and agree with the key. The method rows are judgement, and the total does not
  separate right from wrong: R31 scored 55 with H and C wrong; R23 scored 76 with both T2 numbers wrong.
- T2 (grade): the rubric has no row for correct governing equations. R23 used the wrong radial exponent in equilibrium
  (R23 GRADE.md "(a)") and still received 30/30 constitutive, 20/20 geometry/branch (including "5 for the leading
  tensile branch" for a wrong eigenvalue) and 15/15 verification. A correct answer with terse justification and no
  verification can score below that (30 numbers + 20 + 20 + 1 ≈ 71 < 76).
- T3: every row is a judgement (e.g. R10's C was capped to 12/30 because the grader read its "conditional" heading as
  a blanket refusal, R10 GRADE.md "C: decision"). An automatic checker could test only the three labels, and the
  statement hands them over (line 3). Deductions rest on points beyond the paper: all four graders took 2 points on
  "C: amplitude" for treating contour agreement as proof of path independence (R19, R27, R06, R10 GRADE.md).

## 4 For outside solvers

- Received: `TASK.md` = header (`setup_probe_runs.sh` head_common + rules) + `T?/TASK_PROBLEM.md`; in the "with"
  condition also `knowledge/` (ENTRY, EXTRACT, REPORT, hrr.py, test_hrr.py). No absolute paths in any of these.
- Output: free-form `ANSWER.md` + `PROCESS.md`. There is no machine-readable result; add a JSON block
  (`{"p":…, "H":…, "C":…}` for T1, `{"p":…, "R":…}` for T2) so `check_P015.py`-style checks can run.
- Budget and tools: 45 min stated, 50 min hard timeout (`launch_run.sh`, 3000 s); python3 with numpy/scipy; no web.
  Three of twelve runs hit the limit (exit 124); two of them (GLM, T1-without and T2-with) left no answer.
- Grading: blind Codex grader with `GRADER_INSTRUCTIONS.md`, `RUBRIC.md`, `REFERENCE.md`, `reference.py`; the grader
  may rerun the submission's code. Blinding is partial: all five with-unit answers mention the knowledge unit.
- Record integrity: `runs/P015-T1-without-dsaz/ANSWER.md` was overwritten at 18:30:27Z, 7 minutes after that run's
  timeout (exit 124 at 18:23:29Z, `probes-logs/P015-T1-without-dsaz.exit`), by the run's own helper session
  (its PROCESS.md, 18:30:10 lines). The graded copy is `grading/R31/submission/ANSWER.md` (sha256 ca82d4…, listed in
  `grading_SHA256SUMS_P015.txt:113`). Publish the graded copy; mark the run-folder copy as post-deadline. No other
  P038/P059/P015 probe run has files newer than its exit time (checked for all 36).
- Knowledge unit defect (not used by T1–T3): the frozen `hrr.py` plane-stress shooting fails at non-integer n
  (n = 6.5 returns I ≈ −1.5×10¹¹); it works at n = 5, 6, 7, 9.
- Public: `T1/TASK_PROBLEM.md`, `T2/TASK_PROBLEM.md`, a restated T3 once approved, the header text, `GRADER_INSTRUCTIONS.md`,
  an output schema, and the scores. Knowledge unit ENTRY/REPORT/hrr.py/test_hrr.py as the "with" input; `EXTRACT.md`
  reproduces the paper's content and falls under the source-paper rule (author decides).
- Hidden: `T{1,2,3}/REFERENCE.md`, `RUBRIC.md` (contains key numbers), `reference.py`, `reference_results.json`,
  `__pycache__/`; `BUILDER_NOTES.md` (key values), build `PROCESS.md`, `scratch/`; all `grading/R*`; all `runs/P015-*`
  (answers contain the key to 10 digits); `LABELS_coordinator_only.json`; this audit's scripts and outputs.
  Never shipped: `paper.txt`, `pages/`.

## 5 Fixes

1. [Left out for this release (8 November 2026): the audit's summary of what the paper supports for a restated T3, the basis of a new problem that has not been run.]
2. T2/RUBRIC.md:8 "5 for the leading tensile branch with finite tip energy" → "5 for the leading tensile branch with
   finite tip energy, only if p is within the 0.003 band"; add a row "governing equations (equilibrium + compatibility)
   correct: 10", taken from the 30 constitutive points. Reason: R23 (wrong equilibrium) scored 76.
3. T1, T2 statements: append "Also report the numbers as JSON: …" (fields above). Reason: automatic checking.
4. Release metadata for T1: "Source: Hutchinson (1968), Sec. 5, Eqs. 26, 31–34, Table 1; new n = 6.5". For T2:
   "Variant of the Appendix (p. 31) V-notch method; key computed by us (two routes) and checked against the paper's 90°,
   n = 3 exponent 0.225 and Fig. 9". Reason: open-book rule, variant labelling (Teng decides on variants).
5. Report correctness separately from the rubric total (numbers pass/fail; method rows as a profile). Reason: totals of
   55 and 76 for wrong answers.
6. Optional, T1/TASK_PROBLEM.md:30: delete "Cartesian differentiation of displacement must include rotation of the polar
   basis." if the author wants the J pitfall to be part of the test (it then tests more, and is harder).
7. Publish `grading/R31/submission/ANSWER.md` as the record of T1-without-dsaz, not the run-folder copy.

## 6 Not verified

- Session logs (`probes-logs/*.jsonl`) were not reread for rule compliance (no web, no outside files); the pilot summary
  says they were audited.
- The builder's T3 scale-ratio arithmetic was not rechecked (it describes invented data).
- Uniqueness of the T2 branch is shown only within the scanned window of (s, f″(0)); no proof.
- The V1, V2 and A1 matched-test problems from the same paper were not examined.
- Whether `EXTRACT.md` may be published.
