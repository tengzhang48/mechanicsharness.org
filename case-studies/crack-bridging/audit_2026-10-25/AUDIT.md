# P038 — Bao & Suo 1992, crack bridging: probes T1–T3 (audit 2026-10-25)

Verdict: **T1 READY. T2 READY AFTER FIXES (label as a variant; optional discriminating third output). T3 RESTATE.**
Overall: RESTATE (T3 only); T1 and T2 can go in once the fixes below are made.

Source paper: G. Bao and Z. Suo (1992), "Remarks on crack-bridging concepts", Appl. Mech. Rev. 45, 355–366.
Text layer and 300-dpi page images: `<probe-runs>/P038-build-codex/paper.txt`,
`pages/page-01..12.png` (page-01 = p. 355). Equations below were read from the page images (pp. 357, 359, 361, 362, 364).

Files (B = `<probe-runs>/P038-build-codex`, R = `<probe-runs>`):
- Statements, keys, rubrics: `B/T{1,2,3}/TASK_PROBLEM.md`, `REFERENCE.md`, `RUBRIC.md`, `reference.py`, `reference_results.json`; `B/BUILDER_NOTES.md`, `B/TASK.md`.
- Header and setup: `<corpus-pilot-folder>/probes/setup_probe_runs.sh` (but see §4: the runs got an older header), `launch_run.sh`, `GRADER_INSTRUCTIONS.md`, `build_grading.py`, `LABELS_coordinator_only.json`, `scores_P038.json`.
- Runs: `R/runs/P038-T{1,2,3}-{with,without}-{glm,dsaz}/`; blind grades `R/grading/R{24,08,30,03,26,36,15,21,16,34,18,33}/GRADE.md`.
- Knowledge unit ("with" runs): `<corpus-pilot-folder>/knowledge_v1/P038_bao_suo_1992/`.
- Lead (not a finding): `<study-folder>/library/problems_open_book/check_sections/evans_hutchinson1984_bao_suo1992.md` §P038.
- This audit: `recompute_P038.py` (+`_results.json`), `check_P038.py` (+`check_P038_selftest.{json,txt}`, `runs_extracted.json`), `t2_discriminating_output.py` (+`.json`), `t2_delta5mm_exact.py` (+`.txt`), `ref_copy/` (builder scripts rerun on copies), `T3_TASK_PROBLEM.md.proposed`, `T3_KEY_RUBRIC.md.proposed`.

## 1 Statement

**T1 "A finite flaw under opening tension"** (`B/T1/TASK_PROBLEM.md`, 227 words). The paper's own problem: a central slit 2a in an
infinite plane-strain plate, constant closing traction on a < |x| < a+L (rectilinear law, p. 359 Fig. 1/§2.2), finite-crack
relations Eqs. 45–46 (p. 362), with nonzero front work handled as the paper prescribes ("rising part … by letting K_tip = K0", p. 359;
"subjected to G_tip = Γ0", p. 361). Model, geometry, loading, conventions (full gap, per-end L, area counted once, lines 20–24) and
tolerance (line 35) are all stated; only the method is open. No cue and no source naming. The numbers (E = 130 GPa, ν = 0.28,
a = 2.40 mm, 160 MPa, δc = 9 µm, Γ0 = 22 J/m², lines 4–22) are invented, but as parameters of a stated mathematical problem, not as
measurements; a/(δ0E′/σ0) = 0.302 and Γ0/(σ0δ0) = 0.015 lie inside the range the paper discusses (Figs. 10–11, p. 362; "Γ0/σ0δ0 ≪ 1",
p. 359). No grader contested any item. Minor: the event asked (first connector failure) is not the paper's strength σ_max; the
paper warns σ̄ need not be monotone when G_tip ≠ 0 (p. 361). For these values Σ rises monotonically up to the event (checked, §2),
so the question is well defined.

**T2 "Separation between unequal strips"** (`B/T2/TASK_PROBLEM.md`, 325 words). A **variant** of the paper's splitting beam (§3.4,
p. 364, Eq. 51: equal arms, rectilinear law, E′): unequal arms, linear softening law (the paper's Eq. 21 with s < 0, Fig. 7, p. 359),
Γ0 > 0. The beam idealization is fully prescribed (lines 17–19: linear elastic, small rotation, no stretching or shear, longitudinal
moduli), as are the moment convention (lines 13–15), the gap convention (line 30), the event (lines 37–38) and tolerance (line 45). Well
posed; only the method is open; no cue. Invented content is geometry and parameter values, not data. The paper gives no
solution for this geometry (only the scaling Eq. 53 with a1 from Suo et al. 1992a, not in the paper), so the key is ours.
The 20 mm starter cut (line 11) does not enter the answer (pure couples, zero shear); a harmless distractor.

**T3 "Choosing a fracture description"** (`B/T3/TASK_PROBLEM.md`, 333 words) — **lab narrative with invented data; cues the decision.**
- Invented specimens and measurements presented as observations, none in the paper: "The stated curves come from slow monotonic
  opening tests" (12–13); A "bulk response … is verified linear elastic", connector spacing 2 µm (18–20); B "measured plastic strain
  of 4% across most of the ligament", "A local opening test gives" (30–32); C "their measured opening curve" (42), "In-situ
  observations reliably bound every tunnel-face opening by δ ≤ 3.00 µm" (49–50).
- [Left out for this release (8 November 2026): the audit's account of what the paper supports for a restated T3, the basis of a new problem that has not been run.]
- Minor: B's plate width is undefined ("crack of total length 8 mm and a remaining ligament of 10 mm", 27–28) and "nominal applied
  tension" is not defined (gross or net); the answer does not depend on it.

Paper erratum found (minor, does not affect any key): p. 361 right column, "subjected to G_tip = Γ0 in (35), and δt ≤ δ0 in (36)"
— the equation numbers are swapped ((35) defines δt, (36) defines G_tip).

## 2 Key

| Item | Key | Provenance | Recomputed (this audit) |
|---|---|---|---|
| T1 Σ | 133.730286053811 MPa (`T1/REFERENCE.md:3`, `RUBRIC.md:6`) | paper Eqs. 45–46 with K_tip = √(E′Γ0), δ(a) = δc; builder's computation | **yes**, two routes: Eqs. 45–46 (checked on page image p. 362) and an Abel/Paris double-integral opening (not the builder's log kernel): 133.73028605381 / 133.73028605382 MPa, rel. diff 8e-14; key error 4e-16 |
| T1 L | 4.08560913824553 mm per end | same | **yes**: 4.0856091382455 / 4.0856091382458 mm, rel. diff 6e-14 |
| T1 branch | Σ monotone from initiation (20.288 MPa) to the event | builder | **yes**: Σ and δ(a) strictly increasing on (a, c*]; global branch max 161.01 MPa at L ≈ 240 mm (after the event, irrelevant) |
| T2 M | 0.212625787376721 N·m each strip | our derivation: beam ODE δ'''' = −C t(δ), δ(0)=δ′(0)=0, δ″(0)=√(2CΓ0), δ(L)=δc, δ‴(L)=0; M also from energy Cm²/2 = Γ0 + σbδc/2 | **yes**, own shooting (DOP853 + 2-D root): 0.2126257873767214 (rel. 2e-15); energy route identical |
| T2 L | 10.8001023708973 mm | same | **yes**: 10.800102370897 mm (rel. 9e-15); closed-form route 10.800102370897 |
| T2 δ(L/2) | 38.2188979260933 µm | same | **yes**: 38.21889792610 µm (rel. 2e-13) |
| T2 path | L and M increase with rear opening up to the event; profile monotone | builder | **yes** (continuation over 50 rear openings; all profiles monotone, 0 ≤ δ ≤ δc) |
| T3 | qualitative labels + diagnostics (ℓm = 0.8 mm, W = 180, 4500, 426.67 J/m², cap work 3.0984 J/m² = 0.726% of full, 25.82% of Γ0) | builder's judgement; numbers are arithmetic on invented data | arithmetic **yes** (all match); judgement C **beyond the paper** (below) |

Builder scripts rerun on copies (`ref_copy/`): all three reproduce their saved `reference_results.json` exactly.

The paper's front-work condition and the T2 front condition agree: for equal arms, C = 2·12/(E′h³) gives Cm²/2 = 12m²/(E′h³), the
paper's Eq. 51 left side (p. 364). T2 uses E instead of E′ on purpose (line 19).

[Left out for this release (8 November 2026): the audit's account of what the paper supports for a restated T3, the basis of a new problem that has not been run.]

## 3 Checker

No `checks.py` exists; all three are rubric-graded by blind Codex graders (`GRADER_INSTRUCTIONS.md`). I wrote `check_P038.py`, which
implements the only rows a machine can decide: T1 Σ, L (bands of `T1/RUBRIC.md:14-15`, pass = ±0.3%) and T2 M, L, δ(L/2) (bands of
`T2/RUBRIC.md:26-35`, pass = ±0.5%). Results (`check_P038_selftest.txt`):

- (a) Key-built results and my independent recomputation: all pass, 50/50 and 40/40.
- (b) Headline numbers extracted from the 8 saved T1/T2 answers (`runs_extracted.json`, file:line given): pass/fail and band points
  agree with every grader's numeric rows (R08, R24, R30: 50/50; R03: 0/50; R15, R21, R26, R36: 40/40). No saved automatic flags exist
  to compare (there is no `auto/` for these probes).
- (c) Mutants (all fail at least one item, as expected):

| Mutant | Values | Checker |
|---|---|---|
| T1 seen R03 (wrong opening formula) | Σ 158.525, L 34.324 | fail both; 0/50 |
| T1 Γ0 omitted | 131.835, 6.391 | fail both; Σ −1.4% → 12/25, L +56% → 0 |
| T1 one-face displacement taken as gap | 154.101, 16.186 | fail both; 0/50 |
| T1 total 2L reported | 133.730, 8.171 | Σ pass, L fail; 25/50 |
| T1 plane-stress E | 130.501, 3.616 | fail both; 17/50 |
| T1 Γ0 doubled | 135.496, 3.479 | fail both; 17/50 |
| T1 small-scale-bridging LEFM (no zone) | 165.384, – | fail; 0 |
| T2 Γ0 omitted | M 0.20555, L 12.023, δ 31.013 | fail all; 10.5/40 |
| T2 averaged equal arms | 0.26558, 12.070, **38.219** | M, L fail; **δ(L/2) passes** |
| T2 only upper-arm compliance | 0.23014, 11.236, **38.219** | M, L fail; **δ(L/2) passes** |
| T2 moment per width (N) / sum of moments | 8.505 / 0.42525 | M fail; 25/40 |
| T2 one-arm half of the gap | δ 19.109 | δ fail; 30/40 |

- (d) Tolerances vs physics and key precision: keys are exact to ~1e-13, so ±0.3% and ±0.5% are pure grading choices, well above
  numerical error. T1: Σ alone is weakly sensitive to Γ0 (omitting it moves Σ only 1.4%); L carries the discrimination — both are
  asked, fine. T2: **δ(L/2)/δc depends only on Γ0/(σbδc)**, not on the arm stiffnesses (proved by the mutants above: identical to
  1e-13), so the 10-point δ(L/2) row cannot detect the decisive unequal-arm error. A discriminating replacement is δ at a fixed
  distance, e.g. δ(x = 5.00 mm) = 31.7700837 µm (two routes, `t2_delta5mm_exact.txt`, `t2_discriminating_output.json`: shooting
  and closed form agree to 2e-14); the averaged-arm and one-arm mutants then miss by −23% and −9%.
- (e) Requested items no automatic check covers: T1/T2 "brief method" and "one check" (rubric method/verification rows); all of T3.
  T1's verification row credits any executed elasticity check: R03 got 15/15 for a correct K check while its opening (the decisive
  relation) was wrong (`R/grading/R03/GRADE.md` row Verification) — minor rubric weakness.
- T3 rubric decidability: the labels are decidable, but the rubric says a different label "earns full credit when its scope and
  physical reason are equivalent" (`T3/RUBRIC.md:5-6`) and every item is full/half/zero on "correct reasoning" (lines 9-11), so no
  item is decidable by a script. Graders split on the same kind of content (A "workable next calculation": R33 0/7, R18 3.5/7,
  R16/R34 7/7). T3 scores: 96, 89.5, 88.5, 60.5 (`scores_P038.json`).

## 4 For outside solvers

- **What a solver received:** one file `TASK.md` = common header + rules + `---` + `TASK_PROBLEM.md` (verified: the problem part of all
  12 run `TASK.md` files is byte-identical to the current `TASK_PROBLEM.md`; run hashes in `runs_SHA256SUMS_P038.txt` all OK). The
  header the P038 runs got is **not** the current `setup_probe_runs.sh`: the kill-process rule (script lines 19, 26) was added at
  17:04 UTC, after the P038 runs were built at 16:27 (`R/runs/P038-T1-without-glm/TASK.md`, no "pkill" line). Use a run's
  `TASK.md` as the record of the header. "With" runs also got `knowledge/` = ENTRY.md, EXTRACT.md, REPORT.md, crack_bridging.py,
  zone.py, test_crack_bridging.py (no CORRECTIONS.md or FREEZE.md); the unit contains Eqs. 45–47 and names the paper.
- **Output schema:** none; free-form `ANSWER.md` plus `PROCESS.md`. A release needs a results block (e.g. the JSON schema in
  `check_P038.py`) for T1/T2.
- **Budget and tools:** 45 min stated, 50 min hard timeout (`PROBE_PILOT_SUMMARY.md:19`); python3 with numpy/scipy; no web; one run
  per cell; solvers GLM-5.3 (opencode) and DeepSeek-V4.1-Flash (Cline, Azure) (`launch_run.sh`).
- **Grading:** fresh Codex session per answer, 25-min limit, given TASK_PROBLEM, RUBRIC, REFERENCE, reference.py,
  reference_results.json and the submission (`build_grading.py`); labels blind (`LABELS_coordinator_only.json`), but 6 of 12 P038
  answers (every "with" run) name "Bao & Suo" or the knowledge unit (grep of `R/runs/P038-*/ANSWER.md`;
  `PROBE_PILOT_SUMMARY.md:101` says 4 of 12), so blinding is partial for those.
- **Machine dependencies:** `reference.py` files have no absolute paths and only print JSON. The setup/launch/grade scripts hard-code
  `<volume>/...` and the agent CLIs.
- **Public:** T1 and T2 `TASK_PROBLEM.md`, the restated T3 problem (once approved), the solver header text, `GRADER_INSTRUCTIONS.md`,
  an output schema.
- **Hidden** (hold key values or answers): `T{1,2,3}/REFERENCE.md`, `RUBRIC.md` (key numbers at T1:6-7, T2:28-30, T3:23-24,50-52),
  `reference.py`, `reference_results.json`, `BUILDER_NOTES.md` (L/a = 1.70, 23.67%, 0.0434 rad, 0.726%, 25.8%), all
  `grading/R*/` folders, all run folders (their ANSWER.md hold the key), this audit's `check_P038.py`, `recompute_P038_results.json`,
  `runs_extracted.json`, `check_P038_selftest.*`, `t2_*` outputs, `T3_KEY_RUBRIC.md.proposed`. Never shipped: `paper.txt`, `pages/`,
  the knowledge unit's EXTRACT.md (paper extraction).

## 5 Fixes

1. [Left out for this release (8 November 2026): the audit's account of what the paper supports for a restated T3, the basis of a new problem that has not been run.]
3. **T2:** label it in the problem record as a variant of p. 364 (key computed by us, now checked by two independent routes). Optional,
   to make the third number discriminate the unequal-arm adaptation: `T2/TASK_PROBLEM.md:42` "3. $\delta(L/2)$ in µm." →
   "3. $\delta$ at $x=5.00$ mm, in µm."; key 31.7700837 µm (this audit, two routes); rubric row `T2/RUBRIC.md:30` changed accordingly.
4. **T1:** none required. Optional (open-book purity): state the inputs also as the paper's groups, a/(δcE′/σb) = 0.302474 and
   Γ0/(σbδc) = 0.0152778, so the parameters are visibly inside the paper's range.
5. **Release:** add an output block (JSON with the field names of `check_P038.py`) to T1 and T2; record the exact header used
   (run `TASK.md`), not the current `setup_probe_runs.sh`.

## 6 Not verified

- The T3 diagnostic arithmetic was rechecked, but no claim in T3 (A, B or C) was checked against a full simulation (none exists).
- The knowledge unit's code and claims were not audited (only its file list and that it contains Eqs. 45–47).
- Run session logs (`probes-logs/*.jsonl`) were not reread for rule compliance; I rely on `PROBE_PILOT_SUMMARY.md:20` for "no web use".
- The paper's Table 1 values and Eqs. 47–54 beyond those cited above were not checked.
- Whether the T2 beam idealization (simple beam theory, no root rotation) is physically adequate is not at issue (the statement
  prescribes it); its error against a 2-D elastic solution was not estimated.
