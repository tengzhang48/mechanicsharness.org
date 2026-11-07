# B — Suddenly pressurized crack in a poroelastic solid (multiphysics round, problem B)

Verdict: READY AFTER FIXES (three statement additions; KEY.md body must be rewritten to the corrected values; two bands
widened). The checker already uses the corrected key and does separate the drained-tip value from the crack-scale value.

Source paper: E. Detournay & A. H.-D. Cheng (1991), "Plane strain analysis of a stationary hydraulic fracture in a
poroelastic medium", Int. J. Solids Struct. 27(13), 1645–1662, Eqs. 1–2, 8, 13–39, p. 1655 parameters
(`<documents>/crack_corpus_pdfs/Plane strain analysis …pdf`; page images `work/pg-05..08.png`). Cross-check: Yu, Bouklas,
Landis & Huang (2020) J. Appl. Mech. 87, 031005, Eq. 12 (`…/Poroelastic Effects on the Time- and Rate-Dependent Fracture
of Polymer Gels.pdf`, text `work/yu2020.txt`:327–331).

Files: `multiphysics_round/TASK_HEADER.md`, `B/PROBLEM.md`, `B/KEY.md`, `B/checks.py`, `B/FE_CHECK.md`,
`B/CASE_STUDY_GPT_CHALLENGE.md`, `B/reproduce_yu2020.{py,txt}`, `reproduce_keys.py` (B block); bundle
`assets/solution-evidence-2026-10-25/{poroelastic/wilson_rerun_2026-10-25.txt, check-outputs/poroelastic_closed_forms.txt,
check-outputs/solution_formula_checks.txt}`; runs MP-06, 10, 13, 21, 22, 24, 27, 30, 33, 37, 38, 39, 41, 46, 50. Audit
here: `audit_checks.py` (+ `.out`; `test_proposed.out` = same script on `checks.py.proposed`), `mutants_result.json`,
drafts `PROBLEM.md.proposed`, `checks.py.proposed`.

## 1 Statement

- Model checked against the paper and by derivation: the plane-strain constitutive law (PROBLEM.md:11–12) follows from the
  3D Biot law with e_zz = 0 (σ_zz eliminated; I re-derived 2Gζ = α(1−2ν)S + α²(1−2ν)²p/(ν_u−ν)); undrained ζ = 0 gives
  the paper's Eq. 16. Loadings = paper Eqs. 9–10; parameters = p. 1655 ("ν = 0.2, ν_u = 0.4 and α = 0.89",
  `work/paper_text/detournay_cheng1991.txt`:690). No invented narrative.
- **statement — order of limits for K(0⁺) (PROBLEM.md:26, 29).** "K … at the crack tips … just after loading (t = 0⁺)"
  admits two readings: tip limit first (0.807593, the key) or undrained crack first (1). The key itself says so
  (KEY.md:24–28); ROUND_RECORD.md:133–134 and the case page ("The next version of the problem states that the tip limit
  is taken first") promise the fix, but `B/PROBLEM.md` is unchanged (and the 24 Oct Kimi K3.1 trials received this same
  unrevised text: `runs/kimi-k31-2026-10-24/workspace/PROBLEM.md` identical). 11 of 15 runs gave the crack-scale value (10
  as 1.0, MP-13 as √π in √L units). The literal text leans to the tip reading, so these remain mostly physics misses, but
  the gap must be closed. `gel/G2/PROBLEM.md` already has the right wording.
- **statement — face shear not prescribed (PROBLEM.md:23–24).** Only σ_n and p are given (GPT MP-37 "Questions"; the FE
  check and all factorizations assumed σ_xy = 0, FE_CHECK.md "Limits").
- **statement — `p_faces_0` (PROBLEM.md:40–41).** "immediately next to the crack faces at t = 0⁺" conflicts with the
  boundary value p = 0; the intended value is the paper's Eq. 17 (the paper itself notes it "is not in equilibrium with the
  imposed boundary condition p = 0", p. 1649). MP-37 flagged the literal trace 0; graders marked (c) contested for it.
- minor — cueing/recall: the statement copies the paper's notation and example values, so a recalled answer gives the
  paper's 0.75 and 1.354, both wrong under the corrected key. That works as a recall trap; record it.

## 2 Key

| Item | Key (checks.py:6–7) | Origin | Status here |
|---|---|---|---|
| K1_0 | 0.807593 | corrected key (KEY.md:16–19): Wiener–Hopf factorization (Codex, Grok, GPT ×3) | not re-derived; supported by FE 0.798–0.801 (FE_CHECK.md table, independent of the kernel), the 25 Oct Wilson rerun (0.80759300, two scales), and Yu 2020 Eq. 12 at ν_u = ½ (factorization reproduces [1/(2(1−ν))]^0.735 within 0.2 %, `reproduce_yu2020.txt`) |
| K1_inf, V1_0, V1_inf | 1, 0.6π, 0.8π | Eqs. 26, 15, 19 (image-checked p. 1649–1650) | recomputed |
| V1_short_coeff | 2.440220 | corrected: 1.354055 faces (Eq. 18) + 1.086165 tips (factorizations) | faces term recomputed; **tip term not verified**: FE only extrapolates to 2.2–2.4 (FE_CHECK.md "Reading") |
| p_faces_0 | 0.749064 | Eq. 17 | recomputed |
| K2_0, K2_inf, V2_inf | 0, −0.33375, −0.838805 | Eqs. 38, 39, 33 (image p. 1651–1652) | recomputed |

- **statement-level defect in the key document (KEY.md:52–64, 70–79, 91, 95).** The correction block (lines 11–34) sets
  0.807593 and 2.440220, but the "Reference answer" still states K(0⁺) = 0.75√(πL) as "the decisive point" and gives Eq. 18
  as the small-time result; the (c) table keeps 0.75 and 1.354 in the value column. Graders coped (MP-37 GRADE.md:5 "The
  corrected key … controls grading"), but a released key must say one thing.
- minor — KEY.md:61 writes Eq. 27's prefactor as √(πG)/[√8(1−ν)]; dimensionally it must be G√π/[√8(1−ν)] (plane-strain
  opening 4(1−ν)K√(r/2π)/G). The print (`work/pg-06.png`) is ambiguous about the radical's extent. Found by the MP-37 grader
  (GRADE.md:27); confirmed here.
- minor — KEY.md:21–22 attributes the form [(1−ν_u)/(1−ν)]^0.743 to Yu et al.; Yu's Eq. 12 is for ν_u = ½ with α ≈ 0.735
  ("immersed"); 0.743 is this study's factorization (also noted by the Codex page review).
- minor — at the paper's values η = α(1−2ν)/[2(1−ν)] = 0.33375 and η_t = (ν_u−ν)/(1−ν_u) = 0.33333 differ by 0.12 %, so
  the 1 % checks on K2_inf, V2_inf cannot tell the two formulas apart (mutant "wrong eta" passes). Only the `formulas`
  string (grader) can.

## 3 Checker

(a) Key → `checks.py` copy: all nine `*_ok` and both own-checks pass; labels `tip_limit_key`, `faces_plus_tips_key`.
(b) Rerun on 15 saved runs: identical to saved `auto/*.json` for all 15. The saved flags were regenerated at 18:13 on
8 Oct, the same minute `B/checks.py` and `KEY.md` were last modified, and every B grading folder holds the current KEY.md
(`cmp` identical), so **the corrected key is what both layers used**.
(c) Lead confirmed: the checker separates the drained-tip value from the crack-scale value. K1_0 = 1 → fails,
`K1_0_matches = outer_limit_1`; 0.75 → fails, `paper_heuristic_0.75`; √0.75 = 0.866 → fails; 0.808, 0.807 (GPT) → pass.
Other mutants (`mutants_result.json`): MP-30 doubled volumes, ν↔ν_u swap, p sign, K2_inf = 0 (seen ×7), −α (MP-46),
+η, V2 with ν_u, faces-only 1.354 (MP-39), 2/√π (seen ×5) all fail ✔. Pass although wrong: "η_t for η" (above).
(d) Bands: 1 % on K1_0 fails a resolved FE answer (0.798 fails; FE gave 0.7978–0.8013) — **grade** (right method,
rejected). 1 % on V1_short_coeff fails any FE route (2.2–2.4) and is tighter than the evidence for the tip term. K2_0 band
0.01 absolute is fine (exact 0).
(e) Uncovered: `formulas` strings; whether K1_0 is reasoned with a drained tip (grader layer). The own-checks are weak by
design (V ratio passes MP-30's doubled volumes; Griffith check just restates K1_inf).

## 4 For outside solvers

- Received: `TASK_HEADER.md` (50 min) + `B/PROBLEM.md`, byte-identical for all 15 runs (`audit_checks.out` §4). The
  24 Oct Kimi K3.1 and Opus trials used a different prompt (`runs/kimi-k31-2026-10-24/prompt.txt`) with the same
  PROBLEM.md — release one header.
- Schema: ten keys (PROBLEM.md:36–43) + `reproduce_command`. Tools: python3, numpy, scipy (the corrected K1_0 needs a
  Wiener–Hopf factorization or a resolved FE; both are feasible with scipy but long for 50 min).
- `B/checks.py` has no absolute paths; `reproduce_yu2020.py` imports MP-37's code by absolute path; `fe_check_K0.py` needs
  `<home>/miniconda3/envs/fenicsx-env`.
- Public: `TASK_HEADER.md`, `B/PROBLEM.md` (proposed), `GRADER_INSTRUCTIONS.md`.
  Hidden: `B/KEY.md`, `B/checks.py`, `reproduce_keys.py`, `B/CASE_STUDY_GPT_CHALLENGE.md`, `B/DIMENSIONAL_ANALYSIS.md`,
  `B/WHY_DRAINED_TIP_MISSED.md`, `B/FE_CHECK.md`, `B/fe_check_*`, `B/reproduce_yu2020.*`, `corpus-runs/…/review_B_K10/`,
  the 15 runs and grades, the bundle's `cases/Poroelastic_Crack_Case_Study_2026-10-23.html`, `runs/kimi-k31-*/reference/`,
  and `gel/G2/` (same model, states the order of limits and its key).

## 5 Fixes

1. `B/PROBLEM.md:23–24` add "σ_xy = 0" to both loadings (reason: unprescribed shear; all references assumed it).
2. `B/PROBLEM.md:26–27` "Let K be the mode I stress intensity factor at the crack tips" → "For t > 0 let K(t) =
   lim_{r→0} √(2πr) σ_yy(x = L + r, y = 0, t) … K(0⁺) = lim_{t→0⁺} K(t): the limit r → 0 is taken first" (wording of
   `gel/G2/PROBLEM.md`).
3. `B/PROBLEM.md:40–41` define `p_faces_0` = lim_{y→0⁺} lim_{t→0⁺} p(x, y, t), |x| < L, and say the face value is 0.
   (Draft of 1–3: `PROBLEM.md.proposed`.)
4. `B/KEY.md:46–99` rewrite the reference answer to the corrected values (tip-first K1_0 = n₁₁/(1+η_t) = 0.807593;
   A = 1.354055 + 1.086165) and move the paper's 0.75 and Eq. 18 to a "paper value, superseded" note; `:61` √(πG) → G√π;
   `:21–22` say the 0.743 exponent is ours, Yu's Eq. 12 is ν_u = ½, α ≈ 0.735.
5. `B/checks.py:16` K1_0 band 1 % → 2 %, V1_short_coeff 1 % → 5 % (still excludes 0.75, 0.866, 1, 1.354, 2/√π; admits
   FE 0.798); `:19, :22` widen the label windows to match. Draft `checks.py.proposed`; on the 15 runs it changes no flag
   (`test_proposed.out`).
6. Before release, an independent check of the tip volume term 1.086165 (only AI factorizations share it; FE is ±10 %).

## 6 Not verified

- K1_0 = 0.807593 and the tip term 1.086165 were not re-derived here (no independent Wiener–Hopf solution); I rely on
  the cited FE table, the archived factorization rerun and the Yu 2020 comparison. Atkinson & Craster (1991) was not read.
- The FE logs were not rerun (fenicsx env, long).
- Blind grades were not regraded.
