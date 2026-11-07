# E — Grain-boundary diffusion wedge in a thin film (multiphysics round, problem E)

Verdict: READY AFTER FIXES (key values all reproduced to 10⁻⁵ or better; part (c) must be defined; the round's
"key checks" count must include λ₁; part (a)'s rubric must say whether the amplitude c₁ is required)

Source paper: H. Gao, L. Zhang, W. D. Nix, C. V. Thompson & E. Arzt (1999), "Crack-like grain-boundary diffusion wedges
in thin metal films", Acta Mater. 47(10), 2865–2878; Eqs. 1–21 (pp. 2868–2869), Table 1 (p. 2869), Appendices A–B
(`<source paper PDF, not included>`; page images
`work/pg-04.png` = p. 2868, `work/pg-05.png` = p. 2869).

Files: `multiphysics_round/TASK_HEADER.md` (budget 170), `E/PROBLEM.md`, `E/KEY.md`, `E/checks.py`, `reproduce_keys.py`
(E block), `summary_data.py`; bundle `assets/solution-evidence-2026-10-25/{grain-boundary_check.py,
check-outputs/grain-boundary_check.txt}`; runs MP-02, 03, 05, 07, 09, 15, 18, 20, 32, 40, 44, 47. Audit here:
`key_recompute.py` (+ `.out`), `audit_checks.py` (+ `.out`), `mutants_result.json`, `test_proposed.py` (+ `.out`),
drafts `PROBLEM.md.proposed`, `checks.py.proposed`.

## 1 Statement

- Model matches the paper line by line (image-checked): S(z,ζ) = Eq. 1, traction = Eq. 2, μ = Eq. 3, flux = Eq. 4, mass
  balance = Eq. 5, BCs = Eqs. 8–10, t₀ = Eq. 12. I re-derived ∂σ/∂τ = ∫S σ‴ dξ with this t₀, so λf = −∫S f‴ (Eq. 13)
  follows and the stated sign is right (MP-18's "unphysical sign" and "π normalization" challenges are wrong). No invented
  narrative. The uniform initial traction is incompatible with σ(0,t) = 0 at t = 0 (MP-02 challenge); that is the paper's
  formulation and the modal solution handles it.
- **statement — (c) "On what time scale … does the traction relax?" (PROBLEM.md:31) has no defined measure.** MP-03 (c)
  graded contested ("e-folding and 95% threshold times differ"); the case page calls it "a statement ambiguity". The band
  0.05–0.15 passes the e-folding time 1/λ₁ = 0.0852 and the half-time at z = h 0.0953, and fails the 95 % time at z = h,
  0.291 (`key_recompute.out`), all legitimate readings.
- **statement — (a) "Describe how σ_gb(z,t) evolves" (PROBLEM.md:25) has no stated content.** Graders added the late-time
  amplitude c₁ (not asked, not in the key) and counted its omission or misvalue as decisive: MP-02 (a) partly
  ("missing first-mode amplitude", GRADE.md:27), MP-05 D2 (GRADE.md:30); MP-09's c₁ = 1.547 only a side error
  (GRADE.md:33). Either ask for c₁ or say it is not required.
- minor — (b) K∞ is a handbook value (edge crack, 1.1215); the header allows "your own knowledge", so a recalled value
  is a legitimate answer, but it then tests recall, not the model (see §3).
- Cueing: t₀ = 4πkTh³/(E*δ_gbD_gbΩ) and the kernel are the paper's own; recall possible (KEY.md:87).

## 2 Key

| Item | Key (checks.py:4, KEY.md) | Origin | Recomputed (`key_recompute.out`, paper's collocation, N = 30…480, Richardson) |
|---|---|---|---|
| λ₁ | 11.74265 | Table 1 | 11.742668 (runs MP-02/05/09 agree: 11.742668) |
| λ₂ | 264.6348 | Table 1 | 264.63424 |
| f₁(0.5)/f₁(1) | 0.8020 | ours (KEY.md:42) | 0.801957 |
| K∞/(σ₀√(πh)) | band 1.09–1.15; handbook 1.1215 (KEY.md:50–55) | handbook, "not printed in the paper" | **1.121522 from the model's own kernel** (N = 480: 1.12127; extrapolated 1.121520–1.121522). The key's "1.117, 1.119, 1.120 … converges towards it" is now confirmed. |
| relaxation time | 0.05–0.15, "t₀/λ₁ ≈ 0.085" | p. 2869 "of order t₀/λ₁" | 1/λ₁ = 0.08516 |
| (not in key) c₁ | — | graders: 1.531049 | projection 1.52887 (N = 120), 1.52996 (240), 1.53050 (480); extrapolated **1.531048** |
| (not in key) t_½ at z = h | — | — | 0.09527 (full modal sum, N = 480) = ln(2c₁)/λ₁ = 0.09530 |

- **Lead on c₁ refuted.** "Three graders give 1.531 where a projection gives 1.529 and is still rising": the projection
  rises because it converges as ~1/N; at N = 480 it is 1.5305 and both Richardson pairs give 1.53105. The graders' value
  is right. The case page's "our projection 1.529 at N = 120" should give the converged value.
- Paper text vs Table 1 (decay ratios "25 … 1000" vs 22.5 … 709): recorded in KEY.md:35–37, confirmed on p. 2869.
- KEY.md:43 "agrees with Fig. 4a and Fig. B1" is a figure reading; the value itself is computed.

## 3 Checker

(a) Key → `checks.py` copy: class `exact`, all `*_ok` true. (b) Rerun on the 9 runs with results: identical to the saved
`auto/*.json` (MP-07 no output and no auto file; MP-15, MP-32 no results.json). (c) Mutants (`mutants_result.json`):
- **grade — λ₁ is left out of the round's "key checks" (lead confirmed).** `summary_data.py:22` counts only keys ending in
  `_ok`; E's λ₁ result is `lambda1_class`. Mutant "λ₁ × 4π, everything else right" → class `wrong` but key checks **4/4**.
  In the saved runs MP-40 (λ₁ 1.1 % off, `approximate`) shows 3/4, and MP-44's single pass is a typed-in K.
- **minor — K check passes recall (lead confirmed).** MP-18 (λ₁ = 188) and MP-44 (λ₁ = 88) pass `K_inf_ok` with a typed
  1.1215 (also MP-40, MP-47). The value is right, so this is not a wrong pass; it shows K∞ cannot evidence the
  computation.
- minor — `mode1_shape_ok` band 2 % passes the paper's Appendix B crack-shape estimate 0.8111 (1.1 % off); computed runs
  give 0.80196 (MP-40's under-resolved 0.8001 is 0.24 % off).
- (c)-band: passes any λ₁ in 6.7–20 via 1/λ₁; fails the 95 % reading (0.291) — consequence of the undefined measure.
- Fail ✔: MP-20 (self-consistent but wrong, own-checks pass), Griffith K = 1.0, K in σ₀√h units, K at finite time
  0.71, λ₂ = λ₃, t = 1/λ₂.
(d) Bands vs precision: λ₁ 1 %, λ₂ 2 % are wide against the key (10⁻⁵) and still separate MP-40's 1.1 %/2.7 % errors;
K band 1.09–1.15 admits N = 30 collocation (1.1175) and excludes 1.0 — fine.
(e) Not covered: `method`; (a)'s description; the "what does the boundary become" part of (b).

## 4 For outside solvers

- Received: `TASK_HEADER.md` (170 min) + `E/PROBLEM.md`, byte-identical for all 12 runs (`audit_checks.out` §4). Three
  runs ended without results for harness reasons (MP-07 Cline compaction stop, MP-15 text-only turn, MP-32 provider
  quota; case page Appendix B, HARNESS_NOTES.md correction 2) — an outside harness needs a long-run policy.
- Tools: numpy/scipy suffice (the paper's Gauss–Chebyshev collocation reproduces Table 1 in < 1 s at N = 120).
- `E/checks.py` has no absolute paths; `grade_v2.py` gives E a 1800 s rerun timeout (line 49).
- Public: `TASK_HEADER.md`, `E/PROBLEM.md` (proposed), `GRADER_INSTRUCTIONS.md`.
  Hidden: `E/KEY.md`, `E/checks.py` (key values line 4), `reproduce_keys.py`, the runs and grades, the bundle's
  `cases/Grain_Boundary_Diffusion_Case_Study_2026-10-25.html` and `assets/solution-evidence-2026-10-25/grain-boundary_check.*`.

## 5 Fixes

1. `E/PROBLEM.md:31` "(c) On what time scale, in units of t₀, does the traction relax?" → "(c) Give the time t_½, in
   units of t₀, at which the traction at z = h has fallen to σ₀/2." Key 0.09527 (ours, `key_recompute.out`), band 1 %.
   (Minimal alternative: ask for the e-folding time of the late decay; then (c) only repeats 1/λ₁.)
2. `E/PROBLEM.md:25–26` add "At late times σ_gb ≈ c₁σ₀f₁(z/h)e^(−λ₁t/t₀) with f₁(1) = 1; give c₁" and a report key
   `mode1_amplitude_c1`; key 1.5310 (ours; graders 1.531049). Or, if c₁ is not wanted, `E/KEY.md` must say "(a) does not
   require the amplitude" so graders stop scoring it. (Draft of 1–2: `PROBLEM.md.proposed`.)
3. `E/checks.py`: add `lambda1_ok` (exact band) so `summary_data.py:22` counts it — or change `summary_data.py:22` to count
   `lambda1_class == "exact"`; shape band 2 % → 0.5 %; time check → t_½ ± 1 %; add `c1_ok` ± 1 %; own check
   t_½ = ln(2c₁)/λ₁ (draft `checks.py.proposed`; `test_proposed.out`: key passes; c₁ = 1, 1.227, 1.547, the B6 shape and
   the old 1/λ₁ answer fail).
4. `E/KEY.md:50–55` replace "handbook value, not from this paper" support with the model's own limit 1.121522 (ours,
   N ≤ 480 + Richardson); note that a recalled 1.1215 passes.

## 6 Not verified

- c₁ and t_½ come from one discretization family (the paper's collocation) plus extrapolation; agreement with three
  graders' independent computations (1.531049, 1.53072–1.53105) is the only cross-check.
- Eigenvalues λ₃–λ₅ and Appendix A's "extending method" were not recomputed.
- Blind grades were not regraded; MP-47 has no grade.
