# X — Energy release for periodic cracks in an elastic film on a foundation (Xia–Hutchinson 2000, §4.2)

Verdict: READY AFTER FIXES

Source paper: Z. C. Xia and J. W. Hutchinson (2000), "Crack patterns in thin films", J. Mech. Phys. Solids 48, 1107–1131,
Eqs. 8-12, 28, 30 (copy: `reference/frozen/paper.pdf`, `paper.txt`).

Files (`<pilot-folder>/`):
- Statement: `frozen_inputs/TASK.md` (sha e88f5297…) + `frozen_inputs/PROBLEM.md` (sha ff71c35e…); identical in S01-S04
  and D01 (sha256 check). Controller phase messages during the 50-minute budget (TASK.md:7). No knowledge entry.
- Key: `reference/frozen/REFERENCE.md`, `reference.json` (authored by GLM-5.3, checked by Grok-4.6, KEY_FREEZE.json).
- Grading: `RUBRIC.md` (7 items, Grok reviewer), numeric agreement at 1e-6 absolute (KEY_FREEZE.json "tolerances") recorded
  in `NUMERICAL_EVIDENCE.json/.csv`; no checks.py.
- Runs: `submissions/S01` (GPT-6.1), `S02` (DeepSeek), `S03` (Kimi, nothing delivered), `S04` (Kimi subscription),
  `followup_runs/D01_…`, `continuations/K01_…`; corpus copies under `<solver-runs>/`.
- This audit: `audit_xh.py`, `audit_xh.out`, `audit_xh_output.json`.

## 1 Statement

Complete and faithful to the paper: plane-stress film with prestress σ₀I on a local foundation f = −ku (PROBLEM.md:8-37
= XH2000 model, Eqs. 4-7), the energy functional (l.39-51, which also resolves the paper's ambiguous Eq. 6, REFERENCE.md:26-27),
G per advancing crack area and G = (K_I² + K_II²)/E (l.54-58), the two geometries (l.62-76 = XH2000 Fig. 3 and §4.2:
"advancing midway between previously formed infinite cracks … spaced a distance 2H apart"), ν = 0.3 and the four sampling
locations (l.11, 90-92). ℓ is defined by the problem ("decay length of the wake field", l.80), so the normalization is
unambiguous. No invented data.

Issues:
- **statement — coordinator text shipped to solvers.** PROBLEM.md:3-4: "Draft for Teng's review. This is the proposed
  source-faithful reproduction task; source solutions and the proposed geometry extension are excluded." It tells the
  solver the task reproduces a published source (a recall cue) and is irrelevant to the problem.
- **statement — accuracy not stated.** Nothing in TASK.md or PROBLEM.md says how accurate the sixteen values must be, yet
  the numeric record uses 1e-6 absolute (KEY_FREEZE.json "tolerances"). The method is open (PROBLEM.md and TASK.md:3), and
  S02's 2D FEM values (0.2-1.6 % off) were counted 0/16 (NUMERICAL_EVIDENCE.csv S02 rows; R01 grade.json item 5).
- **minor — recognizable source.** The closed forms are printed in the paper (Eqs. 28, 30); the paper is freely hosted
  (README.md: harvard.edu link). Recall cannot be controlled.

## 2 Key

| Item | Key | Origin | Recomputed here |
|---|---|---|---|
| ℓ = √(Ēh/k) | REFERENCE.md:45 | XH2000 Eq. 12 (with Eqs. 8-9) | yes (strip ODE Ēu″ = (k/h)u) |
| γ_A = tanh(H/2ℓ) | REFERENCE.md:90-93 | XH2000 Eq. 28 | **yes, two routes**: by hand (integral of the functional over one period of the exact wake) and numerically (solve_bvp on the half-strip, quadrature of the stated functional, no tanh): max difference 2.2e-16 |
| γ_B = 2 tanh(H/2ℓ) − tanh(H/ℓ) | REFERENCE.md:105-112 | XH2000 Eq. 30 and G = 2G′(H) − G′(2H) (p.1120) | yes, same numerical route |
| k_I = √((1−ν²)γ), K_II = 0 | REFERENCE.md:119-123 | XH2000 Eq. 28 second form | yes |
| 16 sampled values | REFERENCE.md:145-150, reference.json | closed forms at ℓ/H = 0.5…2 | equal to my closed-form values exactly (max difference 0) |
| limits H/ℓ → 0, ∞ | REFERENCE.md:135-143 | series of the closed forms; XH2000 §4.2 ("H exceeds approximately 3l") | by inspection |

The key author is a model (GLM-5.3) with a model checker (Grok-4.6) (KEY_FREEZE.json); my two-route recomputation and the
paper's printed Eqs. 28 and 30 now back it independently.

## 3 Checker

No checks.py. The automatic numeric layer (16 values, 1e-6 absolute) was re-implemented in `audit_xh.py`:
(a) key: 16/16. (b) saved runs: S01 16, S02 0, S04 16 — equal to NUMERICAL_EVIDENCE.json; D01 16 and K01 16 (not in
that file; README.md: "matched the sixteen requested values"). (c) mutants:

| Mutant | 1e-6 abs | 0.5 % rel |
|---|---|---|
| k_I with Ē (no 1−ν²) | 8/16 (k_I fail) | 8/16 |
| omit foundation energy (reference.json control C1) | 8/16 | 8/16 |
| Case B per period instead of per tip (×2) | 8/16 | 8/16 |
| Case B spacings swapped (control C3) | 8/16 | 8/16 |
| ℓ/H read as H/ℓ | 4/16 (only ℓ/H = 1) | 4/16 |
| decay length 2ℓ | 0/16 | 0/16 |
| +0.5 % uniform bias (a converged-FEM-type error) | 0/16 | 16/16 |

(d) Band: 1e-6 is right for the closed form but not for an open-method task with no stated accuracy; every physically
wrong mutant above already fails at 0.5 %, so a band of ~0.1-0.5 % relative separates right from wrong while admitting
numerical methods. (e) Not covered automatically: the fields (Q1), the limits and σ_c (Q3), checks (Q5) — rubric only.
Rubric (RUBRIC.md:15-38): items 1-5 are decidable against REFERENCE.md; 6-7 are judgments. Item 5 does not state the
numeric band, so "right/partly" for numbers depends on the reviewer (R01 used 0/16 at 1e-6 and 3/16 at 1e-3).

Issues:
- **grade — tolerance not stated and too tight for the open method.** KEY_FREEZE.json fixes 1e-6 absolute after the
  statement left accuracy open; an approximately right numerical answer counts 0/16.

## 4 For outside solvers

- Received: TASK.md + PROBLEM.md only; 50 min including reporting, with controller phase prompts at minutes 13, 15, 45, 50
  (TASK.md:7); local Python/SymPy/SciPy; no web.
- Output: results.json with `rows` of ell_over_H, gamma_A, gamma_B, kI_A, kI_B and a `checks` list (TASK.md:9).
- Grading: numeric comparison (1e-6) + one blind Grok reviewer with RUBRIC.md and REVIEW_TASK_TEMPLATE.md.
- Machine paths: runner scripts (`runner.py`, `subscription_runner.py`) and README.md list `<volumes>/...` roots;
  solvers need none.
- Public: frozen_inputs/TASK.md, PROBLEM.md (after fixes), RUBRIC.md (holds no values), REVIEW_TASK_TEMPLATE.md.
- Hidden: reference/ (all versions), KEY_FREEZE.json, NUMERICAL_EVIDENCE.*, submissions/, reviews/, followup_runs/,
  continuations/. Never distribute `reference/frozen/paper.pdf` or `reference/huang_etal_2003/*.pdf`.

## 5 Fixes

1. PROBLEM.md:3-4 delete "Draft for Teng's review. This is the proposed source-faithful reproduction task; source
   solutions and the proposed geometry extension are excluded." Reason: coordinator text and a recall cue.
2. PROBLEM.md:90-92 append: "Report each value to a relative accuracy of 0.1 % or better." (or another stated band) and use
   the same band in the checker; keep the 1e-6 comparison only as a note. Reason: §1, §3.
3. Add `checks.py` that reads results.json `rows`, compares the 16 values with the stated band, and fails missing values;
   `audit_xh.py:compare()` is a template.

## 6 Not verified

- Whether S02's FEM values converge to the key under refinement (its errors 0.2-1.6 % are consistent with discretization,
  but I did not rerun its code).
- The Grok grades (R01-R05) and the K01/D01 follow-up reports beyond their numbers.
