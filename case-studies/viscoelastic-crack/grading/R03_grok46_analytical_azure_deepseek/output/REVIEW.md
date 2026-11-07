# Review of P01_deepseek (preliminary)

**Verdict (preliminary):** fail. Numeric accuracy fails public tolerances. Analytical recovery is not demonstrated for the four core closed-form anchors; the submission is a finite-difference mixed-boundary solver with a correct constitutive reduction.

This file will be completed after executing `reference/score_submission.py`, reproducing `python solve.py case.json ...`, and checking source equations. Unresolved items below are marked inconclusive.

**Disclosure.** The user requested analytical weighting after solver runs began. The original frozen task (`agent/PROBLEM.md`) explicitly allowed analytical or numerical methods. Both rubrics are reported. An alternative valid numerical solution remains eligible under the original rubric; mesh agreement cannot establish analytical recovery under `reference/SCORING_ANALYTICAL_V2.md`.

Reviewer: Grok 4.6 (external). Candidate: DeepSeek-V4.1-Flash Azure/OpenCode closed-book submission.

## Four analytical items (preliminary)

1. **Stress and zone selection — 0 / 15 (inconclusive pending source line cites).** Candidate reduces to a harmonic mixed strip BVP and cancels the front singularity by a three-term least-squares SIF fit on a finite-difference grid (`candidate/solve.py` `sif()`, `root()`; `ANSWER.md` §§1,5). No closed-form stress field analogous to Hui–Long eq. (12)/(17) is derived.
2. **Material history — 15 / 15 (preliminary).** Overstress convolution from the SLS ODE with decay ahead (`DERIVATION.md`; `solve.py` `build_avg_matrix`). Eligible integral representation.
3. **Separation and loading — 0 / 15 (preliminary).** Critical separation and K=0 are imposed on the discrete system; no analytical selection equation (paper eq. 51 / a(ε) relation).
4. **Speed range and limits — 7.5 / 15 (preliminary).** Slow-speed `a→√δ̄_c` is argued analytically; finite `V_crit` is a discrete SIF interpolation (`ANSWER.md` §3; hardcoded `13.44` in `solve.py`).

## Category scores (preliminary; incomplete arithmetic)

Original rubric (max 10/20/20/20/20/10): inconclusive pending executed scorer and field-by-field tally.
Analytical v2 (max 10/60/5/5/10/10): inconclusive pending same.

## Gates (preliminary)

- Numeric accuracy: fail (submitted report: 8 failed checks, including two ε values and six field components). To be re-executed.
- Energy consistency: likely pass (spatial integral present; residuals in submitted report pass). To be re-executed.
- Scientific support: inconclusive.
- Reproduction/access: inconclusive pending live `solve.py` run. Closed-book claim supported only by `evidence/RUN.json` / transport metadata (network none; model relay). Unlogged actions are not proven absent.

## Defects / limitations (preliminary)

- `critical_speed` is a literal `13.44`, not recomputed by the entry point.
- Tabulated `_GUESS` supplies reported `epsilon` values to printed precision.
- Organizer key is a Codex reproduction of Hui–Long, not independent ground truth.

Completed assessment will replace this stub.
