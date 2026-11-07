# Scoring plan

The task's accuracy requirements are public in `agent/CONTRACT.md`. All main physical inputs and scored sampling locations are public. The paper and computed answers remain private. There is no hidden mesh requirement or undisclosed accuracy cutoff.

## Required gates

Report each gate separately:

1. Numeric accuracy: 0.5% relative error in epsilon, a, R and finite terminal speed; absolute error at most 0.005 in each normalized field component; bulk-loss error at most 0.01 times max(1, reference bulk loss).
2. Energy consistency: absolute residual |bulk loss + 1 - R| divided by max(1, R) at most 0.01, and an actual spatial integral with tail and refinement evidence.
3. Scientific support: correct model and conditions, justified steady branch and speed range, and checks that use the reported fields. A reviewer can credit a correct alternative method.
4. Reproduction and access: code regenerates the results within the allowed environment, and the run respected the closed book rules.

An overall study pass requires all four. The JSON scorer implements the numeric gates and energy residual only. It leaves scientific review and study pass unset. It cannot establish that the submitted energy number came from an independent spatial integral; code and convergence evidence must be reviewed.

## Review score

| Category | Points | Evidence required |
| --- | --- | --- |
| Model and conditions | 10 | Correct geometry, moving origin, relaxed preload, cohesive law and energy definition |
| Derivation and steady branch | 20 | Loading relation, zone selection, speed range and limiting behavior derived from the model |
| Scalar results | 20 | Published-case outputs meet public tolerances; partial credit reflects accuracy |
| Stress and strain fields | 20 | Field samples and equilibrium, compatibility and constitutive checks |
| Spatial dissipation | 20 | Actual field integral, refinement and wake truncation, with energy consistency |
| Reproduction and reporting | 10 | Runnable code, clear unresolved quantities and evidence matching claims |

The total cannot override a failed gate. Do not penalize a formulation merely for differing from the reference. Correct conformal mapping, boundary integral, finite element, spectral and other methods are eligible. Do not reward or deduct points for quadrature count, agent verbosity or source recognition alone.

Apply these category anchors consistently:

- Model and conditions: 2 points each for the strip/kinematics, moving-coordinate convention, relaxed far-ahead history, cohesive conditions including separation, and the supplied free energy. Give 1 for an incomplete but correct treatment and 0 for a missing or incorrect one.
- Derivation and steady branch: 5 points each for stress/zone selection, material-history strain construction, loading/separation relation, and speed range with limiting behavior. Give 5 for a justified correct result, 2.5 for a substantial correct but incomplete derivation, and 0 for a missing or incorrect result.
- Scalar results: 5 points for an accurate finite terminal speed. Divide the remaining 15 equally among the nine positive outputs epsilon, a and R at the three speeds. Award each numerical item's share only if it meets the public tolerance. Give an absent or unresolved output 0.
- Fields: divide 20 equally among the 180 requested component values. Missing, duplicate or inaccurate samples receive 0 for the affected values.
- Spatial dissipation: 5 points for accurate bulk-loss numbers, split equally across the three speeds; 10 for a demonstrated spatial integration (5 if implemented but not successfully converged, 0 if only subtraction or unsupported numbers); and 5 for numerical evidence, comprising 2 for ordinate/spatial refinement, 2 for wake-tail control and 1 for the energy residual. A subtraction-only calculation still fails the integration gate regardless of its numeric points.
- Reproduction and reporting: 4 points for the entry point reproducing the JSON, 4 for runnable verification code with matching recorded outputs, and 2 for claims/limitations supported by the work. Use half of an item's points when evidence is substantial but incomplete, and 0 when absent or incorrect.

The scientific-support gate requires the core stress/zone, strain-history, separation/loading and speed-range claims to be correct and supported, as well as actual spatial energy integration. Report disagreements between reviewers by item instead of silently averaging contradictory physical conclusions.

This is a proposed rubric. Before launch, have a domain reviewer apply it to the reference and substantive alternative/incorrect derivations, resolve disagreements, then freeze its version. Log subsequent corrections and apply them to every submission, reporting which scores changed.

## Evaluator checks

`test_scoring.py` checks wire-format behavior using a copy of the numerical key as a synthetic exact payload. It also injects an unbounded-speed claim, a small-zone length approximation, missing relaxed strain history, zero bulk loss, missing data, duplicate cases/samples, booleans and nonfinite JSON. These are evaluator mechanics tests, not independent correct/wrong scientific submissions and not agent evidence.

Published Figure 9 scalar values provide an external source check for the key. The local history ODE, field identities, elastic control and power integration provide additional numerical consistency checks. A complete independently reviewed positive control and a complete solver using an inappropriate small-scale approximation remain desirable before a scientific study. The synthetic mutations do not replace that review.
