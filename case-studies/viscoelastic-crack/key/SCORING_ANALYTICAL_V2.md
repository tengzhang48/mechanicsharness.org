# Analytical-weighted assessment, version 2

Recorded on 2026-10-22 after the user clarified “We will rate analytical solution higer”. This is a user-directed supplementary assessment of the existing pilot submissions. It was introduced after launch and after GPT's delivery, so it is not the original frozen rubric or a preregistered outcome. Retain the original rubric and scores alongside this assessment. Do not change the active solvers' instructions or send them feedback.

## Weights

| Category | Maximum | Required evidence |
| --- | ---: | --- |
| Model and conditions | 10 | Correct geometry, moving coordinate, relaxed history, cohesive law and supplied free energy |
| Analytical derivation and steady branch | 60 | Derived stress/zone solution, strain history, loading/separation relation, speed range and limits |
| Scalar numerical verification | 5 | Required scalar outputs meet the original public tolerances |
| Field numerical verification | 5 | All requested field components meet the original public tolerances |
| Spatial dissipation and numerical evidence | 10 | Actual spatial integration, refinement, wake-tail control and energy consistency |
| Reproduction and reporting | 10 | Reproducible code and verification, accurate claims and limitations |

Analytical content accounts for 70 points; numerical verification for 20; reproduction and reporting for 10.

## Analytical anchors

Model and conditions retain the original five 2-point anchors in SCORING.md.

Award up to 15 points for each of the following four derivations:

1. **Stress and zone selection:** an analytical field representation solving the mixed strip boundary problem, its branch and boundary behavior, and removal of the cohesive-front singularity. Explain how the loading depends on zone length.
2. **Material history:** an analytical strain/displacement history satisfying the creep law and relaxed far-ahead condition, with the relevant homogeneous mode excluded.
3. **Separation and loading:** derive the critical-separation equation and the relations selecting zone length and loading at a prescribed speed. Justify the physical steady branch.
4. **Speed range and limits:** derive the finite-endpoint condition where applicable, explain admissibility on either side and the limiting zone length/loading/energy behavior, and establish the slow-speed limit.

For each item give 15 for a justified correct analytical result, 7.5 for substantial correct but incomplete analytical progress, and 0 for a missing, incorrect or exclusively discretized numerical result. Explain the item-level judgment with submitted formulas and file locations. A matching table, a mesh computation or a curve fit does not itself demonstrate an analytical derivation.

Exact integral representations, convergent series, transforms and implicit scalar selection equations are eligible analytical solutions when their derivation, branches and conditions are justified. Numerical quadrature or root finding to evaluate such a derived representation does not turn it into an exclusively numerical solution. Accept different valid analytical methods and notation; do not require the source paper's particular construction.

## Remaining anchors and gates

Use the original public accuracy thresholds without alteration. Scale the original scalar and field category points by one quarter, and spatial-dissipation points by one half. Keep the original reproduction/reporting anchors unchanged. The three bulk values therefore receive 2.5 points in total, demonstrated spatial integration 5, and refinement/tail/residual evidence 2.5. Do not silently replace an actual integral by energy subtraction.

Retain the original four gates and their separate results: numeric accuracy, energy consistency with a real spatial integral, scientific support, and reproduction/access. Add a separate analytical-recovery status: **supported**, **partial**, **not demonstrated**, or **inconclusive**. Supported analytical recovery requires all four core derivations above to be correct and justified. A weighted total cannot establish that status by itself or override a failed original gate.

An external reviewer must assess the analytical items. Mechanical reweighting of an external review is allowed only where its item-level evidence unambiguously establishes these analytical anchors; disclose that the weighting is organizer arithmetic. If the original review does not distinguish an analytical representation from a numerical approximation, request a focused external follow-up using the supplied reference solution. Do not author a substitute scientific verdict locally.
