# Analytical solution is the intended benchmark target

The user clarified on 2026-10-22 at approximately 16:41 UTC: “We are not wanting the numerical solutions, right? We want to check the analytical solution”. This clarification arrived after the solver runs and Grok review had started.

The intended assessment is whether an agent derives an analytical solution of the stated finite-strip crack problem from the model. Derived integral representations, convergent series and implicit scalar selection equations can constitute an analytical solution. Numerical quadrature and root finding may evaluate those relations and check the result. A discretized field solution by itself does not demonstrate analytical recovery.

The assessment should identify the submitted derivations of:

1. The strip stress field, analytic branch, mixed boundary conditions and cancellation of the cohesive-front singularity.
2. The strain/displacement history with the relaxed far-ahead state.
3. The loading and critical-separation relations selecting the cohesive-zone length.
4. The admissible speed range, finite terminal speed and endpoint and slow-speed limits.
5. The local dissipation formula from the supplied free energy, with numerical integration serving as a verification check.

Each finding needs an actual formula and justification in the submission, with an external review against the supplied source solution. Numerical agreement alone is insufficient. A different exact analytical representation is eligible; reproducing the paper's notation or method is not required.

## Existing-run limitation

The issued problem explicitly says “The analytical or numerical method is yours to choose”, and the frozen rubric accepts boundary-integral, finite-element and spectral methods. That scope is broader than the clarified analytical target. Preserve those original inputs and scores. Report analytical recovery separately, naming the evidence and remaining gaps, rather than retrospectively claiming an analytical-only rule was given to the solvers. No feedback or reference material is to be injected into the active solver runs.

The ongoing Grok review already requests a source-based check of the four core derivations. Inspect its finished review for sufficient analytical evidence. If its assessment only establishes numerical agreement, obtain a focused follow-up from Grok using the same supplied solution and this clarification. Do not create an additional blind solver phase.

A future analytical-only benchmark must explicitly require these derivations in its solver-facing prompt before launch. The clarification here is a recorded reporting focus, not a revision of the frozen runs or a new preregistration.

## User-directed weighting

The user subsequently requested a higher rating for the analytical solution. [SCORING_ANALYTICAL_V2.md](private/SCORING_ANALYTICAL_V2.md) records 70 points for analytical content, 20 for numerical verification, and 10 for reproduction/reporting. The exact rubric hash and timing are retained in reviews/grok46_gpt61_2026-10-22/ANALYTICAL_RUBRIC_V2_RECORD.json. This supplementary weighting does not rewrite the original frozen assessment.
