# What the scores miss in the crack tip pilot

The detailed solutions separate several abilities that one total compresses: reducing the material model, closing the mixed crack boundary problem, selecting the physical branch, resolving numerical errors, and checking a solution independently. GPT recovered an explicit analytical construction. Both DeepSeek attempts contain useful material reductions, while their treatment of zone selection and admissibility needs separate attention.

GPT has a completed external Grok review and process audit. Azure DeepSeek’s final analytical assessment is pending; direct DeepSeek has numerical/reproduction results and no final external analytical grade. The observations below assign no new scores. Source arithmetic and qualitative interpretation are identified separately from externally graded findings.

## Useful analytical reductions survive incomplete solutions

Azure DeepSeek derives the SLS overstress equation, excludes the growing far-ahead homogeneous mode, reduces the displacement problem to Laplace’s equation, and writes a history convolution. Direct DeepSeek also derives a harmonic bulk problem and an interface kernel. These are reusable steps even though neither submitted a complete analytical closure of the mixed boundary problem. A failed numerical gate does not erase that partial mathematical structure. See the [Azure derivation](../runs/P01_deepseek/workspace/DERIVATION.md) and [direct solution sections 1 and 2](../runs/P04_deepseek_direct_1017/workspace/ANSWER.md).

The distinction to examine is how a solver moves from this reduction to the stress field, the load–zone relation, and the separation equation. GPT makes those relations explicit and derives the physical branch and endpoint; numerical quadrature evaluates the resulting representation. Its [solution sections 1 to 3](../runs/P03_gpt61_sol/workspace/ANSWER.md) illustrate why analytical recovery can include numerical integration. The [external review](../runs/R02_grok46_analytical_gpt61/workspace/REVIEW.md) supports the four core derivations; its score does not certify every optional asymptotic term or a uniform error theorem.

The original prompt allowed numerical methods. The user introduced the 70-point analytical weighting after launch. These attempts therefore show what each agent produced under that original choice, rather than establishing how the agents would perform under an explicit analytical-only assignment. A later benchmark should state the required analytical deliverables before the run. See the [supplementary rubric](../private/SCORING_ANALYTICAL_V2.md).

## One admissibility claim produces many failed checks

Direct DeepSeek concludes that a pole reaching an analyticity-strip edge implies a terminal speed of 1/lambda = 10. It then marks V=10 inadmissible and omits that case’s fields. Of its 31 failed checks, 22 concern the V=10 case, one concerns the terminal speed itself, and eight concern the remaining cases. Those 22 failures share the submitted admissibility decision; they are not 22 independent physical mistakes. The missing case also removes 60 requested component comparisons, so raw check counts have different denominators across submissions. See [the claim in section 2](../runs/P04_deepseek_direct_1017/workspace/ANSWER.md), [the production entry point](../runs/P04_deepseek_direct_1017/workspace/solve.py) and [the recorded report](../runs/P04_deepseek_direct_1017/SUBMITTED_NUMERIC_REPORT.json).

The published endpoint condition couples the modulus ratio to the critical opening, through (1-d) = (1-lambda) k(V). For the supplied parameters it gives V approximately 13.42357, with a finite steady solution at V=10. The consequential error to investigate is the inference from a particular factorization strip to physical nonexistence. This is more informative than the failure count. The source is main-paper equation 59, PDF page 13; see [the source map and endpoint relation](../private/REFERENCE.md). This comparison is organizer interpretation of the submitted claim against the cited relation, not an external analytical grade for the direct run.

## An absolute tolerance hides substantial relative dissipation errors

At V=0.001 the reference bulk-loss ratio is 0.00729790. Azure DeepSeek reports 0.00892466, about 22.3 percent higher; direct DeepSeek reports 0.01084624, about 48.6 percent higher. Both pass the published bulk-number tolerance of 0.01 max(1,D). Thus a passing numerical item can coexist with a large error relative to a small physical loss.

The numerical inequalities are even loose enough that a hypothetical zero bulk loss at this speed would pass both the bulk-number tolerance and the normalized energy residual, if the other scalars were exact. That would not establish a full study pass: the spatial-integral and scientific-support gates still apply. The useful addition is to report relative loss error alongside the frozen absolute pass/fail result, while retaining the existing grades. The arithmetic is recorded in [the evidence calculations](QUALITATIVE_EVIDENCE_CALCULATIONS.json).

## An energy condition can help select a state without independently validating it

Direct DeepSeek states that its reported low-speed states are a compromise selected partly to satisfy global energy closure on the production grid. It reports a different critical-separation/bounded-stress root and an energy residual that moves outside one percent after refinement at V=0.1. The independently evaluated spatial integral is valuable; its agreement with energy balance is weaker evidence once the same balance helped choose the state. See [sections 4 and 7 of the submitted answer](../runs/P04_deepseek_direct_1017/workspace/ANSWER.md).

Using energy balance as a solving equation can be mathematically valid when equivalence and convergence are justified. The issue here is whether all defining conditions converge together. A useful follow-up should freeze the state obtained from separation and boundedness, then evaluate energy closure on an independently refined calculation.

## Small discrete residuals do not settle every continuum condition

A fresh offline replay of Azure DeepSeek’s checks.py exits successfully in 35.34 seconds and reproduces the reported small Laplace and boundary residuals. At V=0.001 it also gives relative constitutive residuals of about 0.339 for S1 and 0.197 for S2, while the reported relative Laplace residual is 2.15e-8. These residuals have different definitions and normalizations; their ratio is not a field-error estimate. Together they show why checking only the imposed discrete bulk equation would overstate the verification evidence. See [the replay output](grok46_azure_deepseek_2026-10-22/verification_replay_20261022T205012Z/verification_output.txt), [its isolated command record](grok46_azure_deepseek_2026-10-22/verification_replay_20261022T205012Z/MANIFEST.json), and [the check definitions](../runs/P01_deepseek/workspace/checks.py).

A separate source-based check is inexpensive: a = (2/pi) acos(exp(-epsilon/2)). Azure’s submitted a and epsilon disagree with this necessary load–zone relation by about 2.45 percent in a at V=0.001 and 1.52 percent at V=0.1. Its production code estimates the front singular coefficient over the fixed window 0.02 < X < 0.4. Whether that finite-window fit explains the discrepancy is a hypothesis to test by changing the window and mesh independently. The arithmetic uses main-paper equation 17 on PDF page 7 and is recorded separately from the grades in [the calculations](QUALITATIVE_EVIDENCE_CALCULATIONS.json).

## Reproduction and parameter transfer are different claims

The delivered entry points reproduce their submitted results. Azure’s terminal speed is a literal constant in solve.py. Direct DeepSeek stores its two states and terminal speed as constants and treats every speed outside its two stored states as inadmissible. Such case-specific code can satisfy reproduction for the supplied inputs without establishing a solver for other parameters or speeds. The original fixed task does not justify a retrospective penalty for unrequested cases. See the [Azure entry point](../runs/P01_deepseek/workspace/solve.py) and [direct entry point](../runs/P04_deepseek_direct_1017/workspace/solve.py).

A useful next test would change d while holding lambda fixed, ask for the analytical existence condition first, and then sample speeds around the predicted endpoint. The source endpoint relation changes character across d=lambda; a cutoff depending only on lambda would miss this. The new task should explicitly request those regimes and an analytical explanation, with its reference and tolerances independently checked before launch.

## Execution failures need their own explanation

Kimi’s retry delivered no solution after repeated output-limit responses. It provides evidence about this CLI/API execution path, with no analytical work available to grade. DeepSeek’s Azure and direct runs both used reasoning, but their routes, revisions, timings and provider token accounting differ. These few runs cannot isolate a causal effect of backend or reasoning-token quantity on mathematical performance. See the [Kimi diagnostic](../runs/P05_kimi_recheck/FAILURE_DIAGNOSTIC.json), [Azure thinking verification](../runs/P01_deepseek/THINKING_VERIFICATION.json) and [direct thinking verification](../runs/P04_deepseek_direct_1017/THINKING_VERIFICATION.json).

Direct DeepSeek’s answer reports approximately six hours of work, while its host timestamps record about two hours and 58 minutes. Resource comparisons should use recorded controller times and API accounting, with provider differences kept visible. The discrepancy is a reporting observation, not evidence of an agent’s intent. See its [execution record](../runs/P04_deepseek_direct_1017/RUN.json) and answer section 9.

## The grader also needs delivery and evidence checks

The two tool-using Grok attempts on Azure DeepSeek consumed their one-hour and half-hour budgets without delivering a final grade. They did complete useful source reads and numerical checks. A successful replay therefore needs a separate delivery check for the final assessment. These attempts are reviewer execution failures, separate from the solver’s mathematical outcome; see the [first incomplete review record](grok46_azure_deepseek_2026-10-22/INCOMPLETE_REVIEW_RECORD.json) and [continuation failure record](grok46_azure_deepseek_2026-10-22/CONTINUATION_FAILURE_RECORD.json).

The later fixed-packet assessment also timed out without a final grade. Its method differs from the tool-using reviews, so any subsequent assessment must retain that provenance. The [delivery audit](grok46_azure_deepseek_2026-10-22/REVIEW_DELIVERY_AUDIT.json) records all three outcomes and unchanged input hashes. Grader delivery, arithmetic consistency, source coverage, and actual executed verification are study outcomes alongside solver performance.

The next study should retain the score and add an explanation of each failure’s cause, useful partial derivations, the independence of its checks, and sensitivity to untested parameter changes. That would make the results more useful for designing agents and benchmarks than a ranking alone.
