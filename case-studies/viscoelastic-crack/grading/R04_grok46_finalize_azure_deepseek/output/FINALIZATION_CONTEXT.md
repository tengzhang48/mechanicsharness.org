Continuation context: original evidence below is reproduced verbatim. Prior Grok findings are preliminary and may be corrected. The tally is mechanical, not an analytical assessment.


# Original rubric

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


# Supplementary rubric

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


# Published solution and source locators

# Reference solution and source provenance

This is an organizer reproduction of the published Hui and Long solution. The paper authors supply the scientific source; Codex implemented the numerical reproduction and checks. Numerical agreement between these routines is not an independent derivation or review. Source hashes and environment versions are recorded in `SOURCE_MANIFEST.json`.

## Source map

The main article is [Hui and Long 2026](<source paper PDF, not included>), DOI [10.1016/j.jmps.2026.106662](https://doi.org/10.1016/j.jmps.2026.106662). The [supplement](<source paper supplement PDF, not included>) supplies derivations. Page numbers below are PDF page numbers, equal to the printed article page numbers.

| Statement or calculation | Main article locator | Supplement locator |
| --- | --- | --- |
| Strip, loading, moving origin and cohesive boundary conditions | pp. 4–5, Figure 1 and equations 1–4 | S2.1–S2.2 |
| Standard-solid creep compliance and modulus ratio | p. 5, equations 5–7 | S3 |
| Critical separation and intrinsic fracture energy | p. 4 and p. 6, equation 10; p. 12, equation 48 | S3–S4 |
| Stress field and cancellation of the front singularity | pp. 6–7, equations 12 and 17 | S2.3–S2.6 |
| Relaxed far-ahead state and energy release rate | pp. 5–6, equation 9 | S2.1 |
| Damage parameter and bound on energy release | p. 8, equations 22–27 | S13 discusses Mode I limits |
| Strain history and relaxed offset | pp. 11–12, equations 38–47 | S3 |
| Critical-separation integral and loading relation | pp. 12–14, equations 49–69 | S3–S7 |
| Finite terminal speed | pp. 14–15, equations 58–59 and 68–70; p. 18, first paragraph | S6 |
| Main input values and rounded scalar checks | p. 19, Figure 9; related definitions pp. 15–16 | S9 |
| Free energy used to define local dissipation | pp. 22–24, equations 77–84 | S11 |
| Spatial energy integral | p. 24, equations 85–86 | S11 |
| Nonuniqueness of local dissipation from creep data alone | p. 25, paragraph before section 5 | S12 |

The primary parameters, speeds and streamline locations are published values. Sample X coordinates, output schema, accuracy thresholds, resource budget and rubric are organizer evaluation choices. The task uses dimensionless ratios because the paper does not prescribe one necessary set of dimensional material constants for this example.

## Scalar solution

Write \(\lambda=G_\infty/G_0\), \(d=\bar\delta_c\), \(\epsilon=\pi l_c/h\), and \(V=\pi vt_c/h\). The published stress field is

\[
T(Z)=\frac{\tau_2+i\tau_1}{\tau_c}
=1-\frac{2}{\pi}\tan^{-1}\sqrt{\frac{1-e^{-Z}}{e^\epsilon-1}},
\qquad Z=X+iY.
\]

Use the branch analytic in the open upper half strip. The front stress condition selects

\[
a=\frac{2}{\pi}\cos^{-1}(e^{-\epsilon/2}),\qquad R=a^2/d.
\]

Equation 51 determines epsilon:

\[
d=\frac{2}{\pi^2}\int_0^\epsilon
\left[1-(1-\lambda)e^{-\xi/V}\right]
\log\frac{1+\alpha}{1-\alpha}\,d\xi,
\quad \alpha^2=\frac{e^{-\xi}-e^{-\epsilon}}{1-e^{-\epsilon}}.
\]

The integrable logarithm is evaluated with `expm1` and `log1p`, and the substitution \(\xi=\epsilon t^2\) regularizes numerical quadrature. Bisection/Brent bracketing solves the equation. As V tends to zero, R tends to 1 and \(\epsilon\to-2\log\cos(\pi\sqrt d/2)\).

For d greater than lambda, the finite upper speed is selected by equation 59:

\[
1-d=(1-\lambda)k(V_*),\qquad
k(V)=\frac{8}{\pi^2}\sum_{n=1}^\infty
\frac{1}{(2n-1)(2n-1+2/V)}.
\]

The implementation sums this series through an organizer-derived digamma identity and checks it against numerical integration of the paper's equation 58. For the main case, admissible positive speeds are \(0<V<V_*\); \(V_*=13.4235715733\), \(R\to5\), \(a\to1\), and \(\epsilon\to\infty\) at the endpoint. The endpoint is a limiting state, not a finite-zone steady solution. This is a statement of nonexistence within the stated quasistatic steady model; it does not prove a dynamic instability or a specific post-terminal motion.

| V | epsilon | a | R |
| --- | --- | --- | --- |
| 0.001 | 0.5447055551 | 0.4488424892 | 1.0072979004 |
| 0.1 | 0.7407809317 | 0.5148121867 | 1.3251579381 |
| 10 | 6.3333056844 | 0.9731620763 | 4.7352221340 |

These values reproduce the independently published Figure 9 column headings within their three-decimal rounding. `published_fig9.json` preserves the transcribed printed values. Their precision supports a rounded source check, not all the digits above.

## Field and energy calculations

Integration by parts of equation 38 gives the numerical convolution

\[
E(X,Y)=\lambda T(X,Y)+(1-\lambda)\int_0^\infty e^{-u}T(X+Vu,Y)\,du,
\quad E=G_\infty(\gamma_2+i\gamma_1)/\tau_c.
\]

An alternative evaluation solves \(A_{,X}=(A-T)/V\) from the relaxed far-ahead state and sets \(E=\lambda T+(1-\lambda)A\). These use different numerical history evaluations but share the same published stress field. This comparison can expose integration errors, not a shared wrong stress formula. Finite differences additionally check equilibrium, compatibility and the local SLS equation.

With the free energy specified in the task, equations 83–86 give

\[
\widehat\phi_D=\frac{|T-E|^2}{1-\lambda},\qquad
\frac{\Gamma_D}{\Gamma_0}=\frac{1}{\pi dV}
\int_{-\pi}^{\pi}\int_{-\infty}^\infty\widehat\phi_D\,dX\,dY.
\]

The code integrates the upper half and doubles it. Gaussian ordinates use \(Y=\pi t^2\) to resolve the crack-line layer. A BDF history solve supplies strain on each ordinate; adaptive quadrature integrates power over X. The integration routine does not take R or R minus 1 as input. Its result is then compared with the separate scalar solution. Convergence is checked at 16, 32 and 64 ordinates. The V=10 wake is also checked with a longer integration domain and tighter history tolerance.

## Reproduction

From the packet root:

```bash
python private/reference_solver.py agent/case.json private/reference_results.json --energy-ny 64
python private/verify_reference.py
python -m unittest discover -s private -p test_scoring.py -v
```

`verification_results.json` records errors, tolerances, source checks and refinements. `reference_results.json` is the computed numerical key. The scripts require NumPy and SciPy; exact versions are in the source manifest. No finite element or commercial solver is required.

## Limits requiring review

The reference has not received a separate scientific review. It reproduces a published analytical solution and checks its numerical consistency. Stress formulas, root selection and field branches should receive domain review before graded runs. The current evaluator checks a finite set of field samples; the derivation and boundary checks remain essential. No uniform error theorem for the reproduction code is claimed.

This reference concerns small-strain Mode III fracture in an infinite strip. It does not supply a finite-strain Mode I solution, nor a transient solution after the terminal speed. Its difficulty for current agents remains unmeasured.


# Candidate analytical reduction

# Derivation notes (working)

## Constitutive form
Standard linear solid: equilibrium spring $G_\infty$ parallel to a Maxwell branch
(spring $G_1=G_0-G_\infty$, dashpot $\eta$, $t_M=\eta/G_1$).
$$(1+t_M\partial_t)\tau_i=(G_\infty+t_MG_0\partial_t)\gamma_i .$$
Creep compliance $C(t)=\frac1{G_\infty}[1-(1-\lambda)e^{-t/t_c}]$, $\lambda=G_\infty/G_0$,
$t_c=t_M/\lambda$; free energy as given.

## Steady moving frame
$x=x_{lab}-vt$, stationary fields, material derivative $\dot{}=-v\partial_x$.
Eliminating the overstress gives
$$(G_\infty - v t_M G_0 \partial_x)\nabla^2 w = 0 .$$
With decay ahead ($\nabla^2 w\to0$ as $x\to+\infty$) this forces
$$\boxed{\nabla^2 w=0}$$ everywhere. So $w$ is harmonic; viscoelasticity enters only
through the traction/opening relation on $y=0$.

## Crack-plane operator
Half strip $0<y<h$, $\hat w(q,y)=\hat u(q)\sinh(q(h-y))/\sinh(qh)$, so
$\hat w_{,y}(0)=-q\coth(qh)\hat u$.
With $U=G_\infty w(x,0^+)/(\tau_c h)$, $X=\pi x/h$, $K_0=1/(\lambda V)$, $V=\pi vt_c/h$:
$$\hat T_2(Q)=a\,2\pi\delta(Q)-\ell(Q)\hat U(Q),\qquad
\ell(Q)=\frac{\pi}{\lambda}\,Q\coth(\pi Q)\,\frac{iQ-K_0\lambda}{iQ-K_0}.$$
Limits: $K_0\to\infty\Rightarrow \ell=\pi Q\coth(\pi Q)$ (relaxed, $G_\infty$);
$K_0\to0\Rightarrow \ell=\frac\pi\lambda Q\coth(\pi Q)$ (glassy, $G_0$).

## Local (bulk) form used numerically
$$\nabla^2 W=0,\qquad W(X,\pi)=a,\quad X=\pi x/h,\ Y=\pi y/h,$$
$$T_2=\pi W_{,Y}+S_2,\quad T_1=\pi W_{,X}+S_1,$$
$$(\partial_X-K_0)S_1=\tfrac{\pi(1-\lambda)}{\lambda}W_{,XX},\qquad
(\partial_X-K_0)S_2=\tfrac{\pi(1-\lambda)}{\lambda}W_{,XY},$$
with $S_i\to0$ as $X\to+\infty$, i.e.
$S_i(X,Y)=-\tfrac{\pi(1-\lambda)}{\lambda}\int_0^\infty W_{,Xi}(X+u,Y)e^{-K_0u}\,du$.
Boundary data on $Y=0$: $W=0$ ($X\ge0$), $T_2=1$ ($-\epsilon<X<0$), $T_2=0$ ($X<-\epsilon$);
far field $W\to(a/\pi)Y$ (ahead), $W\to a$ (behind).
Unknowns $a,\epsilon=\pi l_c/h$ from: bounded $T_2$ at $X=0$ and $2U(-\epsilon)=\bar\delta_c$.

# Previous Grok preliminary findings

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


# Executed commands

=== 2026-10-22T19:00:32Z score_submission.py submitted results ===
CMD: python /workspace/reference/score_submission.py /workspace/reference/reference_results.json /workspace/candidate/results.json /workspace/audit/submitted_score_report.json
{
  "all_required_numeric_gates_pass": false,
  "scientific_review": "pending",
  "study_pass": null,
  "scope": "JSON accuracy only; derivation, integration provenance, reproducibility and access audit require review"
}
EXIT:1
=== 2026-10-22T19:08:34Z candidate solve.py reproduce ===
CMD: python solve.py case.json /workspace/audit/reproduced_results.json
CWD: /workspace/candidate
wrote /workspace/audit/reproduced_results.json
EXIT:0
=== 2026-10-22T19:10:11Z solve.py finished ===
=== 2026-10-22T19:16:30Z score reproduced_results.json ===
CMD: python /workspace/reference/score_submission.py /workspace/reference/reference_results.json /workspace/audit/reproduced_results.json /workspace/audit/reproduced_score_report.json
{
  "all_required_numeric_gates_pass": false,
  "scientific_review": "pending",
  "study_pass": null,
  "scope": "JSON accuracy only; derivation, integration provenance, reproducibility and access audit require review"
}
EXIT:1


# Mechanical numeric tally

{
  "basis": "Mechanical count from Grok executed frozen scorer; no analytical judgment",
  "all_checks": 251,
  "failed_checks": [
    {
      "name": "V=0.001 epsilon",
      "passed": false,
      "actual": 0.577,
      "expected": 0.5447055550934302,
      "absolute_error": 0.03229444490656974,
      "tolerance": 0.0027235277754671513
    },
    {
      "name": "V=0.001 X=0.0 y/h=0.1 T1",
      "passed": false,
      "actual": -0.2010047360823951,
      "expected": -0.21875036317973948,
      "absolute_error": 0.01774562709734437,
      "tolerance": 0.005
    },
    {
      "name": "V=0.001 X=0.0 y/h=0.1 E1",
      "passed": false,
      "actual": -0.20052351775387467,
      "expected": -0.21823543828839817,
      "absolute_error": 0.0177119205345235,
      "tolerance": 0.005
    },
    {
      "name": "V=0.1 epsilon",
      "passed": false,
      "actual": 0.7686,
      "expected": 0.7407809316699946,
      "absolute_error": 0.02781906833000536,
      "tolerance": 0.003703904658349973
    },
    {
      "name": "V=0.1 X=-2.0 y/h=0.1 T1",
      "passed": false,
      "actual": -0.27722808531132476,
      "expected": -0.27217093759850347,
      "absolute_error": 0.005057147712821286,
      "tolerance": 0.005
    },
    {
      "name": "V=0.1 X=-2.0 y/h=0.1 E1",
      "passed": false,
      "actual": -0.29468711699114897,
      "expected": -0.2890637938540239,
      "absolute_error": 0.005623323137125047,
      "tolerance": 0.005
    },
    {
      "name": "V=0.1 X=0.0 y/h=0.1 T1",
      "passed": false,
      "actual": -0.1791389180331695,
      "expected": -0.19215423942861287,
      "absolute_error": 0.013015321395443369,
      "tolerance": 0.005
    },
    {
      "name": "V=0.1 X=0.0 y/h=0.1 E1",
      "passed": false,
      "actual": -0.1469096223975444,
      "expected": -0.15755108930188164,
      "absolute_error": 0.010641466904337232,
      "tolerance": 0.005
    }
  ],
  "field_component_count": 180,
  "field_component_pass_count": 174,
  "scalar_checks": [
    {
      "name": "critical speed",
      "passed": true,
      "actual": 13.44,
      "expected": 13.42357157330645,
      "absolute_error": 0.01642842669354927,
      "tolerance": 0.06711785786653225
    },
    {
      "name": "V=0.001 epsilon",
      "passed": false,
      "actual": 0.577,
      "expected": 0.5447055550934302,
      "absolute_error": 0.03229444490656974,
      "tolerance": 0.0027235277754671513
    },
    {
      "name": "V=0.001 a",
      "passed": true,
      "actual": 0.44940200380307205,
      "expected": 0.4488424891698805,
      "absolute_error": 0.0005595146331915379,
      "tolerance": 0.0022442124458494027
    },
    {
      "name": "V=0.001 R",
      "passed": true,
      "actual": 1.009810805111082,
      "expected": 1.0072979004210714,
      "absolute_error": 0.002512904690010531,
      "tolerance": 0.005036489502105357
    },
    {
      "name": "V=0.1 epsilon",
      "passed": false,
      "actual": 0.7686,
      "expected": 0.7407809316699946,
      "absolute_error": 0.02781906833000536,
      "tolerance": 0.003703904658349973
    },
    {
      "name": "V=0.1 a",
      "passed": true,
      "actual": 0.5152246957660931,
      "expected": 0.5148121867403258,
      "absolute_error": 0.00041250902576728254,
      "tolerance": 0.002574060933701629
    },
    {
      "name": "V=0.1 R",
      "passed": true,
      "actual": 1.3272824356363158,
      "expected": 1.3251579380817804,
      "absolute_error": 0.002124497554535365,
      "tolerance": 0.006625789690408902
    },
    {
      "name": "V=10.0 epsilon",
      "passed": true,
      "actual": 6.34,
      "expected": 6.333305684439911,
      "absolute_error": 0.006694315560088704,
      "tolerance": 0.03166652842219956
    },
    {
      "name": "V=10.0 a",
      "passed": true,
      "actual": 0.9729338032023915,
      "expected": 0.9731620763214321,
      "absolute_error": 0.00022827311904061176,
      "tolerance": 0.004865810381607161
    },
    {
      "name": "V=10.0 R",
      "passed": true,
      "actual": 4.733000927069349,
      "expected": 4.735222133951204,
      "absolute_error": 0.002221206881855231,
      "tolerance": 0.02367611066975602
    }
  ],
  "bulk_checks": [
    {
      "name": "V=0.001 bulk dissipation",
      "passed": true,
      "actual": 0.008924659279008866,
      "expected": 0.007297898854430374,
      "absolute_error": 0.0016267604245784918,
      "tolerance": 0.01
    },
    {
      "name": "V=0.1 bulk dissipation",
      "passed": true,
      "actual": 0.3262022658682324,
      "expected": 0.32515792310981306,
      "absolute_error": 0.0010443427584193588,
      "tolerance": 0.01
    },
    {
      "name": "V=10.0 bulk dissipation",
      "passed": true,
      "actual": 3.738413275977389,
      "expected": 3.7352219340297856,
      "absolute_error": 0.003191341947603288,
      "tolerance": 0.03735221934029786
    }
  ],
  "energy_checks": [
    {
      "name": "V=0.001 energy balance",
      "passed": true,
      "error": 0.0008775364925665932,
      "tolerance": 0.01
    },
    {
      "name": "V=0.1 energy balance",
      "passed": true,
      "error": 0.0008138205848896749,
      "tolerance": 0.01
    },
    {
      "name": "V=10.0 energy balance",
      "passed": true,
      "error": 0.001143534301268622,
      "tolerance": 0.01
    }
  ],
  "reproduction_match": "Grok retained shell output reports submitted==reproduced True; prior independent organizer replay also matches. See actual tool evidence and commands.log."
}
