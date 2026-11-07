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
