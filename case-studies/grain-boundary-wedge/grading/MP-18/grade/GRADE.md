# Grading profile

## 1. Answers

| Part | Judgment | Reason |
|---|---|---|
| (a) | **Partly right** | The eigenmode expansion, exponential relaxation, and reduced eigenproblem are correct, but both requested rates are wrong: 188.1648 instead of 11.74265 (16.02 times the reference), and 1089.2984 instead of 264.6348 (4.12 times the reference). Neither the exact nor approximate grading band admits these values. |
| (b) | **Right** | A traction-free, crack-like grain-boundary wedge with its tip at the interface is correct, and 1.1215 is in the 1.09–1.15 intensity band. The additional numerical arguments for this number are defective; see the side errors below. |
| (c) | **Wrong** | The formula t_relax = t0/lambda1 is correct, but 0.00531 t0 is outside 0.05–0.15 t0; the reference reciprocal is 0.08515965 t0. |

No answer is contested. The objections to the stated model do not withstand the checks in section 9.

## 2. Supporting state

**Wrong.** The requested intermediate quantity is the slowest traction mode at half depth, normalized by f(1). The submission gives −0.74478; the reference is approximately +0.8020, within a 2% band. This is a different, sign-changing discrete mode, not the physical first mode. The nondimensional time and the equation s_tau = integral S s''' are derived correctly overall, despite the dimensional slip identified below.

## 3. Method

**PV hat-function collocation with finite-difference third derivatives, eigenvalue filtering, and Richardson extrapolation; `general_or_wrong` as implemented.** Rows are pointwise collocation equations, rather than Galerkin test-function integrals. `solve_final.py:39` forms pi times `M_stmt` times `D3`; the pi restores the correct unit coefficient. `solve_final.py:58` drops the near-zero eigenvalue, then extrapolates the remaining spectrum.

Numerical singular-integral solution is an appropriate route, but this implementation does not impose both interface conditions correctly. Refining an incorrectly constrained matrix and obtaining stable numbers does not validate the continuum problem. The first mode has a singular third derivative at the tip, making correct endpoint treatment especially important.

For an independent audit, `scratch/grade_audit.py` expands f''' in odd Chebyshev polynomials divided by sqrt(1−x²), integrates three times to impose f(0)=f'(1)=f''(1)=0, handles the Cauchy principal-value part analytically, and integrates the regular image term by Gauss quadrature. It uses the stated kernel and stated time unit:

| Number of basis functions | lambda1 | lambda2 | f1(0.5)/f1(1) |
|---|---:|---:|---:|
| 16 | 11.74355207 | 264.70965044 | 0.80195234 |
| 32 | 11.74273323 | 264.65473674 | 0.80195620 |
| 64 | 11.74266127 | 264.63697838 | 0.80195654 |
| 96 | 11.74266295 | 264.63480392 | 0.80195653 |

These results support the key. No error in the key was found.

## 4. Errors

The list separates failed requested quantities from implementation and reasoning errors, and groups repeated occurrences of each error. Superseded exploratory answers are not additional final answers.

### Decisive: 6

| ID | Evidence | Error and effect |
|---|---|---|
| D1 | `results.json:2`; `ANSWER.md:9`, “188.2” | Wrong leading decay rate; it fails the key by a factor of 16.02. |
| D2 | `results.json:3`; `ANSWER.md:10`, “1089” | Wrong next decay rate; it fails the key by a factor of 4.12. |
| D3 | `results.json:4`; `ANSWER.md:11`, “−0.7448” | Wrong first-mode shape, including its sign. |
| D4 | `results.json:6`; `ANSWER.md:136`, “5.3×10⁻³ t0” | Wrong relaxation time, inherited from D1; neither proposed normalization puts it in the reference band. |
| D5 | `d3_uniform.py:51`, `:57`, `:62–67`; `ANSWER.md:49` | The last two rows use node indices shifted from the interior row coordinates and do not impose f''(1)=0. The claimed ghost elimination is also incorrect. With h=1/(N+1), q=x(2−h−x) is a discrete null vector although q''(1)=−2. Even the cubic f=1−(1−x)³, which satisfies all three required BCs, gives D3 f=5.5 rather than 6 in the last row for every tested N. This changes the eigenproblem. |
| D6 | `ANSWER.md:52–55`, “a discretisation artefact”; `solve_final.py:58` | Recognizing a defective null mode does not justify deleting it and declaring the rest of that matrix the physical spectrum. The continuum third derivative with the three BCs has no quadratic null space. The BC defect must be repaired. The claim of no physical mode below 188 is refuted by the independently constrained solution above. The code's “x-independent” null-mode description (`solve_final.py:26`) is also false: the null vector tends to 2x−x². |

### Side: 9

| ID | Evidence | Error and effect |
|---|---|---|
| S1 | `ANSWER.md:143–146`, “a positive wedge increase[s] the traction” | The sign test uses positive constant u', not a positive wedge opening. With the tip anchored, positive u=1−x has negative u'. The negative integral of S against 1 therefore produces stress relief for a positive wedge, not the claimed sign defect. The final solver nevertheless uses the stated sign. |
| S2 | `ANSWER.md:102–114`, “one real ambiguity”; “Both are internally consistent” | There is no pi ambiguity with the specified t0. The two dimensional coefficients multiply to one, as both the key and the solver's final derivation state. Dividing the operator by pi requires changing the time unit. The headline code uses the correct coefficient, so this erroneous alternative is a side error. |
| S3 | `ANSWER.md:46`, `:85`, “< 1×10⁻¹⁰” | The named 800-node validation is overstated. Rerunning `phys_oper.py` on the copy gives maximum absolute difference 5.7263e−10 and relative difference 1.1339e−9. The kernel agreement is good but does not support the reported tolerance. |
| S4 | `ANSWER.md:53`, `:88`, “sigma_min/sigma_max ≈ 10⁻¹⁰ at N=160” | Direct SVD of the actual final D3 gives approximately 8.7e−18, with rank 159. Rank deficiency is real, but the stated singular-value ratio is inaccurate. |
| S5 | `ANSWER.md:90`, “same leading nonzero” | The full chain spectrum retains a slower nonzero rate: at N=160 its first two rates are about 0.22494 and 60.65628 in the old /pi convention. The scripts track the mode near 60 and disregard the slower mode, rather than demonstrate agreement of the leading nonzero eigenvalue. Both routes also share incorrect endpoint treatment. |
| S6 | `ANSWER.md:91`, “a real mode would survive a grid transform”; `eta_transform.py:31–35`, `:70–71` | The alleged coordinate-transform confirmation changes the operator even at p=1: the first derivative uses the first unknown instead of its right neighbor. At N=160, p=1 introduces a growing eigenvalue of approximately +628701.5 in the old convention. Filtering to positive decay rates hides it. Disappearance of a slow mode under this altered, unstable discretization is not evidence against a real mode. |
| S7 | `kinf.py:15–19`, `:53–56`; `ANSWER.md:128–130` | K is extracted from stress ahead of the tip (z>h), not residual normal traction on the traction-free faces (z<h). The code samples z=1−r, where the continuum face traction is zero. Its saved tails approach small values, not 1.1215. This same face/ahead-of-tip confusion does not justify the proposed zero-intensity reading. The final handbook factor remains correct. |
| S8 | `kinf.py:79–88`; `ANSWER.md:126–128` | The “singular column” is exactly M times a vector of nodal samples; it is already in the span of the hat columns. The N-by-(N+1) augmented system is underdetermined, and [−col,1] is an exact null vector. A nonzero minimum-norm coefficient C and a tiny residual do not demonstrate a physical singular density or verify K. |
| S9 | `submission/PROCESS.md:19`, “fits stable to 4 dp”; `ANSWER.md:75` | The pasted lambda1 intercepts 188.1675, 188.1648, and 188.1634 are not stable to four decimal places. The final answer's wider ±0.5 uncertainty is retained; this unsupported precision claim is not another numerical answer. |

### Presentation: 5

| ID | Evidence | Error and effect |
|---|---|---|
| P1 | `ANSWER.md:32`, “u-hat_t = −(Omega/2) j_z” | For u-hat=u/h, this intermediate equality needs an additional 1/h. The subsequent h³ coefficient and reduced eigenproblem are correct. |
| P2 | `ANSWER.md:39`, “a smooth PV kernel” | S has a Cauchy pole, not a smooth kernel. Principal-value handling is correctly recognized elsewhere. |
| P3 | `results.json:7`; `solve_final.py:91`, “Galerkin” | The delivered hat matrix is used at pointwise collocation nodes; no Galerkin test-function projection is performed. |
| P4 | `ANSWER.md:15`, `:89–100`; `submission/PROCESS.md:18` | The spectrum's cutoff and pi convention are mixed. Multiplying the pasted /pi rates by pi gives approximately −3312.8, −1089.9, and −188.3; only two are in (−1200,0), and only two are in (−1100,0). The third mode is not “below ~1200.” The reported three-mode count applies to the old convention's wider physical interval. |
| P5 | `kinf.py:45–46`, “integral = pi” | M already contains 1/pi, so solving M a=pi gives the raw integral pi². The later division by pi in the residual makes the collocation residual internally consistent, but the stated raw-integral normalization is wrong. |

## 5. Checks

**11 reported checks; 7 could fail in a relevant verification role.** This counts C1–C8, the separate rate-convergence/Richardson comparison, and the two equilibrium/intensity checks claimed through `kinf.py`. External automatic checks and superseded exploratory runs are not added to the count. “Could fail” credits a falsifiable comparison; it does not mean it passed or excludes every systematic error.

| Check | Tests the deciding quantity? | Could fail? | Assessment |
|---|---|---|---|
| C1: analytic kernel versus direct PV quadrature | Tests the elasticity matrix, not the full eigenproblem | **Yes** | A different kernel implementation can disagree. Good agreement, but the claimed tolerance fails on rerun. |
| C2: D3 on cubic/quadratic polynomials | Tests bulk derivative consistency; misses the deciding BC rows | **Yes** | Wrong bulk stencils can fail. Interior agreement is real; the last row gives −159.5 on x³ and −115.5 on x² at N=10, so this does not verify the claimed boundary elimination. |
| C3: S(0,y)=0 | Tests the pointwise surface-kernel identity | **Yes** | A wrong image term can fail this known identity. It is true for the stated kernel, but does not impose interface conditions or verify eigenvalues. |
| C4: D3 rank/null space | Diagnoses an algebraic defect, not continuum accuracy | **No** | The expected rank-deficient result accepts a matrix that has lost a BC; interpreting it as a passed solution check is backwards. |
| C5: full discrete-spectrum enumeration | Tests whether a smaller nonzero mode in that matrix was missed | **Yes** | It could expose wrong mode selection within that matrix. It cannot exclude modes missing from the continuum discretization, and the stated interval/count is inconsistent. |
| C6: chain versus D3 rate | Tests a selected rate, not the full leading spectrum | **Yes** | These different derivative approximations can disagree. This is a limited consistency comparison, not independent validation; the claimed leading-nonzero agreement is false. |
| C7: coordinate transform | Does not preserve the stated operator/BCs | **No** | Even p=1 changes a row and introduces instability; discarded growing and marginal modes prevent the claimed diagnostic from validating the solution. |
| C8: first-mode shape versus N | Tests convergence of the selected discrete shape | **Yes** | Nonconvergence could fail it. The supplied numbers reproduce, but convergence is to the wrong mode/problem. |
| Rate convergence and Richardson fits, `ANSWER.md:61–76` | Tests numerical stability of the selected rates | **Yes** | Resolution or fit instability could fail the comparison. Stable intercepts cannot cure the BC defect, and four-decimal stability is not shown. |
| Equilibrium K-limit/literature comparison, `kinf.py:50–69` | Samples face traction instead of the field defining K | **No** | It prints small tails and then prints the reference constant without an agreement test. It cannot verify the stated intensity. |
| Singular-density amplitude and residual, `kinf.py:71–88` | Does not identify a unique singular amplitude | **No** | The extra column is dependent, and the existing square matrix already fits the equilibrium RHS. A tiny residual is an algebraic consequence, not a singularity check. |

## 6. Claims and reproducibility

**Partly supported.** `AUTO_CHECKS.json` records a successful reproduction and 5/5 printed numerical claims. A fresh `python3 solve_final.py` run in `scratch/reproduction/` exits successfully and produces a `results.json` byte-for-byte identical to the submitted file. The five reported values and the main convergence table therefore reproduce. Reproducibility does not establish correctness: the automatic rate, shape, and time checks fail, while the K band and the two internal ordering/reciprocal checks pass.

The reproduction assigns `K_inf_norm = 1.1215` at `solve_final.py:79`; it does not calculate K or execute C1–C7. The separately delivered files support some limited comparisons, but not the stronger validation claims: the C1 tolerance, C4 singular-value ratio, C5 interval/count, C6 leading-mode claim, C7 preservation of the operator, and purported K confirmations have the mismatches explained above. The correct handbook value for K is credited without granting a numerical verification that the code does not provide.

## 7. Questions and missing information

**Four entries** appear in `ANSWER.md:143–150`:

1. Kernel sign: worth checking, but its test confuses opening with opening slope and does not identify missing information.
2. Initial condition versus surface BC at t=0: sensible. Read the initially uniform interior stress with the surface condition imposed for t>0; an endpoint incompatibility/transient does not alter the eigenvalues.
3. The 1/pi convention: the stated t0 already resolves it; no clarification is needed.
4. Finite depth/half-space and cohesive versus crack interpretation: the geometry and crack-like limit are specified. No cohesive law is supplied, so this does not establish an alternative model.

## 8. Source

**None named.** Neither the final answer nor its supporting files name Gao, Zhang, Nix, Thompson, Arzt, their paper, or another source author/title for this problem. `kinf.py:66` says only “Literature check -- single edge crack of depth h in a half plane”; that is not a named source. The author list in `KEY.md` is not credited to the solver.

## 9. Challenges

| Claim | Argument and checkability | Assessment |
|---|---|---|
| The stated kernel has an unphysical sign (`ANSWER.md:143–146`). | It computes integral S against a positive constant slope and observes that it is negative. **Checkable:** evaluate that PV integral and reconstruct an opening satisfying u(1)=0. | The computed integral is negative (−2.36917 at x=0.25 and −1.12361 at x=0.5), but the positive-slope test corresponds to negative opening. For positive opening u=1−x, the given traction law relaxes stress. The objection is based on the wrong physical test and does not warrant a contested grade. |
| A pi factor is ambiguous (`ANSWER.md:102–114`). | It correctly derives coefficient one, then offers its other scripts' /pi operator as an equally consistent convention. **Checkable:** nondimensionalize the supplied equations and time scale. | The fixed time unit yields coefficient one. Even dividing both the reference and submitted rates by pi leaves their factor-of-16 and factor-of-4 discrepancies and the wrong shape. The alternative does not rescue (a) or (c). |
| A cohesive/zero-intensity reading could instead give K_inf=0 (`ANSWER.md:128–130`). | It associates traction-free conditions including the tip with zero intensity. **Not checkable as a distinct alternative model as stated:** no cohesive traction law or additional tip condition is specified. | Traction-free crack faces are compatible with nonzero intensity ahead of the tip. The key's 1.1215 limit is appropriate to the stated geometry. The speculative alternative does not warrant a contested grade; see also S7. |

There is also a **valid, checkable criticism of the delivered draft `march.py`**, at `ANSWER.md:159–160`: differentiating its primitives gives a kernel different from the stated S. This diagnosis is supported by `diag_kernel.py` and `diag_kernel2.py` and is not counted as an error. It concerns an internal exploratory implementation, not a defect established in the problem, key, or published paper. The final code restores the stated kernel but retains the boundary-discretization defect. There is no verified challenge requiring an alternative answer grade.

```json
{"answers": {"a": "partly", "b": "right", "c": "wrong"},
 "supporting_state": "wrong",
 "method": "PV hat collocation with finite-difference D3, null-mode filtering and Richardson extrapolation",
 "method_class": "general_or_wrong",
 "errors": {"decisive": 6, "side": 9, "presentation": 5},
 "checks_reported": 11,
 "checks_could_fail": 7,
 "claims_supported": "partly",
 "questions_listed": 4,
 "source_named": "none",
 "challenges": [{"claim": "Stated kernel has an unphysical sign", "checkable": true}, {"claim": "Fixed t0 leaves a pi normalization ambiguity", "checkable": true}, {"claim": "Cohesive interpretation could require K_inf=0", "checkable": false}]}
```
