# Grading profile

## 1. Answers

- **(a) Right.** The surface first relaxes, followed by exponential modal decay dominated by the first mode; `lambda1 = 11.7427` and `lambda2 = 264.6342` satisfy the key's exact-solution bands. The extra spectral-asymptotics claim is wrong but does not enter these computed rates.
- **(b) Right.** The limiting boundary is a traction-free, crack-like diffusion wedge, with `K_inf_norm = 1.1215`, inside the key's 1.09–1.15 band. The two final-reproducer estimates are 1.121522 and 1.121516.
- **(c) Contested only as to the definition of relaxation; right against the key.** The delivered `relaxation_time_t0 = 0.0852` agrees with `1/lambda1 = 0.0851595` using the unrounded eigenvalue and the 0.05–0.15 band. Against the solver's alternative 95% criterion, approximately 0.27 t0 is the appropriate mid-boundary threshold time. Both conventions are identified, so this is not a numerical error; see item 9.

## 2. Supporting state

**Right for the quantities requested.** The reproduced first mode has `f(0.5)/f(1) = 0.80196`, reported as 0.8020 and matching the key within 2%. Its other displayed values, 0.21316 at x = 0.1 and 0.99619 at x = 0.9, are consistent with the stated profile. The nondimensional equations `s = 1 - A[w]`, `w_tau = -s_xxx`, and the eigenproblem `A[f'''] = -lambda f` have the correct signs and time normalization. The three traction boundary conditions are built into the reconstruction of f. The singular final-state density and its K normalization are implemented correctly in the final reproducer.

## 3. Method

**Chebyshev singular-integral collocation; method class: right.** The solver subtracts the principal-value singularity, integrates in the variable `xi = 1 - t^2`, and represents `g = f'''` with the integrable tip factor `x/sqrt(1-x)`. Integrating g with `f(0) = f'(1) = f''(1) = 0` gives a generalized eigenvalue problem. This is the careful singular-integral numerical route called for by the key. The final crack calculation and a separate differentiation-based traction pencil provide useful verification. This is not an unnecessarily general method. The erroneous explanation of spectral growth is not used to obtain the answers.

## 4. Errors

**Decisive: 0.** No identified defect changes a requested numerical answer. The following side errors include defects in retained development files; those files are distinguished from the working final reproducer.

1. **Side — false spectral law and false comparison.** `ANSWER.md:29–30,125–128` says the rates follow `n^4.5` and quotes supposedly theoretical ratios 22.6, 5.14, 2.83, 2.17, 1.85. The actual ratios are 22.5361, 5.1295, 2.8327, 2.1657, 1.8472, but `(n+1)^4.5/n^4.5` gives **22.6274, 6.2003, 3.6494, 2.7296, 2.2715**. Only the first is close. Moreover, the interior principal symbol of `-A d^3/dx^3` is `pi |k|^3`: the Hilbert-transform part has order zero, and endpoint square roots do not turn it into an operator of order 4.5. `solve.py:159` prints the asserted agreement without computing a comparison. This is wrong reasoning accompanying right eigenvalues.
2. **Side — overstated residual bound and grid coverage.** `ANSWER.md:115–116` claims an independent 2000-node grid and a bound of about `7e-7` everywhere. `solve.py:167–177` evaluates just seven spatial points; 2000 is the quadrature-node count. An actual independent 2000-point spatial grid on `[0.0005,0.9995]`, with the same 2000-node quadrature, gives a residual **5.1704e-6 at x = 0.0005**. The residual is still small, but the claimed global bound is false.
3. **Side — imposed surface value presented as emergent verification.** `ANSWER.md:92–95,133–134` calls `w(0) = -1/4` emergent from the crack solve. `solve.py:68–70` and `crack.py:39–42` explicitly add that value as a least-squares boundary row. Its small row residual is a consistency check, not an independent discovery or independent confirmation of the kernel. The analytic surface-limit argument itself is correct.
4. **Side — circular decay-rate verification.** `ANSWER.md:117–118` says the fitted slope could expose a spurious lambda1. `solve.py:258–261` constructs the time history with the very same eigenvalues it compares against. Any inserted positive slow rate with nonzero fitted amplitude will reproduce itself. This can inspect whether a mode participates in the reconstruction, but cannot establish that its decay rate solves the physical problem.
5. **Side — overstated uniform-grid agreement.** `ANSWER.md:104–106` says the rates are unchanged for uniform collocation; `submission/PROCESS.md:5` claims five-digit agreement. The delivered `dbg13.py` study gives, at uniform `Nv = 48`, **11.7427495, 264.72523, 516.3776**, rather than the corresponding Chebyshev values 11.7426681, 264.63422, 1357.4493. Uniform `Nv = 32` does give a useful close comparison, 11.7426229 and 264.63578, but does not justify blanket agreement across the supplied study or absence of spurious modes.
6. **Side — excluded endpoint included in development checks.** `check1.py:7–13` sends x = 1 into formulas defined for interior points, and reports `max rel err ... = nan` at both N = 40 and 80. `dbg10.py:11–25` similarly includes x = 1 in its supposedly interior collocation set and produces `A finite: False`. The final reproducer excludes the endpoints correctly. These retained checks cannot support successful global validation as written.
7. **Side — development K estimator does not take the tip limit.** `crack.py:61–70` averages `s(1+d)*sqrt(2d)` at d = 0.02, 0.05, 0.1, 0.2 instead of extrapolating to d = 0. Its printed “kappa(fit ahead)” is approximately **1.167424**, while its density estimate is 1.121522. This estimator is inaccurate for K; `solve.py:73–80` uses smaller distances and extrapolation and fixes the problem for the delivered answer.
8. **Side — unfinished mode-shape helper.** `eig_g.py:76–79` promises a normalized mode evaluator, uses the identically zero expression `lam1 - lam1` to select a column, and returns `None`. The final reproducer computes the shape with its separate working `f_of` function.
9. **Side — eigenvalue/vector correspondence broken in a development interface.** `eig.py:45–51` removes infinite eigenvalues without removing the matching columns of V. At N = 24 it returns 22 eigenvalues alongside 25 eigenvector columns, so indexing the returned arrays together is unsafe. The final reproducer keeps the full eigenvalue array when selecting its matching vectors.

**Presentation: 3.**

1. `ANSWER.md:66` mixes dimensional factors into its written K-coefficient relation. For the declared dimensionless density `w = v/sqrt(1-x)`, the code correctly uses `K/(sigma0 sqrt(pi h)) = pi sqrt(2) |v(1)|`; that is the relation the prose should state.
2. `submission/PROCESS.md:6` calls the fitted logarithm `log10` while quoting slope `-lambda1`. The code uses the natural logarithm, correctly. A base-10 slope would be `-lambda1/ln(10)`, approximately -5.0998.
3. `gbcore.py:13` describes `diff_mat` as a “General-nodes” matrix, but lines 15–17 use Chebyshev–Lobatto weights. It is correct on the grids actually used; on nodes `[0,0.2,0.5,1]`, differentiating `x^2` has a maximum error of 0.4. The documented domain is too broad.

## 5. Checks

Counting the **seven check families explicitly reported in ANSWER.md**, five provide substantive checks that could fail. A printed comparison need not be an automated assertion to receive credit.

| Reported check | What it tests and its limitations | Could fail meaningfully? |
|---|---|---|
| Kernel versus analytic `A[1]` | Reproduced agreement for constant density; checks the constant PV/log term and image-term integral, not the nonconstant subtraction machinery or the eigenvalues directly. | Yes |
| Crack limit | Density and extrapolated ahead-tip K agree with each other and the external edge-crack value. They directly test answer (b). The imposed surface value is not independent evidence. | Yes |
| Eigenvalue convergence | The basis-size table reproduces, and `dbg13.py` supports quadrature stability over 200–1600 nodes. Directly tests resolution of lambda1 and lambda2. Uniform-grid evidence is less strong than claimed. | Yes |
| Independent traction pencil | The N = 40, 60, 96 results and extrapolations reproduce for both rates and the half-mode shape. It uses different unknowns and differentiation, though it shares the kernel. | Yes |
| Off-collocation eigenpair residual | Checks the eigen-equation at seven points with finer quadrature and could reveal collocation artifacts. The claimed global bound fails the broader audit. | Yes |
| Late-time slope from mode expansion | Confirms participation of the chosen slow mode and arithmetic of the expansion; its rate is inserted into that expansion. It cannot reject an incorrect participating rate. Threshold times are approximation outputs, not independent rate checks. | No, for the deciding rate |
| Spectral asymptotics | No comparison is computed in the reproducer, and the claimed theoretical ratios are false. | No, as delivered |

The transient threshold values also depend on truncation and initial-condition fitting: the ten-mode fit has about 0.103 maximum initial-profile error on its 19 fitting points. A 16-mode, 399-point fit moves the two thresholds from about 0.1005/0.2705 to 0.1017/0.2717. This supports the answer's approximate 0.10/0.27 statements, not four-decimal accuracy for the physical transient.

## 6. Claims and reproduction

**Partly supported overall; all five required numbers are supported.** `AUTO_CHECKS.json` reports exit 0, 5/5 numerical claims printed, and all seven automatic value/principle checks passed. My rerun of `python3 solve.py` on `scratch/submission` also exits 0 and produces a `results.json` identical to the submitted file. The pasted kernel, crack, convergence, independent-pencil, and late-time numerical outputs reproduce. Harmless temporary divide-by-zero warnings in the interpolation utilities do not spoil those final values.

The automatic layer does not establish the spectral-asymptotics assertion, the global residual bound, independence of the surface/decay checks, or the uniform-grid generalization. Those mismatches are documented above. The quadrature/uniform-grid study resides in `dbg13.py`, rather than being rerun by the main reproduction command. Evidence is retained in `scratch/solve-output.txt`, `scratch/grid-check-output.txt`, `scratch/audit-output.txt`, `scratch/check1-output.txt`, and `scratch/dbg10-output.txt`.

## 7. Questions and missing information

**Five numbered questions/assumptions are listed** (`ANSWER.md:132–144`). The surface integral's limiting interpretation, location/definition of K, relaxation criterion, incompatible initial/surface corner data, and stated modeling assumptions are sensible. They are resolved or explicitly adopted; none is presented as missing information that blocks a solution. In particular, the endpoint and initial-data notes explain standard limiting interpretations without disputing the computed eigenvalues.

## 8. Source naming

**No source paper or author is named.** The nearest attribution is “handbook edge-crack in half-space, uniform face pressure” (`ANSWER.md:91`), which gives neither a handbook title nor an author. The Gao et al. paper named in the key is not named in the submission. Accordingly, `source_named` is `none`.

## 9. Challenges

There is one modest, checkable **definition challenge**, not a challenge to the governing equations or a published result. The solver writes: **“The problem does not define the measure, so both are given”** (`ANSWER.md:139`). Its argument distinguishes the dominant-mode e-folding time, 0.0852 t0, from the mid-boundary 95% threshold, approximately 0.27 t0. The distinction can be checked by rerunning the eigenvalue calculation and threshold reconstruction; the latter's precision limitations are stated above. This is a defensible ambiguity in an operational relaxation criterion, so it is not counted as an error and (c) is marked contested as to definition: **right against the key's e-folding convention; right, approximately, against the solver's explicitly stated threshold convention**.

The two endpoint concerns are clarifications rather than competing solution claims: `S(0,xi) = 0` pointwise does not commute with the limiting singular integral (`A[1](0+) = -4`), and the uniform initial traction conflicts with the zero surface value at the initial corner, requiring an immediate surface transient. Both are checkable and correctly explained, but neither is argued to change the key's eigenvalues or crack limit. No challenge to a named paper is made, and no error in the key's required results was found.

```json
{"answers": {"a": "right", "b": "right", "c": "contested"},
 "supporting_state": "right",
 "method": "Chebyshev singular-integral collocation",
 "method_class": "right",
 "errors": {"decisive": 0, "side": 9, "presentation": 3},
 "checks_reported": 7, "checks_could_fail": 5,
 "claims_supported": "partly",
 "questions_listed": 5,
 "source_named": "none",
 "challenges": [{"claim": "The relaxation criterion is unspecified; e-folding and 95% threshold times differ", "checkable": true}]}
```
