# Grade

## 1. Answers

| Part | Judgment | Reason |
|---|---|---|
| (a) | **Partly right** | Correctly describes modal relaxation with a dominant late exponential, but λ₁ = 0.59877 and λ₂ = 135.9101 are wrong: the key gives 11.74265 and 264.6348 (errors −94.90% and −48.64%). Neither exact nor approximate grading bands are met, and the discretization does not substantiate the claimed evolution. |
| (b) | **Partly right** | Correctly identifies a traction-free, crack-like equilibrium wedge and a square-root tip singularity, but K∞/(σ₀√(πh)) = 3.11614 is outside 1.09–1.15; the reference is about 1.1215. The numerical wedge is not actually traction-free. |
| (c) | **Wrong** | 1.67008 t₀ is outside 0.05–0.15 t₀ and is about 19.61 times the reference 0.08516 t₀. Taking the reciprocal of the submitted λ₁ is internally consistent but inherits its error. |

There is no substantiated challenge requiring a contested grade.

## 2. Supporting state

**Wrong.** The requested intermediate mode value is 0.74025 rather than 0.8020 (−7.70%, outside the 2% band). Normalizing the numerical eigenvector to f(1) = 1 is done correctly, but its interface curvature is wrong: at n = 26, the delivered second-difference operator gives f″(1) ≈ −2.16185 instead of zero. The equilibrium solve has max |Ew∞ − 1| ≈ 0.54625 at n = 30. Its fitted coefficient and opening therefore do not establish the required equilibrium state. The additional asserted rates 767, 2280, 4942 also disagree with the key's 1357.44, 3845.2, 8328.1.

The stated t₀, dimensionless variables, and transcription of the three traction boundary conditions are correct; their implementation and the diffusion sign are not.

## 3. Method

**Hybrid finite-difference / polynomial singular-integral collocation, followed by least-squares wedge solution and tip fitting; general_or_wrong as delivered.** `solve.py` integrates derivatives of global Lagrange polynomials on uniform nodes, combines that matrix with a three-point Laplacian, keeps positive nearly real eigenvalues, and extrapolates them against 1/n. For K it switches to interior Gauss nodes, imposes a zero-mean opening, solves by least squares, and fits a square-root coefficient on s > 0.85.

Singularity subtraction and numerical eigenmode solution are appropriate ideas. This implementation is not a valid solution of the stated three-boundary-condition eigenproblem: it changes the elastic normalization, omits f″(1) = 0, treats the tip endpoint inconsistently, and discards problematic modes without establishing their irrelevance. It is neither the key's converged singular-integral method nor its crack-shape approximation.

For the submitted w and τ, the given conservation law yields wτ = −g″, hence gτ = A₁D²g and λf = −∫S f‴. The code's rate matrix −ED² has this decay sign, despite the opposite evolution stated in the explanation. Thus the textual sign error is distinguished from the numerical defects below.

## 4. Errors

Counts refer to distinct defects listed below, not repeated occurrences or each downstream wrong number. File references are inside `submission/` unless stated otherwise. Errors in unused exploratory files are classified as side errors because `solve.py` does not import or execute them.

| ID | Class | Quote or location and error |
|---|---|---|
| D1 | decisive | `solve.py:118`: **“(2.0 / np.pi)”** multiplies elasticity, although the defined w makes its coefficient exactly 1 (`ANSWER.md:15`). The mpmath reference at `solve.py:150` repeats this same erroneous factor. This changes eigenvalues and the equilibrium opening; it is not a physical validation of normalization. |
| D2 | decisive | `solve.py:127`: only the Neumann reflection is imposed at the interface. No constraint enforces g″(1) = 0. The returned slow mode has curvature about −2.16, so it is not an admissible eigenfunction. This also permits nonzero opening rate at the blocked interface. |
| D3 | decisive | `solve.py:114`: the principal-value logarithmic contribution is simply set to zero at a = 1. This finite-part endpoint prescription is not the interior-limit traction required by the stated problem. For a smooth opening with nonzero endpoint derivative, the interior logarithmic term diverges rather than tending to zero; the physical tip needs its singular endpoint behavior treated consistently. |
| D4 | decisive | `ANSWER.md:60` calls the excluded mode an **“over-damped gauge direction”**; `solve.py:172` filters complex and negative rate eigenvalues. At n = 26 the negative rate is about −1.6011×10⁶, has a nonzero initial-state coefficient (≈0.038), and is not a zero constant-opening gauge. Under gτ = −Ag it grows. Other eigenvalues have imaginary parts up to about 4.49×10⁵. Filtering does not establish a physical relaxation spectrum. |
| D5 | decisive | `solve.py:185` / `ANSWER.md:108`: **“robust monotone ... extrapolation”** is unsupported. The finest-grid λ₁ jumps from 0.62736 to 0.85117 and is rejected by a data-dependent filter. Merely replacing `i*(1/n)` by `i/n` changes n = 30 to 0.51326. Global uniform-grid polynomial differentiation is unstable here; the asserted ±0.05 uncertainty and continuum rate are not established. |
| D6 | decisive | `solve.py:206`: the overdetermined least-squares result is treated as an exact solution of **“A1 w∞ = 1”** (`ANSWER.md:69`) without inspecting its residual. Max |Ew∞ − 1| ≈ 0.54625 contradicts the claimed cancellation of traction and affects the tip estimate. |
| D7 | decisive | `ANSWER.md:76`: **“u = (2K/E*)√(ε/2π)”** is wrong for the specified half-opening u: the coefficient is 4K/E*. With the stated w, K/(σ₀√(πh)) = cπ/√2, not cπ√2 (`solve.py:214`). This factor can also be derived directly by evaluating the given kernel just ahead of the tip. The trial conversions in `t_K.py:41` do not supply a valid normalization either. |
| D8 | decisive | `solve.py:209` / `ANSWER.md:109`: a fixed-window square-root fit is presented as a limiting K **“stable to the digits quoted”**, without convergence or fit-window validation. Using the same code, n = 30 gives 3.11614 for cutoff 0.85 and 3.18672 for cutoff 0.95; n = 60 at cutoff 0.85 gives 3.27923. These are materially different estimates, not five-digit agreement. |

| ID | Class | Quote or location and error |
|---|---|---|
| S1 | side | `ANSWER.md:17,24` / `solve.py:22,26`: wτ = **“g_ss”** and gτ = −A₁D²g have the wrong signs. The final rate extraction happens to use −ED², the correctly signed rate matrix, so this textual defect is not charged again as a numerical cause. The same sign mistake is implemented in `t_time.py:40`, which evolves the rate matrix as the generator, and in the opening evolution stated in `t_be.py:8`, `gb6.py:5`, and `gb7.py:7`. |
| S2 | side | `ANSWER.md:21`: **“w_t = g”** is not the diffusion relation, with either dimensional or dimensionless time. Eliminating w requires its time derivative to be proportional to −g″. |
| S3 | side | `ANSWER.md:26` / `solve.py:29`: **“symmetric kernel”** is false for S. For example S(0.25,0.75) = −2.25, while S(0.75,0.25) = 0.75. Self-adjointness of an appropriately integrated elastic operator requires its domain and endpoint conditions; it does not follow from this premise, nor does self-adjointness of two factors imply that their product is self-adjoint. Repeated symmetry claims in `gb5.py` and `gb7.py` do not verify it. |
| S4 | side | `ANSWER.md:89`: cₖ = ⟨fₖ,1⟩/⟨fₖ,fₖ⟩ assumes ordinary L² orthogonality, which is not established for this product operator and fails for the nonsymmetric discrete matrix used. Expansion coefficients require left eigenvectors, the appropriate weighted inner product, or solving the eigenvector system. No coefficients are computed by the final solver. |
| S5 | side | `ANSWER.md:74`: **“rigid blocked interface”** changes the stated physical description. The substrate has the same elastic constants and forms part of the elastic half-space. Blocking diffusion/sliding is not infinite elastic stiffness. The final code nevertheless uses the supplied elastic kernel. |
| S6 | side | `ANSWER.md:72,75`: the displayed w ≈ −0.433 + 0.701√(1−s) is claimed to give u(h−ε) ∝ √ε. A zero-mean gauge can be used for traction calculation, but this opening does not close at the tip without an additive shift to enforce u(h) = 0. `t_equil.py:20` likewise says it fixes w(0) = 0 although its code fixes the mean. This shift would not change the fitted square-root coefficient. |
| S7 | side | `results.json:7` / `solve.py:231` / `t_final.py:11`: validation **“to 1e-13”** is contradicted by the actual maximum error 7.39×10⁻¹². The ≈7×10⁻¹² claim in `ANSWER.md:42` is reproducible. |
| S8 | side | `gb.py:23`: moving the two source derivatives onto S yields Sξξ, not Sss. They differ by 16(s−2ξ)/(s+ξ)⁴. The accompanying argument about vanishing at infinity is inapplicable to this integral on [0,1]; boundary terms must be checked at the actual endpoints. This auxiliary operator is not the key's eigenproblem. |
| S9 | side | `gb.py:30,32`: **“uniform traction is a steady state”**, f′(0) = 0, and automatic interface derivative conditions are false for the stated domain. A nonzero uniform traction violates g(0) = 0; the eigenmode condition there is f(0) = 0. S and its derivative do not vanish at s = 1 for general ξ. A constant derivative-null vector outside the boundary-condition domain proves no physical steady mode. |
| S10 | side | `gb.py:67`: the regular diagonal is 1/(2b³), not 5/(4b³)−3/(4b⁴). At b = 0.25 the code gives −112 instead of 32. |
| S11 | side | `gb.py:129,135`: replacing a hypersingular diagonal by its regular tail does not implement a Hadamard finite-part integral. The ordinary elastic kernel is also sampled directly at coincident collocation/quadrature nodes, leaving infinities in `E_gu`. These matrices cannot serve as an independent singular-integral solution. |
| S12 | side | `gb2.py:37`: J₂ = 1/(1−a)−1/a is wrong, and the labeled ordinary integral across its pole diverges. Its finite part is −1/(1−a)−1/a; the coefficient of the subtracted linear term in the cubic denominator is instead +1/(1−a)+1/a. At a = 0.25 the code's −2.66667 is neither quantity. |
| S13 | side | `gb2.py:15,151`: subtracting only constant and linear Taylor terms leaves a simple principal-value pole in the cubic-kernel remainder, despite calling it **“regular”**. The code then applies ordinary Gauss quadrature rather than treating that remaining pole analytically. Together with S12 this invalidates the auxiliary finite-part method. |
| S14 | side | `derive.py:12,55`: the proposed Sss itself is incorrect; its difference from the actual second derivative is 12ξ(s−5ξ)/(s+ξ)⁵. The process log acknowledges that an earlier derivative was wrong, and the later expression in `gb.py:54` is correct; the obsolete wrong formula remains delivered. |
| S15 | side | `derive.py:13,58`: its purported algebraic expansion is also unequal to the already incorrect preceding expression, by 12ξ(s−ξ)(s+ξ−1)/(s+ξ)⁶. The printed symbolic difference is nonzero. |
| S16 | side | `derive.py:33`: evaluating antiderivative endpoints without resolving real logarithmic branches does not return the real finite-part integral. For n = 2 at a = 1/4 it returns 1.12361 + πi; the intended integral is real. |
| S17 | side | `gb5.py:10,17`: the claimed weak form drops a derivative. Testing ∂s²∫Sg′ requires either ∂sS with v′ or S with v″, not the displayed v′Sg′ form. Thus −DᵀWD does not represent that stated transport operator and is not automatically symmetric/negative semidefinite. |
| S18 | side | `gb5.py:44`: **“pa = Pq[i, :]”** evaluates the basis at the i-th quadrature point, not the collocation node a. The singularity-subtraction value must be φⱼ(a), which is the nodal identity row. |
| S19 | side | `gb7.py:14`: an interior-node Lagrange basis does not **“automatically”** vanish at s = 0, and dropping one basis function is not generally projection of its constant combination. No such boundary/gauge elimination is actually applied in the function. |
| S20 | side | `gb7.py:16,58`: the Galerkin form interchanges A₁D² with D²A₁ without establishing commutation, omits the appropriate quadrature inner product, and forcibly symmetrizes a collocation matrix. Those operations change the operator rather than merely correcting quadrature asymmetry. They supply no valid alternative calculation of the reported rates. |
| S21 | side | `gb7.py:50,61`: the documented generalized mass/stiffness construction never builds M and returns only K. Consequently it does not implement the promised M c′ = −K c / generalized eigenproblem. |
| S22 | side | `gb6.py:66`: the cell kernel integrals divided by h are used as the opening matrix, but no adjacent-opening difference is applied. The reconstruction requires (wⱼ₊₁−wⱼ)/h. The delivered M therefore does not annihilate constant openings and is not the elasticity map described. |
| S23 | side | `gb6.py:52,61`: the stated principal-value cell logarithm has reversed sign: it should be log((a−lo)/(hi−a)). The centered cells happen to make both versions zero; this is incorrect reasoning yielding the right diagonal logarithm. |
| S24 | side | `gb6.py:89`: on the stated cell-centered grid, Dirichlet reflection g₋₁ = −g₀ gives first diagonal −3, not −2, while zero-flux reflection at the other end gives last diagonal −1 and adjacent coefficient 1, not −2 and 2. These are not the claimed ghost-cell boundary conditions. |
| S25 | side | `gb6.py:93`: the documented return of diffusion and elasticity matrices is actually `(D2, None)`. The proposed finite-volume transport calculation is incomplete. |
| S26 | side | `t_conv.py:13`: reflecting both ends imposes Neumann conditions at the free surface, replacing the required g(0) = 0. Its convergence sequence tests a different problem. |
| S27 | side | `t_grad.py:28`, `t_grad2.py:20`, `t_be.py:32`: the reflected last Laplacian row has its sign reversed. It should be −2gₙ + 2gₙ₋₁, divided by the squared last spacing. These graded/time diagnostics are not independent confirmations of the final operator. |
| S28 | side | `t_id.py:32,42`: the smoothness-selected `best` candidate is calculated but never used; the code records the first positive eigenvalue instead. It does not perform the physical-mode selection it describes. |
| S29 | side | `t_expm.py:36`: assigning rejected eigenvalues zero does not **“kill”** their mode contributions; it makes them constant in time. Combined with the generator sign mistake, this produces a zero fitted slope for n = 12–28, rather than independent confirmation of the reported slow decay. |
| S30 | side | `t_equil.py:23,27`: it locates a σ ∼ K/√(h−z) stress singularity on the relaxed boundary behind the tip, where it also correctly states σgb = 0. The tensile stress-intensity field is ahead of the crack tip (z > h), while the square-root opening is behind it. |
| S31 | side | `t_time.py:42,47`: the time check does not examine integrator success before evaluating its dense interpolant through τ = 5. The solver actually fails at τ ≈ 0.00126127 with **“Required step size is less than spacing between numbers.”** Its enormous extrapolated values and fitted λ = −0.5461 do not validate relaxation. |

| ID | Class | Quote or location and error |
|---|---|---|
| P1 | presentation | `ANSWER.md:60`: **“smallest magnitude eigenvalue ... large negative stiff mode”** confuses magnitude with algebraic ordering; the large negative eigenvalue does not have the smallest magnitude. |
| P2 | presentation | `ANSWER.md:91`: **“surface traction g(1,τ)”** is interface traction. The free surface is s = 0, where traction is prescribed zero. |
| P3 | presentation | `ANSWER.md:54,95`: τrel = 1/λ₁ is dimensionless but is assigned units of t₀. Write τrel = 1.67 or trel = 1.67 t₀. |
| P4 | presentation | `ANSWER.md:83`: **“dimensionless intensity ≈ 3.12 σ₀√(πh)”** attaches dimensional units to a dimensionless ratio. |

**Totals: 8 decisive, 31 side, 4 presentation.**

## 5. Checks

The following counts distinct check families; repeated versions on different grids are not additional checks. Auxiliary scripts with explicit validation/monitoring claims are included even when the final reproduction does not run them. “Could fail” does not mean that a check passed or that it validates the full physical problem.

| Check | Tests the deciding quantity? | Could fail? Assessment |
|---|---|---|
| 1. Elasticity versus mpmath (`solve.py:132`, `ANSWER.md:41`) | Only local interior elasticity, not rates, interface BCs, or K | **Yes.** A quadrature/interpolation error could disagree with the smooth-function reference. Error 7.39×10⁻¹² reproduces. Both sides share the incorrect 2/π, and the test uses Gauss nodes, so it does not validate uniform-grid conditioning or the tip. |
| 2. Mesh/1/n convergence (`solve.py:165`, `t_dir.py`, `t_final.py`, `t_grad*.py`, `t_id.py`) | Yes, proposed rates and mode shape | **Yes.** The final λ₁ sequence actually fails at n = 30; rejecting that point is not proof of convergence. Graded versions also change the interface Laplacian sign. No error bounds support the claimed precision. |
| 3. Real/positive spectrum and discarded “gauge” (`ANSWER.md:56`, `solve.py:172`) | Yes in principle | **No as implemented.** Filtering to nearly real positive eigenvalues forces the displayed property. The full spectrum contradicts it, and there is no valid gauge projection test. |
| 4. Smoothness/constant overlap (`ANSWER.md:63`, `t_lam1.py`) | An eigenvector diagnostic, not its admissibility or rate accuracy | **Yes.** Roughness and overlap could reveal a localized/oscillatory candidate. A smooth vector can still violate f″(1) = 0; the ≈0.92 overlap is not proof of a physical eigenmode. `t_id.py` does not use its computed smoothness selector. |
| 5. Clean square-root tip / stable K (`ANSWER.md:70,109`, `solve.py:209`) | Intended to test K | **No as delivered.** A least-squares coefficient is always returned, without residual, window variation, refinement, or independent K comparison. The grader's variations disprove the asserted stability. |
| 6. Direct time integration as “independent confirmation” (`t_time.py:1`) | Yes, late decay slope | **Yes, and it fails.** The copied script prints g(0.5,0.2) ≈ 1.37×10¹¹⁰ and fitted λ₁ = −0.5461, not the claimed decay. Its generator has the wrong sign and it ignores integrator status. |
| 7. Backward-Euler opening monitor (`t_be.py:10`) | Intended to measure late decay | **Yes as a diagnostic.** A printed history/slope could disagree with relaxation. No successful output is supplied or invoked by `solve.py`; its diffusion sign, interface row, and unpinned surface in the evolution prevent accepting it as confirmation. |
| 8. Discrete exponential / fitted slope (`t_expm.py:1`) | Yes, intended late decay slope | **Yes, and it fails.** The rerun returns zero slope at n = 12–28 while printing rates 6.06×10⁴ to 2.18×10⁶. At n = 32 it prints rate 0.24141 and slope 0.22066. These do not confirm 0.59877; the eigenvalue substitution is invalid. |
| 9. Symbolic kernel differentiation / expansion (`PROCESS.md:2`, `derive.py:54`, `gb.py:50`) | A local operator ingredient | **Yes.** A symbolic difference can be nonzero; the delivered expansion check indeed is. The later Sss formula in `gb.py` is correct, but the earlier expressions are not, and this does not validate moving source derivatives to observation derivatives. |
| 10. Claimed high-precision finite-part primitives (`gb2.py:18`) | A local integration ingredient | **No delivered comparison can fail.** No mpmath check for J₂/J₃ is supplied; the assertion is unsupported, and J₂ is wrong. |
| 11. A@1 relative error ≈10⁻¹¹ (`PROCESS.md:2`) | Derivative-null preservation, not the physical relaxation mode | **Yes for arithmetic/derivative implementation.** It could expose failure to annihilate a constant. Any correctly differentiated constant passes regardless of the kernel, and the nonzero constant traction violates the surface condition, so it cannot establish a physical steady mode or the requested rates. |

**11 reported check families; 8 could fail in the limited senses identified above.** The two passing principle checks in `AUTO_CHECKS.json` are external checks, not additional solver validations. Eigenvalue ordering is enforced by sorting, and trel = 1/λ₁ is assigned directly, so neither is independent evidence of correctness.

## 6. Claims and reproduction

**Partly supported.** `AUTO_CHECKS.json` reports exit 0 and five of five scalar claims printed by the reproduction, with every key-based numerical check failing. Running the exact `python3 solve.py` command in `scratch/submission/` confirms that result: the entire parsed `results.json`, including the method string, matches the delivered file. `scratch/solve.stdout` records the output; no submission files were changed.

The numbers are reproducible even though incorrect. The ≈7×10⁻¹² elasticity error and displayed grid/tip coefficients reproduce. The 10⁻¹³ accuracy, real/nonnegative full spectrum, physical gauge exclusion, stable quoted digits, exact traction cancellation, and independent time confirmation do not. The method string also describes the opposite evolution sign from the rate matrix actually diagonalized. Thus scalar reproducibility is fully supported, while the broader numerical and validation claims are only partly supported.

Grader diagnostics are in `scratch/audit.py` and `scratch/audit.stdout`; auxiliary time outputs are in `scratch/t_time.stdout` and `scratch/t_expm.stdout`. `scratch/t_time_status.stdout` records the failed integration's status and actual stopping time.

## 7. Questions

**None listed (0).** There is no request for missing physical information. Diagnostic questions about whether an eigenvector is a gauge mode are numerical investigations, not missing-information questions to the problem setter. The uncertainty caveat is not such a question.

## 8. Source

**None named.** No source paper or author is identified in the submitted answer, results, process log, or code. `t_K.py:33` says **“by the source”**, but supplies no name or citation. The source named in `KEY.md` is the key's attribution, not the solver's.

## 9. Challenges and key assessment

**No challenges.** The submission does not argue that the problem, model, key, or published result is wrong or ambiguous. Its different numbers and remarks about stiffness/discretization are not a challenge with an articulated argument against the stated model. Consequently there are no challenge-based alternative grades.

No error in the key was established. Its diffusion/eigenproblem signs follow directly from the given conservation law; the delivered code instead omits an explicit boundary condition and misnormalizes the kernel. The crack half-opening conversion follows from the given kernel's field just ahead of the tip and contradicts the submission's factor, not the key. No web search, fetch, or external paper access was used.

```json
{"answers": {"a": "partly", "b": "partly", "c": "wrong"},
 "supporting_state": "wrong",
 "method": "hybrid finite-difference/polynomial PV collocation; least-squares wedge and tip fit",
 "method_class": "general_or_wrong",
 "errors": {"decisive": 8, "side": 31, "presentation": 4},
 "checks_reported": 11, "checks_could_fail": 8,
 "claims_supported": "partly",
 "questions_listed": 0,
 "source_named": "none",
 "challenges": []}
```
