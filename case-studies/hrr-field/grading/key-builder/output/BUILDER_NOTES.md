# Builder notes

These probes use the mechanics of corpus P015, Hutchinson's 1968 tensile-crack analysis, to test transfer rather than recall. Work used only the supplied text layer and page images; no network source was consulted. [PROCESS.md](PROCESS.md) records steps with actual UTC `date` output. Exploration is in `scratch/`; each final reference script runs independently of it.

## Deliverables and intended use

Each of `T1/`, `T2/`, and `T3/` contains a short `TASK_PROBLEM.md`, a grader-facing `REFERENCE.md`, a 100-point `RUBRIC.md`, standalone `reference.py`, and generated `reference_results.json`. Give the solver **only its problem statement**, alongside whichever experimental knowledge condition is being tested. References, rubrics, these notes, scripts, and the paper are grader/build artifacts. Solve and score the probes independently.

| Probe | What it tests | New element | Requested result |
| --- | --- | --- | --- |
| T1 | Nonlinear singular-field scaling, angular solution, energy calibration, and consistent normalization | Noninteger n = 6.5, not a tabulated parameter | p, normalized tensile coefficient H, disk-energy coefficient C |
| T2 | Adapting equilibrium/compatibility and the inner constitutive limit to a different geometry | 120° opening, n = 4, plane-strain notch | Stress exponent p and same-radius normal-stress ratio R, with model justification |
| T3 | Recognizing loading-history limits and separating local-field existence from elastic outer calibration | Three engineering requests with scale/history data | Use / replace / conditional recommendations and usable next steps |

Problem statements do not name the source, its authors, or its equation numbers. Numerical problems provide the constitutive law, geometry, units, shear convention, and requested normalization. The general Airy and kinematic identities are mathematical aids; the T2 statement does not provide its reduced equation, select the elastic/plastic limit, give the eigenvalue, or supply the crack exponent. The solver must make those adaptations. T3 supplies evidence rather than predetermined classifications. The statements are roughly 316, 282, and 467 whitespace-delimited words, respectively; their intended computation/writing budget is 45 minutes each.

## T1: parameter transfer

The paper's main power-law analysis makes the inner stress scale as r⁻¹⁄⁽ⁿ⁺¹⁾ and ties its amplitude to a path-independent energy flux. Merely knowing that exponent does not give the angular tensile coefficient or the energy inside a disk. T1 forces an actual angular solution and a consistent J calculation at a noninteger parameter.

The probe deliberately stipulates the **fully plastic limiting law** rather than assigning a physical radius and extrapolating it through an elastic transition or blunted core. E is only part of the reference strain. H and C are dimensionless, so α, σ₀, and ε₀ cancel after calibration. The unknown eigenvalue is solved numerically and then checked against its exact scaling value. This avoids pretending that an arbitrarily chosen finite-distance stress is a verified approximation to the full elastic-plastic problem.

Reference normalization f(0) = 1 is arbitrary; it is not max effective stress = 1. For the latter convention, the additional diagnostic I is approximately 4.8183790924. H = 1.8935697252 and C = 0.4573057840 are invariant under that change. C is U/(JR), with U **per unit thickness**; it is not the paper's energy coefficient Eₙ, which in the source's small-scale-yielding nondimensionalization equals πC. H is an opening normal-stress coefficient, not an effective-stress coefficient.

Checks include collocation versus shooting, Cartesian numerical differentiation of the energy flux on three radii, and exact incompressible linear elasticity at n = 1. The source's rounded n = 3 plane-strain table is only an extra corroboration. The new reference is not interpolated from it.

## T2: geometry adaptation

The appendix demonstrates how a crack angular eigenproblem can be adapted to a V-notch. Its illustrated example has a 90° opening, n = 3, plane strain, and a rounded exponent 0.225. The probe instead uses a 120° opening and n = 4. The **solid angular interval is 240°**, so the free faces lie at ±120°, not ±60°. This angle convention is often the most consequential implementation error.

The complete elastic-plus-power-law relation in the statement requires the solver to choose its inner asymptotic limit. The leading stresses are singular; at n = 4 the nonlinear strains dominate their linear counterparts. Geometry then selects the exponent. Substituting the crack's 1/(n+1) = 0.2 gives the wrong answer; using the elastic-wedge exponent ≈ 0.384269 also gives the wrong answer. The reference p ≈ 0.1561752603 is obtained from a new nonlinear eigenproblem.

An absolute local amplitude would require global loading and matching, so the task asks only for two amplitude-independent quantities. The reference uses a mixed stress/strain collocation system and a separately formulated fourth-order shooting system. Both recover an analytic symmetric elastic-wedge solution when n = 1, and the zero-opening crack limit when the face angle becomes π. They also reproduce the appendix's rounded special case. This checks the new result without treating the source's example as the answer.

The physical sharp-tip singularity would be cut off by a finite notch radius, finite strain, or damage. The problem explicitly asks for the continuum asymptote; an unsolicited demand for those unprovided inputs should not prevent solving the requested exponent and ratio.

## T3: applicability and productive scope recognition

The paper assumes small strain and a stationary sharp crack, neglects damage and geometric nonlinearities, and uses total-deformation plasticity subject to no unloading and nearly proportional pointwise stress histories. It determines the inner exponent locally. An elastic outer-field calibration of the amplitude additionally requires small-scale yielding. These assumptions are distinct.

* **A** supplies the local-history evidence and scale separation needed for a useful plane-strain hardening field and elastic K-to-J calibration. It allows approximation corrections; refusing because a metal is plastic or a tip is singular misses the intended useful regime.
* **B** breaks the stationary, nearly proportional, non-unloading description through reversal, tip advance, and wake contact. A single monotonic curve cannot supply cyclic internal variables. The correct response adapts the analysis and identifies missing cyclic/contact inputs; it need not refuse all stress estimation.
* **C** removes reliable small-scale elastic matching while retaining a possible local energy-calibrated field. The observed J agreement is encouraging, but local history, shape/dominance, and finite-strain/blunting margins are unresolved. A thoughtful answer may restrict the asymptote to a verified annulus and use a full specimen solution at the target. It should not equate large-scale yielding with impossibility of every local hardening asymptote, or treat similar J values as sufficient proof of its validity.

The applicability reference is intentionally qualitative. Code independently audits its length and energy scales and constructs a history counterexample using two routes. It does not solve a three-dimensional cyclic crack or establish an inner annulus for C. This limitation is stated in its reference rather than hidden behind a computed number. The rubric caps an unqualified incorrect main recommendation equally for unjustified application and unnecessary refusal; qualified, productive adaptations remain eligible for full credit.

## Reproduction and grading cautions

From this folder, run the three scripts:

```bash
python T1/reference.py
python T2/reference.py
python T3/reference.py
```

The first two need NumPy/SciPy; T3 uses the Python standard library. They write their result JSON next to their own source and print their checks. No script reads `paper.txt`, a page image, another probe, or `scratch/`. Numeric tolerances in statements and rubrics are intentionally much looser than the achieved internal numerical agreement. Very close numerical agreement does not establish physical accuracy outside the stipulated idealization.

The supplied text layer has damaged equation typography. Signs, the plane-strain effective-stress factor 3/4, and the appendix's notch exponent were verified against the page images. The plane-strain n = 3 table corroborates the normalization; historical table values are rounded and should not be required to match many digits. Do not award extra points for matching source terminology or recycling its code. Score transferred mechanics, adaptation, verification, and scope recognition using each rubric's explicit partial-credit rules.
