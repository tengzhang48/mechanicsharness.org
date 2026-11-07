# Grading profile

## 1. Answers

Use the corrected P374 values in KEY.md: `K1_0 = 0.807593` and `V1_short_coeff = 2.440220`, rather than the superseded paper values still printed later in the key. The bands are 1% for nonzero quantities and absolute 0.01 for `K2_0`.

- **(a) Partly right.** Both limiting volumes and the long-time K are right; `K1_0 = 1` is the undrained outer limit, earning the key's explicit partial credit, while the short-time volume coefficient is wrong.
- **(b) Partly right.** `K2_0 = 0` is right; both long-time answers incorrectly vanish instead of giving negative K and crack closure.
- **(c) Partly right.** Five of nine reported numbers pass the bands; four fail, with the initial loading-1 K additionally receiving the order-of-limits partial credit.

All values below use the requested normalized units.

| Quantity | Submission | Corrected key | Judgment |
|---|---:|---:|---|
| K1_0 | 1 | 0.807593 | Partly: alternative outer limit; outside the key band |
| K1_inf | 1 | 1 | Right |
| V1_0 | 1.884955592 | 1.884955592 | Right |
| V1_inf | 2.513274123 | 2.513274123 | Right |
| V1_short_coeff | 1.128379167 | 2.440220 | Wrong |
| p_faces_0 | 0.749063670 | 0.749063670 | Right |
| K2_0 | 0 | 0 | Right |
| K2_inf | 0 | -0.33375 | Wrong |
| V2_inf | 0 | -0.838805239 | Wrong |

For loading 2, the required general long-time expressions are `K = -eta*sqrt(pi*L)` and `V = -alpha*pi*(1-2*nu)*L^2/(2*G)`, where `eta = alpha*(1-2*nu)/(2*(1-nu))`. The submitted zero formulas fail as well as their numerical evaluations.

## 2. Supporting state

**Right for the requested supporting quantity.** The induced, pre-diffusion pressure is correctly derived from `zeta = 0` and the outer face trace `sigma_kk = -2`: `p_faces_0 = 2*(nu_u-nu)/(alpha*(1-2*nu))`. The distinction between this outer jump value and the imposed face value zero is correctly stated (ANSWER.md:75–79, 166–167).

The initial and final loading-1 ellipse normalizations are correct. The inverted constitutive law, storage coefficient, transport law (T), and harmonic stress-pressure combination in ANSWER.md:33–43 are also algebraically correct. The local one-dimensional strain and erf calculations are valid under their stated frozen-tangential-strain approximation; that approximation does not determine the complete finite-crack volume change. The incorrect loading-2 pressure limit is recorded separately below.

## 3. Method

**Elastic-limit substitution with local one-dimensional diffusion; general_or_wrong.** The solver uses Griffith/Sneddon elastic formulas with undrained or drained Poisson ratio, then integrates a face boundary-layer strain. SciPy BDF supplies a numerical check of the assumed scalar diffusion equation, not a solution of the coupled crack problem.

The elastic volumes and local diffusion are useful building blocks. The universal total-traction rule for K is invalid in a poroelastic body. The short-time calculation omits the reciprocal/global elastic response and the tip contribution required by the corrected key. This is an incomplete and physically incorrect route, rather than a validated alternative or a full transient solution.

## 4. Errors

Counts below refer to distinct erroneous statements or approximations, with repeated formulas and their numerical consequences grouped together. They are errors in the submission; the separate key typo is not charged to it.

| ID | Class | Evidence and consequence |
|---|---|---|
| D1 | Decisive | ANSWER.md:46–51: pore pressure supposedly “does not change the `r^{-1/2}` amplitude,” so K is fixed by face traction alone. This drops the coupled stress contribution and the drained tip layer. It produces `K1_0 = 1` instead of the tip-first 0.807593 and also underlies the loading-2 zero-K claim. The same premise appears in solve.py:12–14 and ANSWER.md:144–145. The value 1 remains valid for the reverse order of limits, but the solver does not identify that distinction. |
| D2 | Decisive | ANSWER.md:109–111: continuously imposed loading-2 pressure “at any fixed material point → 0,” hence zero strain, K, and V. A maintained Dirichlet value `p=1` does not drain away like a finite pulse. This loses the long-time back-stress and gives `K2_inf = V2_inf = 0` (solve.py:61–62). In two dimensions the local large-time pressure limit and the limit at spatial infinity cannot be interchanged. |
| D3 | Decisive | ANSWER.md:91–98, solve.py:54: integrating only the local vertical-strain loss gives `A = 8*(nu_u-nu)/((1-nu)*sqrt(pi)) = 1.128379167`. This is insufficient even for the face contribution: the key's face-only diffusion/reciprocity result is `16*(1-nu_u)*(nu_u-nu)/((1-nu)*sqrt(pi)) = 1.354055001`. The submitted value is 5/6 of that contribution; a correct scalar pressure integral does not supply the missing elastic response. |
| D4 | Decisive | ANSWER.md:155–158 explicitly neglects “near-tip modulation” while describing the result as the correct leading square-root term. The tips contribute at that same order, adding approximately 1.086165 to the face contribution. Omitting them prevents the requested leading coefficient 2.440220, independently of the face-prefactor error D3. The caveat acknowledges incompleteness but does not repair the reported asymptotic coefficient. |
| S1 | Side | ANSWER.md:113–116 calls the combined loading a consistency check and a “classical Rice–Cleary result,” claiming final opening `pi*(1-nu)` and K equal to 1 at all times. Linear superposition of the correct loading-1 and loading-2 limits instead gives normalized `K_inf = 1-eta = 0.66625` and `V_inf = pi*(1-nu)*(1-eta) = 1.674468884`. This additional combined-loading claim is not a requested answer. |
| S2 | Side | ANSWER.md:163–165 gives incorrect answers even for its proposed effective-stress interpretation. That interpretation prescribes total traction `-alpha` while retaining `p=1`, so it is `alpha*loading_1 + loading_2`. Its tip-first initial K would be 0.71875777, final K 0.55625, and final V 1.398008731, rather than 0.89, 0.89, and 2.237. The outer initial K alone can be 0.89; the claimed equality of the limits and omission of loading-2 closure are wrong. |

**Totals: 4 decisive, 2 side, 0 presentation.** The requested normalization and signs are stated clearly; the wrong signs/zero closure are physical errors, not unit or wording errors.

## 5. Checks

There are **four reported checks**: the three numbered checks and the earlier combined-loading consistency check. The boundary value and integral are grouped as components of the single reported boundary-layer check.

| Reported check | Does it test the deciding quantity? | Could fail? |
|---|---|---|
| Elastic COD, near-tip relation, and volume integral (ANSWER.md:120–124) | Tests elastic opening/volume normalization correctly; does not establish poroelastic initial K. The displayed elastic algebra is correct. | Yes: independently mismatched COD or K normalization would disagree. |
| Numerical boundary-layer integral (ANSWER.md:126–142; solve.py:66–91) | Tests the local pressure integral and its `2/sqrt(pi)` factor; does not test the complete coefficient A or tip contribution. | Yes: numerical and analytic integrals are independently obtained and could disagree. The face value alone cannot fail because it is imposed and reset in code. |
| `K1_0 = K1_inf` (ANSWER.md:144–145) | Asserts the deciding equality using the same erroneous traction premise. Both are hard-coded to 1; there is no independent physical comparison. | No. |
| Combined-loading “Consistency” (ANSWER.md:113–116) | Uses the same elastic premise and the assumed zero loading-2 response to construct both sides. It never checks an independently derived coupled response. | No. |

**Two checks could fail.** Their restricted scope explains why they pass despite the wrong crack answers. AUTO_CHECKS.json's additional volume-ratio and Griffith-limit checks pass, but those are grader checks, not extra solver-reported checks.

## 6. Claims and reproduction

**Claims supported computationally: yes.** AUTO_CHECKS.json reports a successful rerun and all nine scalar claims printed. An independent rerun of `python3 solve.py` in `scratch/reproduction/` also exits 0; its generated JSON equals the submitted JSON, including all formula strings. It reproduces the reported numerical pressure-integral check:

- numerical integral: -0.7979468936430866;
- analytic integral: -0.7978845608028655;
- relative difference: 7.81226298680798e-05.

There are no numerical reproduction mismatches. Only the scalar boundary-layer comparison is implemented numerically; the other checks are analytical statements or assertions in ANSWER.md. Reproducing hard-coded formulas does not validate their mechanics. `scratch/review_checks.py` and `scratch/review_checks.txt` retain the independent grading calculations.

## 7. Questions

**Two questions are listed** (ANSWER.md:171–175).

1. Total versus effective normal stress for loading 2: unnecessary as missing information, because TASK_PROBLEM.md:8 and 20 explicitly define total stress and `sigma_n = sigma_yy`. A traction-free face can still experience a pore-pressure-induced deformation; redefining the stress does not resolve the solver's coupling error.
2. Sign/meaning of `p_faces_0`: distinguishing the pre-diffusion outer pressure from the prescribed zero face trace is sensible. The prompt already requests the adjacent pre-diffusion value, however, and the solver's positive answer is correct.

The separate transport and surface-approximation caveats are discussed below rather than counted as additional listed questions.

## 8. Source

**An author pair is named:** “the classical Rice–Cleary result” (ANSWER.md:115); submission/PROCESS.md also says “Rice-Cleary type.” No paper title, year, or verifiable citation is supplied. The combined-loading formula attributed to those authors fails the superposition check above. Griffith and Irwin are also named as elastic building blocks. The submission does not identify the key's Detournay–Cheng paper.

## 9. Challenges and key review

**Loading-2 interpretation.** ANSWER.md:160–165 and 171–173 suggests an ambiguity between total and effective normal stress and claims that the literal total-stress reading gives zero K and V at both limits. Its arguments are the traction-only K rule and pressure draining to zero. This is checkable by solving the stated coupled boundary-value problem, but not well argued: the stress convention is explicit and both physical premises fail. Against the key, (b) and (c) are partly right. If the boundary condition were changed to the proposed effective-stress loading, the delivered three zero loading-2 values would make (b) wrong, while (c) would remain partly right through its loading-1 results; the paragraph's separate alternative estimates also fail the tip-first and final limits. No contested credit is warranted.

A simple check of the pressure-limit premise uses exterior slit elliptic coordinates `xi = Re(acosh((x+iy)/L))`. With pressure 1 on the slit and 0 on a confocal outer ellipse `xi_R`, the harmonic pressure is `p = 1-xi/xi_R`. At a fixed exterior point it tends to 1 as the outer boundary recedes, not 0. This is a reproducible large-domain limit illustrating the nonuniformity of the long-time/infinity limits, not a claimed full transient simulation.

**Transport-law uncertainty.** ANSWER.md:149–154 suggests an alternative “intended” unabsorbed diffusion coefficient that would rescale A by `sqrt((1-nu_u)/(1-nu))`. This is algebraically checkable. However, its own correct transport equation and local compatibility give `c`: substituting `d(sigma_kk)/dt = -alpha*(1-2*nu)/(1-nu)*dp/dt` into (T) yields `(1-nu_u)/(1-nu)*dp/dt = c*(1-nu_u)/(1-nu)*d2p/dy2`. There is no missing choice of transport law in the supplied equations. This tentative ambiguity is rejected; it neither changes the grade nor supplies a valid alternative coefficient.

**Surface-only caveat.** ANSWER.md:155–159 acknowledges that its prefactor is not an exact finite-crack result. That is a useful limitation on its own approximation, rather than a challenge to the problem or key. It qualitatively agrees with the key's need for a tip contribution but supplies no complete competing solution. It does not make (a) or (c) contested. The solver never raises the actual noncommuting limits for initial K.

**Grader correction to the key.** KEY.md:61 prints the opening-to-K prefactor as `sqrt(pi*G)/(sqrt(8)*(1-nu))`. For dimensional opening and distance it must be `G*sqrt(pi)/(sqrt(8)*(1-nu))`, as follows from the correct Irwin relation used by the submission. This is checkable: with L=1, unit opening traction, and G=4, substituting the key's own ellipse gives K=0.886226925 under the printed relation, whereas the corrected relation gives the required `sqrt(pi)=1.772453851`. Credit the submission's correct elastic relation. Grading with the literal versus corrected prefactor leaves all nine normalized numerical judgments and the (a), (b), (c) labels unchanged at the reproduction's G=1. This dimensional typo is distinct from the P374 numerical corrections, which govern this grade.

```json
{"answers": {"a": "partly", "b": "partly", "c": "partly"},
 "supporting_state": "right",
 "method": "Elastic-limit substitution with local one-dimensional diffusion", "method_class": "general_or_wrong",
 "errors": {"decisive": 4, "side": 2, "presentation": 0},
 "checks_reported": 4, "checks_could_fail": 2,
 "claims_supported": "yes",
 "questions_listed": 2,
 "source_named": "the classical Rice–Cleary result",
 "challenges": [{"claim": "Loading 2 may mean effective rather than total normal stress", "checkable": true}, {"claim": "An intended unabsorbed transport law could change the short-time coefficient", "checkable": true}, {"claim": "Grader: the key's COD-to-K prefactor has sqrt(G) instead of G", "checkable": true}]}
```
