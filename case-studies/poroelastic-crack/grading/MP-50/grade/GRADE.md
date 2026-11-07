# Submission grade

## 1. Answers

Use the **corrected** key (`KEY.md:11–33`), including its explicit partial-credit rule for the undrained outer intensity. Nonzero numerical answers have a 1% band; `K2_0` has absolute tolerance 0.01.

| Question | Judgment | Reason |
|---|---|---|
| (a) | **Partly right** | Both volume limits and the long-time intensity are right; `K1_0 = 1` is the reverse-order outer limit, deserving partial credit, while the drained-tip value is 0.807593 and the short-time coefficient is wrong. |
| (b) | **Partly right** | `K2_0 = 0` is right, but the two long-time answers should be negative, rather than zero: (K_2/\sqrt{\pi L}=-\eta), (GV_2/L^2=-\alpha\pi(1-2\nu)/2). |
| (c) | **Partly right** | Five of nine numerical quantities pass; the four discrepancies below exceed the key's bands. |

Here (\eta=\alpha(1-2\nu)/[2(1-\nu)]=0.33375).

| Quantity | Submission | Corrected key | Judgment |
|---|---:|---:|---|
| `K1_0` | 1 | 0.807593 | Partly: valid outer limit; 23.82% above the requested tip limit |
| `K1_inf` | 1 | 1 | Right |
| `V1_0` | 1.8849555922 | 1.8849555922 | Right |
| `V1_inf` | 2.5132741229 | 2.5132741229 | Right |
| `V1_short_coeff` | 1.8054066674 | 2.440220 | Wrong: 26.01% low |
| `p_faces_0` | 0.7490636704 | 0.7490636704 | Right |
| `K2_0` | 0 | 0 | Right |
| `K2_inf` | 0 | −0.33375 | Wrong |
| `V2_inf` | 0 | −0.8388052385 | Wrong |

## 2. Supporting state

**Right for the requested intermediate quantity.** Setting fluid content variation to zero gives (p=- (\nu_u-\nu)\sigma_{kk}/[\alpha(1-2\nu)]); the initial outer face trace is −2, giving the correct `p_faces_0`. This is the undrained value before diffusion, as requested, rather than the imposed face pressure zero at positive time. The volume limits and drained/undrained volume ratio are also correct. Incorrect auxiliary stress and pressure assumptions are recorded below.

## 3. Method

**Elastic traction-only limits plus Betti reciprocity and a one-dimensional diffusion layer; class: general_or_wrong.** The elastic volume formulas and the reciprocity identity are useful and correct. The submission extends an elastic traction theorem to a coupled poroelastic problem, omits the drained tip and parts of the leading volume correction, and replaces loading 2's long-time transient limit with an inadmissible steady-state argument. It never solves or checks those deciding effects. Correct `K1_inf` and `K2_0` therefore do not validate its claimed time-independent intensity theorem.

## 4. Errors

Counts concern distinct errors below; repeated versions in the answer, code, and JSON are counted once. Consequent wrong numbers are not counted again.

### Decisive — 6

1. **Tractions alone are claimed to determine the poroelastic intensity.** “A cracked infinite body is *simply connected*” and its divergence-free stress is “unique” (`submission/ANSWER.md:38–44`); “K1 … at every time” and “K2 = 0 at every time” (`:67–70`; `submission/solve.py:27–28`). The slit exterior is multiply connected, and equilibrium plus tractions does not enforce strain compatibility. With pore-pressure eigenstrain the compatibility equation is
   
   (\nabla^2[(1-\nu)\sigma_{kk}+\alpha(1-2\nu)p]=0).
   
   Thus pressure can alter total stress despite unchanged face tractions. The argument misses the drained tip for `K1_0`, wrongly forces `K2_inf = 0`, and supplies invalid reasoning for the two intensities that happen to be right. The related claim that the combined permeable crack must have classical elastic intensity 1 (`:69–70,141–144`) repeats this error.
2. **The straight-face diffusion profile is substituted for the entire pressure change.** “p^A − p_und = −p_und·erfc” (`submission/ANSWER.md:86–93`) omits the coupled redistribution outside that local layer. Its face contribution lacks the factor ((1-\nu_u)/(1-\nu)): the key's face term is 1.354055, whereas the submission obtains 1.805407. The reciprocity identity itself is not the error.
3. **The leading tip contribution to volume growth is omitted.** The approximation “σ^B_kk ≈ −2 in the thin layer” (`submission/ANSWER.md:89`) is not uniform at the tips. Their singular inner regions contribute another 1.086165 at the same (\sqrt t) order. Correcting the face factor alone would still leave an incomplete coefficient; the total is 2.440220 (`KEY.md:19`).
4. **Loading 2 is assigned a regular decaying steady harmonic pressure that does not exist in this two-dimensional exterior problem.** “the steady pore pressure is the harmonic field p = 1 on the crack, p → 0 at infinity” (`submission/ANSWER.md:101–102`). In elliptic coordinates the angular mean of a bounded harmonic solution is affine in the exterior radial coordinate; it cannot equal 1 at the slit and tend to 0 at infinity. Long-time and spatial-infinity limits do not commute. This invalid steady replacement drives the wrong loading-2 limits.
5. **Local harmonic strain compatibility is treated as a global zero-stress solution.** “∇²p = 0 ⇒ no incompatibility … produces no stress” (`submission/ANSWER.md:102–105`). Local compatibility is insufficient for a globally admissible displacement and the constraints inherited from the infinite transient body. The long-time solution instead has the confining contribution −η and an elastic crack correction, giving `K2_inf = −η`.
6. **Zero stress is equated to zero crack opening.** “harmonic p => no stress => no opening” (`submission/results.json:20`; `submission/solve.py:36`). Even a hypothetical zero-stress state has nonzero eigenstrain (e_{ij}=\alpha(1-2\nu)p\delta_{ij}/(2G)); opening must be established from displacement compatibility and boundary conditions. This additional unsupported inference yields `V2_inf = 0`, rather than closure.

### Side — 3

1. **Wrong square-root branch in the auxiliary reference stress.** `submission/solve.py:84–85` uses `np.sqrt(z*z-L*L)` without selecting the branch asymptotic to z. It breaks reflection symmetry and decay: at ((x,y)=(-100,0.1)) it gives −4.0001000072, while the correct branch gives +0.0001000072. Symmetric integration then spuriously cancels the nonconstant part, explaining the reported −3.999600 at several heights (`submission/ANSWER.md:133–135`). This corrupts validation, not the separately assigned answer numbers. The elliptic-coordinate identity check uses a different, correct branch representation.
2. **The identity check is falsely said to test its reciprocity prefactor.** “would be wrong for a wrong coefficient in (R)” (`submission/ANSWER.md:118–119`). `check_identity()` (`submission/solve.py:87–100`) calculates only a reference stress integral and compares it with (2\pi); the prefactor in (R) never enters that calculation. Changing that prefactor would not change this check.
3. **The discarded large-domain identity results are incorrectly diagnosed as overflow.** “larger mumax overflowed the floats” (`submission/ANSWER.md:117–118`). Rerunning `verify_identity.py` gives finite integrals 2299.37266410 at μmax = 40 and 1179117420968.4712 at μmax = 50, with no overflow warning. The loss of accuracy comes from cancellation in the stress expression amplified by the large Jacobian. Discarding these values was appropriate; the explanation is wrong.

### Presentation — 2

1. **The theorem's name is misspelled:** “Michel's theorem” (`submission/ANSWER.md:37`; also the code) should be Michell's theorem.
2. **Precision is overstated:** “the ratio is 1 to 6 digits” (`submission/ANSWER.md:120`), immediately after reporting 0.999994. Its discrepancy is about (6\times10^{-6}), not agreement to six decimal places.

## 5. Checks

Count the **three numbered check groups**, rather than each quadrature refinement or time sample. **Two could fail as limited computational checks; none establishes the disputed physical quantities.** No assertions are required for a printed discrepancy to be visible.

| Reported check | Tests the deciding quantity? | Could fail? |
|---|---|---|
| (1) Reference-stress squared integral vs (2\pi L^2) (`ANSWER.md:110–120`) | Supports an identity relevant to the initial volume; does not test the poroelastic intensity, full short-time coefficient, or reciprocity prefactor. | **Yes**, for reference-field magnitude, Jacobian, or quadrature errors. |
| (2) Numerical erfc integral vs the formula for A (`:122–139`) | Compares the claimed coefficient with the same assumed face-only pressure model; it cannot detect omitted redistribution/tips. The branch error makes the claimed exact-field validation spurious. | **Yes**, for some factors, integration errors, or disagreement between the programmed coefficient and integral. Agreement does not validate the physical model. |
| (3) (K_1+K_2=1) (`:141–144`) | No independently computed combined poroelastic solution is used. Superposition requires the sum, not that its value equal the elastic traction-only answer. | **No** as a physical validation; the asserted target is the disputed assumption. |

The automatic volume-ratio and Griffith-limit checks pass, but were performed by the grading system and are not additional solver-reported check groups.

## 6. Claims and reproduction

**Overall: partly supported; all nine principal numerical claims are reproduced.** `AUTO_CHECKS.json` records successful reproduction and nine printed quantities. A fresh copy in `scratch/submission/`, run with `python3 solve.py`, exits 0; its entire regenerated `results.json`, including formulas, equals the delivered JSON. All main pasted `[check]` lines also match. Reproduction therefore supports the implementation, including its wrong physics.

The auxiliary script reproduces the first two cited identity refinements. The claimed additional refined A = 1.805231 and face-integral outputs are not printed by either delivered script, and the grid for that refinement is not specified. A finer symmetric grid can reproduce approximately the stated face integrals using the erroneous branch; a correct branch instead gives about −3.602609 at y = 0.01, not −3.999600. The overflow claim is contradicted by the auxiliary rerun. These qualifications prevent full support for all claimed validations.

Evidence is saved in `scratch/submission/rerun_stdout.txt`, `identity_stdout.txt`, their exit-status files, and `scratch/audit_stdout.txt`.

## 7. Questions and missing information

**Five listed question/assumption bullets** (`submission/ANSWER.md:148–163`), counted as listed topics rather than five actual requests for data. Sign and stress convention are sensible things to state, but the problem already specifies total stress and opening-positive volume; they do not justify a traction-only intensity. Prescribed face pressure already establishes permeability, and it does not imply vanishing loading-2 volume. The extra initial volume/short-time loading-2 note and diffusivity note are supplementary statements, not missing information. No clarification was needed to compute the requested values. The actual order-of-limits ambiguity is not identified by the solver.

## 8. Source

The only named method attribution is **“Michel's theorem / Bueckner weight function”** (`submission/ANSWER.md:37–38`). This quote is preserved in the profile. No source paper for this crack problem, or attribution to Detournay and Cheng, is given.

## 9. Challenges and key issues

There is no explicit criticism of a named paper. The statement **“A non-zero K2 would break superposition”** (`submission/ANSWER.md:144`) is an implicit challenge to any nonzero loading-2 result. Its argument is that the combined permeable crack must have elastic net-traction intensity 1. It is **checkable**, by solving the combined coupled problem or evaluating its long-time limit, but is not well argued: pressure produces back-stress, and superposition gives (K_{\mathrm{combined}}=K_1+K_2), not a prescribed value 1. With the key's limits the sum is (1-\eta=0.66625), fully consistent with superposition. The challenge is rejected; its error is included once under decisive error 1. It supplies no accepted alternative grading and no answer is marked contested.

The statement's initial tip-limit ambiguity is real and already recognized in `KEY.md:24–33`. Against the key `K1_0` is partly right; against the reverse-order interpretation that single quantity is right. Question (a) remains partly right either way because A is wrong. The submission gives no computation challenging the corrected tip or volume coefficients.

There is also a **dimensional transcription error in the key**, independent of the submission: `KEY.md:61` writes the opening-to-intensity prefactor as (\sqrt{\pi G}/[\sqrt8(1-\nu)]). It must be (G\sqrt\pi/[\sqrt8(1-\nu)]): substituting the key's elastic near-tip opening otherwise leaves an erroneous factor (G^{-1/2}), whereas stress intensity scales as modulus times opening divided by square root of distance. The solver does not use this relation. Correcting that typo leaves all numerical judgments and letter grades unchanged. The key's historical heuristic 0.75 and face-only coefficient are explicitly superseded and were not used for grading.

```json
{"answers": {"a": "partly", "b": "partly", "c": "partly"},
 "supporting_state": "right",
 "method": "Elastic traction-only limits with Betti reciprocity and 1D diffusion",
 "method_class": "general_or_wrong",
 "errors": {"decisive": 6, "side": 3, "presentation": 2},
 "checks_reported": 3, "checks_could_fail": 2,
 "claims_supported": "partly",
 "questions_listed": 5,
 "source_named": "Michel's theorem / Bueckner weight function",
 "challenges": [{"claim": "A non-zero K2 would break superposition", "checkable": true}]}
```
