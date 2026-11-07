# Grade

## 1. Answers

The P374 correction in KEY.md:11–33 governs the grade, rather than the superseded paper values. Nonzero numerical answers are compared within 1%; K2_0 has the stated absolute tolerance of 0.01.

| Question | Grade | Reason |
|---|---|---|
| (a) | partly | Both volumes and K1_inf are right; K1_0 is the legitimate undrained outer limit but misses the drained-tip limit, and the leading volume-growth coefficient is wrong. |
| (b) | partly | K2_0 = 0 is right; the claimed stress-free steady solution is inadmissible, giving both the wrong K2_inf and the wrong sign and magnitude of V2_inf. |
| (c) | partly | Five of nine numbers pass the key's bands; K1_0 earns the key's explicit partial credit, and the other three discrepancies are wrong. |

| Quantity | Submission | Governing key | Assessment |
|---|---:|---:|---|
| K1_0 | 1 | 0.807593 | partly: reverse order of limits only |
| K1_inf | 1 | 1 | right |
| V1_0 | 1.8849555922 | 1.8849555922 | right |
| V1_inf | 2.5132741229 | 2.5132741229 | right |
| V1_short_coeff | 1.1283791671 | 2.440220 | wrong |
| p_faces_0 | 0.7490636704 | 0.7490636704 | right |
| K2_0 | 0 | 0 | right |
| K2_inf | 0 | −0.333750 | wrong |
| V2_inf | +0.3399549584 | −0.8388052385 | wrong |

The correct general long-time loading-2 results are K2_inf = −η√(πL) and V2_inf = −απ(1−2ν)L²/(2G), where η = α(1−2ν)/[2(1−ν)]. For K1_0, the key distinguishes taking the tip limit first (0.807593) from taking the undrained time limit first (1). The submission uses only the latter and does not identify the ambiguity. This earns partial credit, rather than a contested grade.

## 2. Supporting state

**Partly right.** The requested pre-diffusion face pressure is derived correctly from ζ = 0 and the in-plane trace σ_kk = −2. It is an outer undrained limit; it does not replace the prescribed p = 0 on the drained faces at positive times. The drained/undrained elliptical openings, their integrals, and V1_inf/V1_0 = 4/3 are correct.

The symbolic fixed-e_xx skin derivative α(1−2ν)/[2G(1−ν)] and its diffusivity agree with the given equations. The erfc integral is also right. These local results do not justify freezing the outer displacement response or omitting leading tip effects. The proposed loading-2 pressure and globally stress-free displacement are wrong, as detailed below.

## 3. Method

**Elastic analogies + 1-D skin + stress-free steady ansatz; general_or_wrong.** Elastic analogies correctly supply the volume endpoints and K1_inf, but applying the undrained elastic stress field to the positive-time crack-tip intensity misses the required limit. The short-time calculation uses a valid local consolidation equation, then substitutes a local displacement integral for the complete coupled volume response. The key's diffusion/reciprocity argument and its leading tip correction are absent. The loading-2 argument checks local harmonicity without establishing the full pressure boundary conditions or a single-valued displacement. These are incorrect simplifications; a full transient numerical solver was not required.

## 4. Errors

Counts below treat repeated appearances of the same error in the answer, JSON, and code as one error. The incorrect numbers caused by an identified derivation are not counted again as separate transcription errors.

**Decisive: 5.**

1. **The undrained outer intensity is used as the positive-time tip intensity.** “every field with ζ = 0” and “K of a Griffith crack is independent of elastic constants” (submission/ANSWER.md:6–10,31; submission/solve.py:13–18). Elastic independence is true for the outer problem, but does not eliminate the diffusive tip layer. The value 1 is partly correct under the reverse order of limits; it is not the governing drained-tip value 0.807593.
2. **The face contribution to short-time volume is undercounted.** The frozen-outer-field skin bookkeeping gives A = 8(ν_u−ν)/[(1−ν)√π] (submission/ANSWER.md:22–25,99–101; submission/check2.py:19–45). Keeping normal traction fixed does not establish an unchanged outer displacement field at order √t. Even the key's face-only contribution is 16(1−ν_u)(ν_u−ν)/[(1−ν)√π] = 1.354055; the submitted 1.128379 is only 5/6 of that. Correct local algebra under the imposed constraint does not validate this replacement of the coupled response. The same bookkeeping also produces the unrequested loading-2 coefficient −1.506386 (ANSWER.md:42; solve.py:116–118), rather than even the face-only reciprocity value −1.807663.
3. **Leading tip effects are discarded.** “Tip zones contribute only O((ct/L²)^{3/4})” (submission/ANSWER.md:24,97–101). This misses their effect on the full volume response: the corrected key has an additional 1.086165 multiplying √(ct/L²). Thus repairing the face coefficient alone would still leave the requested leading coefficient wrong.
4. **The proposed loading-2 pressure violates the prescribed pressure boundary condition.** “p_ss = (2/π)Re[arcsin(L/z)] … =1 on the slit” (submission/ANSWER.md:13; submission/solve.py:85–99). The implemented principal arcsine gives p → −1 on both faces for −L < x < 0, and p = 0 approaching the crack center normally. For example, x = −0.6, y = ±10⁻⁹ gives approximately −1, not +1. Only the positive half was checked. Moreover, a bounded two-dimensional harmonic field cannot have constant +1 on the whole finite slit and decay to zero at infinity: under the exterior conformal map to a disk, its angular mean would need to be both 1 at the boundary and 0 at infinity. The long-time local limit is not the assumed global decaying steady field. This invalidates the claimed loading-2 endpoint calculation.
5. **Local strain compatibility is mistaken for global displacement compatibility.** “the strain is compatible because p_ss is harmonic” and “tip closure u_y(±L)=0” (submission/ANSWER.md:14–18; submission/solve.py:40–46,101–113). Harmonicity supplies a local compatibility condition, but the proposed primitive has ∮f′(z)dz = 4iL around the slit. With s = α(1−2ν)/(2G), the displacement u = s(Re f, Im f) changes by 4isL = 1.068i for the supplied constants. It is therefore not single-valued. On the upper face, u_y(L)−u_y(−L) = −0.534, so one integration constant cannot enforce closure at both tips. Doubling the positive-half displacement integral imposes symmetry that this construction does not possess. K2_inf = 0 and positive V2_inf are not admissible solutions. The extra claims that K2 “returns to 0” and V2 approaches +0.3400 (ANSWER.md:41–43) repeat these wrong endpoints; no transient history is computed.

**Side: 4.**

1. **False constraint independence.** “coefficient independent of the lateral strain constraint” (submission/ANSWER.md:23–24). Holding e_xx fixed gives the stated derivative divided by (1−ν); holding σ_xx fixed gives α(1−2ν)/(2G), without that divisor. The delivered algebra uses the fixed-e_xx case, so this extra assertion does not change its numbers.
2. **Inconsistent relative-order estimate.** The claimed tip term of order τ^{3/4}, where τ = ct/L², divided by the leading τ^{1/2} term would scale as τ^{1/4}, not the stated “ratio ~ 2.4·√(ct)/L” = 2.4τ^{1/2} (submission/ANSWER.md:97–98). Both estimates fail to justify the omission; the numerical omission itself is counted above.
3. **Wrong square-root branch in the elastic check on the left side.** submission/solve.py:60–61 uses the principal square root instead of the branch asymptotic to z. At z = −100L + 10⁻⁹iL it gives Re Z ≈ −2.000050, rather than a stress vanishing at infinity. The reported K extraction and face samples use the right side, so their printed numbers survive this defect.
4. **The quantity called Skempton B is not the conventional three-dimensional coefficient.** “Skempton consistency: B = 3(ν_u−ν)/(α(1−2ν)) = 1.12359551” (submission/ANSWER.md:70; submission/solve.py:128–130). This is a rescaling of the in-plane trace used in the problem. For the conventional hydrostatic Skempton definition, plane strain gives σ_kk^(3D) = (1+ν_u)σ_kk^(2D) in the undrained state, hence B = 3(ν_u−ν)/[α(1−2ν)(1+ν_u)] = 0.8025682183. The submitted face pressure remains correct because its custom coefficient is consistently paired with the two-dimensional trace.

**Presentation: 3.**

1. **Incorrect normalization label.** “du_face/(sqrt(c t))” still contains √t (submission/ANSWER.md:80; submission/check2.py:39–40). The actual coefficient per √(ct) has no √t. The final A calculation cancels this bookkeeping factor, so the label does not create an additional numerical error.
2. **Sign typo in the elastic extraction comment.** “constant term −tau in Z” (submission/solve.py:66) should be tau for Z = tau(1−z/√(z²−L²)). The implemented subtraction uses the correct −1 constant.
3. **Overbroad reproduction wording.** “prints all values above and the checks below” after python3 solve.py (submission/ANSWER.md:45–50) does not include the symbolic check. That requires the separately supplied python3 check2.py command (ANSWER.md:71–72). Both commands work.

## 5. Checks

There are **8 distinct reported check groups**: the seven numbered in ANSWER.md:54–85, plus the additional loading-2 short-time “consistency” item labelled [5] in solve.py:116–118. A check counts as capable of failure if its numerical comparison, residual, or independent equation-based algebra can expose an error in the subclaim it actually tests. This does not imply that it tests the complete physical solution. The scripts print comparisons rather than enforce pass/fail tolerances.

| Reported check | What it tests; limitation | Could fail? |
|---|---|---|
| 1. Westergaard stresses and K extraction | Can expose extraction errors in the elastic K1_inf/outer K1_0 and face trace. It does not test the positive-time poroelastic tip limit, and its left branch is wrong. | yes |
| 2. COD volume quadrature | Integrates the assumed ellipse and compares it with the closed-form volume; relevant to the correct volume endpoints, though it cannot establish the assumed elastic analogy independently. | yes |
| 3. Steady-pressure Laplacian, face values, far values | Residuals and boundary samples could reject a bad harmonic candidate. Positive-half-only face sampling misses the actual negative-half and center failures; it does not test global displacement compatibility. | yes |
| 4. Loading-2 displacement/volume integrals | Can expose quadrature or integration-factor errors in the constructed field. Agreement cannot establish its symmetry, single-valuedness, or tip closure. | yes |
| 5. erfc skin integral | Tests the one-dimensional diffusion integral. It does not test outer mechanical feedback or the leading tip contribution to A. | yes |
| 6. B/face-pressure consistency | Merely defines B from the same parameter combination and reevaluates 2B/3. It adds no independent comparison or residual. | no |
| 7. Symbolic constitutive re-derivation | Derives p*, the constrained strain derivative, and storage from the supplied laws, then compares symbolic formulas. It could catch hand-algebra errors, but reuses the disputed skin and outer-field assumptions for A. | yes |
| Extra code [5]. Loading-2 short-time “consistency” | Prints a coefficient and an assumed transient story; no reference, residual, or alternative calculation is compared. | no |

**6 checks could fail on their limited subclaims.** None tests the three deciding omissions: the drained-tip intensity, the complete leading short-time volume response, or the admissibility of the loading-2 steady displacement.

## 6. Claims and reproduction

**Supported as numerical/computational provenance: yes.** AUTO_CHECKS.json reports successful execution and all 9 numerical claims printed. Running python3 solve.py in scratch/submission reproduced the delivered results.json exactly, including its formulas. All printed numerical checks quoted from that script agree with scratch/solve_output.txt. Running python3 check2.py separately also reproduced its quoted symbolic identities and numbers, in scratch/check2_output.txt. There are no numerical output-versus-JSON mismatches. The single-command wording and normalization-label defects are recorded above.

This reproducibility does not certify physical correctness. The four failed key comparisons in AUTO_CHECKS.json are confirmed, and both automatic principle checks pass. Claims that the checks verify “every equation and boundary condition” or a globally exact loading-2 solution are contradicted by the boundary and circulation diagnostics, rather than by an inability to reproduce the calculations. The grader's executable evidence is scratch/grade_diagnostics.py and scratch/grade_diagnostics.txt.

## 7. Questions and missing information

**0 questions listed.** The seven bullets under “Questions / assumptions” (submission/ANSWER.md:87–104) contain conventions, modeling assumptions, and a limitation statement, rather than questions or requests for missing data. The unit load, total-stress definition of K, symmetric face loading, and opening sign are sensible. The tip-order estimate and frozen outer response are not justified. The acknowledgement that A lacks a full transient check is candid; a full transient run was not necessary, but the decisive asymptotic assumptions needed verification.

## 8. Source

The answer names an elastic method attribution: **“from the Westergaard function”** (submission/ANSWER.md:10). It also uses the Biot, Griffith, and Skempton eponyms. No source paper for this poroelastic crack problem is identified, and Detournay and Cheng are not named. The JSON records the actual quoted attribution, with this limited scope.

## 9. Challenges and key issues

**No explicit challenges.** The solver does not argue that the problem, a source paper, or a competing published result is wrong or ambiguous. Its stress-free ansatz is an alternative derivation, not a stated challenge. It is checkable, but the pressure-boundary and displacement-circulation checks reject it; this earns no contested exception. The initial-K order-of-limits ambiguity is already recognized by the key, and the submission receives the prescribed partial credit without having raised that ambiguity itself.

One dimensional typo in the key is worth separating from its numerical corrections: KEY.md:61 writes √(πG) multiplying lim D_n/√ξ. Since the latter has units √length, the factor must be **G√π**, giving K = G√π/[√8(1−ν)] lim D_n/√ξ, with units stress·√length. This typo is in the historical heuristic discussion and does not alter any governing corrected benchmark or this submission's grade. Grading with the literal printed factor or the dimensionally corrected factor yields the same (a), (b), and (c) classifications here; the evaluated submission uses G = 1. No additional substantive error in the corrected key is established by this review.

```json
{"answers": {"a": "partly", "b": "partly", "c": "partly"},
 "supporting_state": "partly",
 "method": "Elastic analogies + 1-D skin + stress-free steady ansatz", "method_class": "general_or_wrong",
 "errors": {"decisive": 5, "side": 4, "presentation": 3},
 "checks_reported": 8, "checks_could_fail": 6,
 "claims_supported": "yes",
 "questions_listed": 0,
 "source_named": "from the Westergaard function",
 "challenges": []}
```
