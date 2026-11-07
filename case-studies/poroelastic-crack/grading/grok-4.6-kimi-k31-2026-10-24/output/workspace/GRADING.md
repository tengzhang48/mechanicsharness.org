# Grading report: suddenly pressurized poroelastic crack

Reviewer: Grok 4.6 (Azure / OpenCode). The two submissions were read independently from `A/` and `B/` against `PROBLEM.md` and the corrected entries in `REFERENCE.md`. Historical paper numbers are recorded only as superseded literature values. Independent arithmetic is in `scratch/`. `A/` and `B/` were not modified.

Corrected numerical key (`ν = 0.2`, `ν_u = 0.4`, `α = 0.89`; `K` in units of `√(πL)`, `V` in units of `L²/G`):

| Quantity | Corrected key | Historical paper (not the key) |
|---|---:|---:|
| `K1_0` (tip first: `r → 0` at fixed `t > 0`, then `t → 0⁺`) | **0.807593** | 0.75 (heuristic, Eq. 29) |
| `K1_0` (undrained first / outer intensity) | **1** (conditional reading of the original wording) | — |
| `K1_inf` | 1 | 1 |
| `V1_0` | `π(1 − ν_u)` = 1.884956 | same |
| `V1_inf` | `π(1 − ν)` = 2.513274 | same |
| `V1_short_coeff` | **2.440220** = 1.354055 (faces) + 1.086165 (tips) | 1.354055 (faces only, Eq. 18) |
| `p_faces_0` | `2(ν_u − ν)/(α(1 − 2ν))` = 0.749064 | same |
| `K2_0` | 0 | 0 |
| `K2_inf` | `−η` = **−0.33375**, `η = α(1 − 2ν)/[2(1 − ν)]` | same |
| `V2_inf` | `−απ(1 − 2ν)/2` = **−0.838805** | same |

Tolerance: 1% on nonzero values; absolute 0.01 on `K2_0`.

The original statement asked for `K` “just after loading (`t = 0⁺`)” without specifying the order of the crack-tip and time limits. An undrained-first value `K1_0 = 1` is therefore scored as **partial / conditional**, not as an unqualified numerical error. It is still a physical miss if the solver never distinguishes an undrained outer field from a locally drained tip.

---

## Submission A

### Overall verdict

**Partial.** Loading 1 volumes, long-time `K`, induced face pressure, and `K2_0` are correct and follow from the stated conditions. `K1_0 = 1` is the valid undrained-first outer intensity, but the write-up never separates that outer field from a drained tip and never produces the tip-first intensity. The short-time volume coefficient retains the complete face-layer (inner + outer) contribution and matches the paper’s incomplete faces-only formula; the tip contribution at the same order is absent. Loading 2 is replaced by a drained elastic crack loaded by the Biot effective face traction `α`, so `K2_inf` and `V2_inf` do not follow from the coupled problem.

### Boundary and limit selection

**Loading 1, face conditions.** The submission imposes the physical conditions of the statement: total normal traction `σ_n = −1` and face pore pressure `p = 0` for `t > 0` (`A/ANSWER.md` lines 46–65; `PROBLEM.md` loading 1). Instantaneous bulk response is taken with `ζ = 0`, which converts the constitutive law to undrained elasticity with `ν_e = ν_u` (lines 33–43). The induced pressure just inside the solid is obtained from `S = σ_xx + σ_yy = −2` on the faces,

```
p_faces(0⁺) = 2(ν_u − ν)/(α(1 − 2ν)),
```

while “the imposed boundary value at the drained face itself is zero for `t > 0`” (lines 59–65). That distinction is correct.

**Loading 1, order of limits for `K`.** The boxed result is

```
K₁(0⁺) = K₁(∞) = √(πL)
```

(`A/ANSWER.md` lines 51–52). Both endpoints are the ordinary Griffith intensity for unit total face traction. There is no drained-tip process zone, no Rice–Simons / Atkinson–Craster argument, and no discussion of `r → 0` versus `t → 0⁺`. The reported `K1_0` is therefore the **undrained outer intensity** (time limit first at fixed distance from the tip, or equivalently `K` read from the undrained Griffith field). Under the original ambiguous wording that number is conditionally admissible. The submission does **not** recognize the distinction and does **not** establish the tip-first value 0.807593.

**Loading 2.** The statement is a coupled poroelastic BVP: `σ_n = 0`, `p = 1` on the faces, zero initial data, fields vanishing at infinity. The submission does not solve that problem. It introduces the Biot effective stress `σ̄_ij = σ_ij + α p δ_ij` and replaces the pair of loadings by a drained elastic crack whose face effective traction is `σ̄_n = −1 + α` (`A/ANSWER.md` lines 197–214):

```
K_{1+2}(∞) = (1 − α)√(πL),
V_{1+2}(∞) = (1 − α) π(1 − ν)L²/G,
```

hence `K2_inf = −α` and `V2_inf = −α π(1 − ν)` in the requested units. Equilibrium is on **total** stress. The effective-stress divergence is `α ∇p`, which is not a homogeneous drained elasticity problem unless `p` is spatially uniform. Loading 2 holds `p = 1` on a finite crack and `p → 0` at infinity, so the long-time pressure is a nonuniform harmonic field and the back-stress coefficient is the plane-strain modulus `η = α(1 − 2ν)/[2(1 − ν)]`, not `α`. The substitution is a wrong reduced problem, not a mere omission of a formula: it is stated as the governing face condition for the long-time calculation (lines 199–214) and is re-imposed as a verification identity in `A/verify.py` lines 142–148 (`K1_inf + K2_inf = 1 − α`).

**Far-field / initial data.** Instantaneous `ζ = 0` in the bulk and vanishing far-field totals are used consistently for loading 1. For loading 2 at `t = 0⁺` the argument is that there is no total traction and no bulk deformation, so a zero-thickness pressure layer has no integrated effect (`A/ANSWER.md` lines 191–195). That is the right `K2_0 = 0`. At long time the far field of loading 2 is implicitly that of a drained Griffith crack with effective load `α`, which is the reduced problem above, not the coupled steady state.

### Question-by-question assessment

**(a) Loading 1**

- `K1_0`: **partial / conditional.** Numerically 1, which is the undrained-first outer value. The write-up equates it to `K1_inf` and never treats a locally drained tip. Not the paper heuristic 0.75, and not the corrected tip-first 0.807593.
- `K1_inf`: **correct.** Drained Griffith value 1, independent of material constants.
- `V1_0`: **correct.** Undrained Sneddon volume `π(1 − ν_u)`.
- `V1_inf`: **correct.** Drained Sneddon volume `π(1 − ν)`.
- `p_faces_0`: **correct.** Undrained Skempton evaluation on the solid side of the faces, with the face itself held at `p = 0`.
- `V1_short_coeff`: **incorrect (incomplete: faces only).** The 1-D face layer is derived carefully: `S + [α(1 − 2ν)/(1 − ν)] p` harmonic across the layer, storage coefficient matching the given diffusivity `c`, `p = p_f erf(y/(2√(ct)))` with drained face `p(0,t) = 0` (`A/ANSWER.md` lines 69–116). Volume change is converted by an undrained Betti identity `V₁ = V₁(0⁺) + Λ ∫_Ω ζ θ_u dA` (lines 117–136). Evaluating `θ_u` at the faces and integrating two layers of length `2L` produces

  ```
  A = 16(ν_u − ν)(1 − ν_u) / ((1 − ν)√π) = 1.354055.
  ```

  `A/verify.py` lines 137–139 split this into an inner (layer) piece `8(ν_u − ν)/((1 − ν)√π) = 1.128379` and an outer undrained piece `8(ν_u − ν)(1 − 2ν_u)/((1 − ν)√π) = 0.225676`. The text at line 187 is therefore accurate about what is retained: face-layer displacement **and** the accompanying undrained outer response to that layer. What is **not** retained is the tip region. The integral uses the face value of `θ_u` along the whole length `2L` (lines 137–167). There is no tip-layer analysis, and no claim that geometric smallness of the tips kills their contribution — the tips are simply omitted. Geometric smallness would not have justified the omission: the tip contribution 1.086165 is the same order in `√(ct/L²)`. The delivered coefficient is exactly the historical paper Eq. 18 (faces only), which the corrected key rejects.

**(b) Loading 2**

- `K2_0`: **correct.** Instantaneous pressure has not entered the bulk; zero total traction is compatible with a null elastic field.
- `K2_inf`: **incorrect.** Reported `−α = −0.89` from the effective-traction reduction. Corrected value `−η = −0.33375`. Consequence: magnitude too large by `α/η = 2(1 − ν)/(1 − 2ν) = 8/3`. Combined hydraulic-fracture intensity would be `1 − α = 0.11` instead of `1 − η = 0.66625`.
- `V2_inf`: **incorrect.** Reported `−α π(1 − ν) = −2.236814` from the same reduction. Corrected value `−απ(1 − 2ν)/2 = −0.838805` (equivalently `−η π(1 − ν)`). Same 8/3 overestimate of closure. The sign (closure) is right; the constitutive coefficient is not.

The long-time loading-2 intensity and volume do **not** follow from the coupled BVP. They follow from a homogeneous drained Griffith problem with face load `α`.

**(c) Numerical evaluation**

The JSON values match the formulas in `A/results.json` and `A/ANSWER.md` §4 exactly (`scratch/check_output.json`). Arithmetic is not the issue.

### Numerical comparison (corrected key)

| Quantity | Reported | Corrected reference | Agrees (1% / 0.01) | Qualification |
|---|---:|---:|---|---|
| `K1_0` | 1.0 | 0.807593 (tip first) | no | agrees with undrained-first value 1; partial/conditional under the original wording |
| `K1_inf` | 1.0 | 1.0 | yes | |
| `V1_0` | 1.884956 | 1.884956 | yes | |
| `V1_inf` | 2.513274 | 2.513274 | yes | |
| `V1_short_coeff` | 1.354055 | 2.440220 | no | equals paper faces-only 1.354055; missing tips 1.086165 |
| `p_faces_0` | 0.749064 | 0.749064 | yes | |
| `K2_0` | 0.0 | 0.0 | yes | |
| `K2_inf` | −0.89 | −0.33375 | no | used `α` in place of `η` |
| `V2_inf` | −2.236814 | −0.838805 | no | used `−α π(1 − ν)` in place of `−απ(1 − 2ν)/2` |

General formulas, relative to the stated normalizations: `V1_0`, `V1_inf`, `p_faces_0`, `K1_inf`, and `K2_0` match. `K1_0` matches the undrained-first normalization, not the tip-first key. `V1_short_coeff` matches the paper faces-only formula, not the corrected faces+tips formula. `K2_inf` and `V2_inf` are the effective-traction formulas, not the coupled drained formulas.

### Verification coverage

`A/verify.py` was copied to `<reviewer-scratch>/grade_A_verify/` and run there; it passed. It was not executed inside `A/`.

What the checks establish:

- Constitutive reduction `ζ = 0` ⇒ undrained Hooke law (`check_undrained_constitutive`).
- Ordinary elastic Griffith volume and a near-tip elastic `K` for total face traction 1 (`check_elastic_crack`). That elastic `K` cannot reject a drained-tip poroelastic intensity.
- Face Skempton pressure and `ζ = 0` on the undrained faces.
- 1-D layer diffusivity identical to the given `c`, heat-equation identity for the erf profile, drained-face value, and the `erfc` integral (`check_layer_diffusion`).
- `∫ S² dA = 2π L²` for the unit Griffith field, and the reciprocal identity reproducing `V1_inf − V1_0` (`check_trace_integral`). This checks the long-time undrained-to-drained volume jump, not the short-time ordering.
- Reconstruction of **their** short-time coefficient from the layer amplitude, including the inner/outer split (`check_short_time_coefficient`). This confirms arithmetic on the assumed face-layer formula. It cannot reject a missing tip contribution.
- Superposition `K1_inf + K2_inf = 1 − α` and the analogous volume identity (`check_loading_superposition`). These encode the reduced effective-traction problem; they would pass for a wrong coupled solution of exactly this type.
- JSON consistency with the same formulas.

Nothing in the suite can fail because of a wrong coupled loading-2 field, a wrong order of the `K` limits, or a wrong short-time tip ordering. The checks are identities, quadratures, and restatements of the candidate’s own reduced model.

### Material errors or unsupported claims

1. **Equating `K1_0` with `K1_inf`** (`A/ANSWER.md` lines 51–52). Supported only if `K` is the outer undrained/drained Griffith intensity. Unsupported as a statement about the physical tip singularity at `t = 0⁺`.
2. **Loading-2 long-time reduction to effective face traction `α`** (lines 197–225; `A/verify.py` 142–148). Demonstrably the wrong face/bulk problem: total-stress equilibrium plus nonuniform `p` does not yield a homogeneous drained crack with load `α`. Consequence: `K2_inf` and `V2_inf` off by 8/3.
3. **Short-time volume from face layers only** (lines 137–185). Internally consistent for that truncation and equal to paper Eq. 18. Incomplete at the retained order `√(ct/L²)` because tips contribute at that order. The write-up does not argue geometric smallness; it simply never treats the tips.
4. Line 261, “There are no unresolved uncertainties,” is overstated relative to the omitted tip limit and the coupled loading-2 field.

---

## Submission B

### Overall verdict

**Partial.** Same loading-1 endpoint volumes, `K1_inf`, face pressure, and `K2_0` as a correct undrained/drained elastic reading. Same conditional `K1_0 = 1` without a drained-tip analysis. The short-time coefficient retains only the local uniaxial-strain opening of the face layers (the “inner” piece of A’s faces formula) and dismisses tips by geometric smallness of an `O(√(ct))` interval — which does not control the order of the volume contribution. Loading 2 is the same effective-traction / globally drained reduction as in A, so `K2_inf` and `V2_inf` are wrong for the same reason.

### Boundary and limit selection

**Loading 1, face conditions.** Total traction `σ_n = −P` and hydraulically open faces `p = 0` for `t > 0` are used (`B/ANSWER.md` lines 23–46). Instantaneous bulk `ζ = 0` gives undrained elasticity. Face pressure just inside the solid is the same Skempton value, with an explicit remark that the imposed face value is zero and that the step creates a surface discontinuity later resolved by diffusion (lines 45–46). That is the correct loading-1 pair.

**Loading 1, order of limits for `K`.**

```
K₁(0⁺) = P √(πL),    K₁(∞) = P √(πL)
```

(lines 49–59). Again the undrained (then drained) outer Griffith intensity. No locally drained tip, no order-of-limits discussion. Reported `K1_0 = 1` is **undrained-first / outer**. Conditionally admissible under the original wording; the distinction is not recognized; the tip-first result is not established.

**Loading 2.** After a correct `K2_0 = 0` (lines 139–142), the long-time argument is (`B/ANSWER.md` lines 144–149):

> In the drained long-time limit, the appropriate traction on the solid is expressed through the effective stress `σ_ij + α p δ_ij`. On a crack face, `σ_n = 0` and `p = P`, so the effective normal stress is `+α P`. … equivalent to `P_o = −α P`. The pore pressure in the surrounding bulk drains toward zero at fixed spatial position in the long-time diffusive limit.

Two substitutions are explicit here: (i) replace total face traction by Biot effective traction `αP`; (ii) assume a **globally drained** bulk, `p → 0` at each fixed point. The coupled problem does not drain to `p = 0` near the crack while `p = 1` is held on the faces; the steady field is harmonic with that Dirichlet data. The resulting `K2_inf = −α` and `V2_inf = −α π(1 − ν)` are the same wrong reduced values as in A. This is a stated change of boundary-value problem, not an unsupported leap from a correct BVP.

### Question-by-question assessment

**(a) Loading 1**

- `K1_0`: **partial / conditional.** Value 1, undrained outer; no tip-first analysis.
- `K1_inf`: **correct.**
- `V1_0`, `V1_inf`: **correct.** Sneddon with `ν_u` then `ν`.
- `p_faces_0`: **correct.**
- `V1_short_coeff`: **incorrect.** The face layer is set up with the right ingredients for 1-D diffusion: `Δσ_yy = 0`, `Δe_xx = 0`, storage matching `c`, erf profile (`B/ANSWER.md` lines 61–93). Opening is then taken as the integral of local `Δe_yy = D Δp` through the layer with the far field unchanged (lines 100–113). That produces

  ```
  A = 8(ν_u − ν) / ((1 − ν)√π) = 1.128379.
  ```

  Terms actually retained: local uniaxial-strain extension of the two face layers, integrated over length `2L`. Terms **not** retained: (i) the undrained outer opening induced by the layer fluid-content change (A’s outer piece 0.225676, which together with this inner piece recovers the paper faces-only 1.354055); (ii) the tip contribution 1.086165. Line 134–135:

  > The neglected tip regions occupy only `O(√(ct))` of the crack length and do not change the displayed leading term.

  Geometric smallness of the tip interval does not establish that its volume is `o(√(ct))`. That claim is unsupported as an order argument and is false for this problem: the tip piece is `O(√(ct/L²))` and comparable to the faces. Even if tips were dropped, the faces contribution itself is incomplete relative to the paper faces formula.

**(b) Loading 2**

- `K2_0`: **correct.** Null bulk field at `t = 0⁺`.
- `K2_inf`, `V2_inf`: **incorrect**, same effective-traction / globally drained reduction as A (`B/ANSWER.md` lines 144–149; `B/verify.py` lines 103–104). Same numerical values `−0.89` and `−2.236814` versus `−0.33375` and `−0.838805`. Sign of closure is right.

**(c) Numerical evaluation**

JSON matches the stated formulas (`B/results.json`; `scratch/check_output.json`).

### Numerical comparison (corrected key)

| Quantity | Reported | Corrected reference | Agrees (1% / 0.01) | Qualification |
|---|---:|---:|---|---|
| `K1_0` | 1.0 | 0.807593 (tip first) | no | agrees with undrained-first value 1; partial/conditional under the original wording |
| `K1_inf` | 1.0 | 1.0 | yes | |
| `V1_0` | 1.884956 | 1.884956 | yes | |
| `V1_inf` | 2.513274 | 2.513274 | yes | |
| `V1_short_coeff` | 1.128379 | 2.440220 | no | inner face-layer only; not paper faces 1.354 nor faces+tips 2.440 |
| `p_faces_0` | 0.749064 | 0.749064 | yes | |
| `K2_0` | 0.0 | 0.0 | yes | |
| `K2_inf` | −0.89 | −0.33375 | no | used `α` in place of `η` |
| `V2_inf` | −2.236814 | −0.838805 | no | used `−α π(1 − ν)` |

Formulas in `B/results.json` are the normalized strings corresponding to those values. They agree with the requested normalizations where the physics is right (`V1_*` endpoints, `p_faces_0`, `K1_inf`, `K2_0`) and encode the reduced models where it is not.

### Verification coverage

A copy of `B/verify.py` was run in `<reviewer-scratch>/grade_B_verify/` (with a copy of `results.json`). It passed. It was **not** run in `B/`, because lines 106–108 write `results.json` in the current directory and would overwrite the submission.

What the checks establish:

- Face undrained identities: Skempton `p0`, `ζ = 0`, and `ν + α(1 − 2ν)β = ν_u`.
- Layer relations `Δσ_xx = −α(1 − 2ν)/(1 − ν) Δp`, `Δe_yy = D Δp`, and layer storage reproducing `c`.
- Riemann-sum `erfc` integral and elliptical crack area.
- Compact `A` reconstructed from `4 D p0` per unit length on two faces — i.e. a restatement of the local-strain formula.
- `V1_inf − V1_0 = π(ν_u − ν)` and `V2_inf = −α V1_inf`, `K2_inf = −α`.

These cannot reject a wrong coupled loading-2 solution (they assume `K2_inf = −α`), a drained-tip `K`, or a tip-scale short-time term. There is no check of the geometric-smallness claim in lines 134–135.

### Material errors or unsupported claims

1. **`K1_0 = K1_inf = 1`** without a tip/time-limit distinction (`B/ANSWER.md` lines 49–59). Same conditional status as A.
2. **Loading-2 effective traction plus globally drained bulk** (lines 144–149). Wrong reduced BVP; `K2_inf` and `V2_inf` do not follow from `σ_n = 0`, `p = 1`, vanishing far field.
3. **Short-time faces truncated to local `Δe_yy`** (lines 100–132). Misses the outer undrained response already present in a complete face-layer calculation.
4. **Tip regions negligible because they occupy `O(√(ct))` of the length** (lines 134–135). Unsupported as an order estimate; false for the retained `√t` term.
5. Line 170, no unresolved uncertainties, is overstated.

---

## Comparison of the two submissions

Do not read this as a model ranking. Both executions share the same physical skeleton and the same two structural reductions:

- Loading 1 is treated as an undrained-then-drained Griffith crack with the correct total-traction / drained-face pair. Both report `K1_0 = 1`, `K1_inf = 1`, `V1_0 = π(1 − ν_u)`, `V1_inf = π(1 − ν)`, `p_faces_0 = 2(ν_u − ν)/(α(1 − 2ν))`, `K2_0 = 0`. Neither analyzes a locally drained tip or the order of limits. Under the original wording, `K1_0 = 1` is partial/conditional for both, not an unqualified miss.
- Loading 2 is replaced, in both, by a drained elastic crack loaded by Biot effective face traction `α` (A via superposition of loadings 1+2 with `σ̄_n = −(1 − α)`; B via `σ̄_n = αP` and an explicit `p → 0` bulk). Both therefore report `K2_inf = −0.89` and `V2_inf = −2.236814` instead of `−η = −0.33375` and `−0.838805`. The error is in the conditions, not in the arithmetic.

They diverge on the short-time volume:

- A retains face-layer diffusion **and** the undrained outer response to that layer (reciprocal identity). Its `A = 1.354055` is exactly the paper’s faces-only coefficient. Tips are omitted without a geometric-smallness argument.
- B retains only local uniaxial-strain face opening, `A = 1.128379`, which is A’s inner piece alone, and asserts that tip intervals of length `O(√(ct))` cannot affect the leading term.

Relative to the corrected key, both miss the tip piece 1.086165; B also misses the outer face piece 0.225676. A is therefore closer to the historical faces-only result and farther from a mere local-strain integral, but neither has a complete `O(√t)` coefficient.

Verification in both cases confirms internal identities of the reduced models and cannot catch the coupled loading-2 error, the `K` limit order, or the short-time tip ordering. B’s script additionally overwrites `results.json` if run in place.

Numerical agreement count is not a substitute for the boundary/limit verdict: five quantities agree unconditionally for both (`K1_inf`, `V1_0`, `V1_inf`, `p_faces_0`, `K2_0`); `K1_0` agrees only with the undrained-first reading; `V1_short_coeff`, `K2_inf`, and `V2_inf` disagree for both, with different short-time numbers.

---

## Reference limitations and problem ambiguities

- The original `PROBLEM.md` does not specify whether `t = 0⁺` for `K` is tip-first or time-first. The corrected key takes tip-first `0.807593`. This review follows the supplied rubric: `1` is partial/conditional, `0.75` would have been wrong, and the distinction must be recognized to earn the tip-first value. That clarification is **not** treated as having been in the solver prompt.
- `REFERENCE.md` keeps Detournay & Cheng’s heuristic `K1_0 = 0.75` and faces-only `V1_short_coeff = 1.354055` on record. Numerical grading uses the corrected column (P374): `K1_0 = 0.807593`, `V1_short_coeff = 2.440220`. No defect was found in those corrected entries; independent evaluation of `η`, `V1_*`, `p_faces_0`, `K2_inf`, and `V2_inf` matches the printed digits (`scratch/check_reference_and_submissions.py`, `scratch/loading2_reduction_note.py`).
- The paper’s long-time loading-2 coefficient is `η`, not `α`. The identity `η π(1 − ν) = α π(1 − 2ν)/2` holds. Submissions that replace `η` by `α` are not matching a hidden paper formula; they are solving a different reduced elasticity problem.
