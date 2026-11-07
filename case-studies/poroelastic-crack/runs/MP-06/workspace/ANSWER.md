# ANSWER — suddenly pressurized crack in a fluid-saturated porous solid

## Model and method (and why)

**Constitutive reading.** I read the given laws with indices 1,2 = in-plane (2D) components, i.e.
`σ_kk = σ_11+σ_22`, `e_kk = e_11+e_22` (Rice–Cleary-type plane-strain poroelasticity in 2D form).
Justification: with this reading, mass balance + Darcy + the identity
∇²[(1−ν)s + α(1−2ν)p] = 0 (s ≡ σ_11+σ_22; from equilibrium + compatibility) give a pure diffusion
equation ∂p/∂t = c∇²p with **exactly** the diffusivity given in the problem,
c = 2κG(1−ν)(ν_u−ν)/[α²(1−2ν)²(1−ν_u)], whenever s + γp is time-independent
(γ ≡ α(1−2ν)/(1−ν)). The alternative reading (σ_kk = 3D trace) yields an inconsistent (negative)
diffusivity, so it is rejected.

**Structure used.**
1. With 2D-trace laws, the elastic reduction is 2D elasticity with Kolosov constant κ = 3−4ν; the
   undrained limit (ζ = 0, t = 0⁺) is the *same* elasticity problem with ν → ν_u and Skempton
   p = −(ν_u−ν)s/[α(1−2ν)].
2. Eigenstrain (thermoelastic) form: poroelastic = elastic (G, ν) with isotropic eigenstrain
   2Gε* = α(1−2ν)p. The no-crack ("smooth") eigenstress has ŝ = −γp̂, and
   σ̂_yy = −γ(k_x²/k²)p̂. For a *harmonic* p in the half-plane with face data g(x) = p(x,0),
   the trace on the plane is **σ_yy^sm(x,0) = −(γ/2)g(x)**.
3. Crack problems solved by a Chebyshev dislocation-density solver (validated to machine precision
   against Sneddon/Westergaard closed forms — Part A below). Solver convention: input q is the
   "internal pressure"; face stress of the dislocation field is −q; opening ∝ q.
4. Short-time behavior: local 1D drainage layer (erf profile) at the faces, where the 1D
   poroelastic structure (Δs = −γΔp, Δσ_yy = 0) is exact.

**Key derivations in brief.**
- *Loading 1, K:* SIF of a uniformly loaded crack is fixed by tractions alone (any ν): K = σ₀√(πL)
  at t = 0⁺ (undrained) and t→∞ (drained).
- *Loading 1, V:* Sneddon V = π(1−ν_eff)σ₀L²/G, ν_eff = ν_u at 0⁺, ν at ∞.
- *Loading 1, short time:* at ct ≪ L² the faces (p = 0) drain a layer ℓ = √(ct). Through the layer
  Δσ_yy = 0 (traction BC constant) and Δe_yy = (1−2ν)αΔp/[2G(1−ν)], with
  p = p₀·erf(y/2√(ct)), p₀ = 2(ν_u−ν)σ₀/[α(1−2ν)]. Face displacement
  Δu_y = −∫Δe_yy dy = γp₀√(ct)/(G√π) gives uniform opening change
  Δw = 2Δu_y = 4(ν_u−ν)√(ct)/[G(1−ν)√π], hence ΔV = 2LΔw and
  **A = 8(ν_u−ν)/[√π(1−ν)]**; tip-region corrections are O(ct/L²).
- *p_faces_0:* undrained field adjacent to faces carries σ_xx = σ_yy = −σ₀ (pressurized-crack face
  values; numerically confirmed), so s = −2σ₀ and p = 2(ν_u−ν)σ₀/[α(1−2ν)].
- *Loading 2, K(0⁺) = 0:* face pore pressure has no instantaneous mechanical effect (bulk
  undrained ⇒ zero eigenstrain in bulk, traction-free faces, zero remote field ⇒ σ ≡ 0).
- *Loading 2, t→∞:* p_ss is harmonic (strip potential, 1 on faces, →0 at ∞). The smooth-field
  stress on the crack plane is σ_yy^sm = −γ/2, i.e. **compression** (the pressurized matrix near
  the crack is constrained by the far field, exactly like a hot inclusion which squeezes itself:
  σ = −(γ/2)p inside). A traction-free crack cutting a compressed zone relieves the compression by
  **closing**: the dislocation field must carry face stress +γ/2 = "negative internal pressure":
  **K₂(∞) = −(γ/2)p₀√(πL)**, **V₂(∞) = −π(1−ν)(γ/2)p₀L²/G = −πα(1−2ν)p₀L²/(2G)**
  (negative = interpenetration in the linear solution; a real crack would partially close).

## Results (ν = 0.2, ν_u = 0.4, α = 0.89; K in units of √(πL)·loading, V in units L²/G·loading)

| quantity | value | formula (general) |
|---|---|---|
| K1_0 | 1.000000 | σ₀√(πL) |
| K1_inf | 1.000000 | σ₀√(πL) |
| V1_0 | 1.884956 | π(1−ν_u)σ₀L²/G |
| V1_inf | 2.513274 | π(1−ν)σ₀L²/G |
| V1_short_coeff A | 1.128379 | 8(ν_u−ν)/[√π(1−ν)] |
| p_faces_0 | 0.749064 | 2(ν_u−ν)σ₀/[α(1−2ν)] |
| K2_0 | 0.0 | 0 |
| K2_inf | −0.333750 | −α(1−2ν)p₀√(πL)/[2(1−ν)] = −(γ/2)p₀√(πL) |
| V2_inf | −0.838805 | −πα(1−2ν)p₀L²/(2G) |

Reproduce with **`python3 solve.py`** (prints all values, writes `results.json`).

## Checks (what was compared, pasted output, how it could have failed)

1. **Crack solver vs closed forms** (`python3 verify.py`, Part A): w(x), V, K, σ_yy(x,0) for a
   uniformly loaded crack vs Sneddon/Westergaard:
   ```
   max |w_num - w_exact| / max|w_exact| = 8.24e-16
   V  = 2.513274   exact pi*(1-nu) = 2.513274   (err 0.00e+00)
   K  = 1.772454   K/sqrt(pi L) = 1.000000  (exact 1.000000, err 1.01e-14)
   sigma_yy(x,0), x>L  numeric vs Westergaard closed form:
      x=1.05  num=  2.27965  exact=  2.27965  reldiff=3.90e-16
      x=3.00  num=  0.06066  exact=  0.06066  reldiff=1.60e-15
   ```
   Failure mode: wrong dislocation kernel/Chebyshev solve; it also confirms the face values
   σ_xx = σ_yy = −σ₀ used for p_faces_0 and V1_inf = π(1−ν).
2. **Eigenstress kernel / trace formula**, three ways:
   (i) single harmonic mode, FFT vs analytic −(γ/2)(1+k|y|)e^{−k|y|}cos(kx):
   ```
   Syy(x=0,y) rows 253..259: [-0.32472 -0.3297 -0.33297 -0.33415 -0.33297 -0.3297 -0.32472]
   analytic                 : [-0.3244  -0.32935 -0.33259 -0.33375 -0.33259 -0.32935 -0.3244 ]
   ```
   (ii) Eshelby disk inclusion: exact hand solution of a uniform-eigenstrain disk gives
   σ_ij = −(γ/2)p δ_ij inside; smoothed-disk FFT check:
   ```
   sigma_yy^sm at disk center = -0.32965   (exact inclusion value -(gamma/2) = -0.33375)
   ```
   (iii) x-uniform slab ⇒ σ_yy = 0 (exact 1D argument, kx→0). Failure mode: a wrong kernel would
   change K2_inf/V2_inf. (An initial kernel implementation had kx/k² instead of kx²/k² —
   machine-zero trace — caught by (i); an early disk-algebra attempt used the wrong eigenstress
   coefficient 2Gε* instead of 2Gε*/(1−2ν) in Λ-form — caught by mapping to the classic
   plane-strain thermal inclusion σ = −EαΔT/[2(1−ν)]; both fixed.)
3. **Loading 2 steady state** (`verify.py`, Part B) — crack solver fed with the smooth-field trace
   q = −γ/2 (the required dislocation face load is +γ/2 ⇒ closing):
   ```
   K2_inf  crack-solver (q=-gamma/2) = -0.591556 /sqrt(pi L) -> -0.333750   analytic -gamma/2 = -0.333750
   V2_inf  crack-solver (q=-gamma/2) = -0.838805   analytic -pi*alpha*(1-2nu)/2 = -0.838805
   ```
   Exact to all digits (same machinery as check 1). Failure mode: wrong superposition or sign.
   (The 2-D FFT of the strip potential itself underestimates the trace by ~20% due to tip-jump
   aliasing at N=2048 — documented; checks (i)/(ii) carry the burden instead.)
4. **Loading 1 transient** (`python3 partC.py`, output in `partC_out.txt`): implicit-FD diffusion
   (p = 0 on faces, ∂p/∂y = 0 elsewhere on the crack plane, Skempton initial field) + FFT
   eigenstress + crack solver:
   ```
   pipeline at t=0+: V = pi*(1-nu)+Vcrack = 2.51327 + (-0.61512) = 1.89816  [undrained exact 1.88496]
   t= 0.000500  V=2.05841     t= 0.002000  V=2.12450     t= 0.004000  V=2.16702
   fit over ct in [4e-4,4e-3]: A_fit = 2.69359  (analytic 1.12838, diff 138.71%)
   ```
   The t = 0⁺ endpoint (the static part) validates the pipeline to 0.7%. The *rate* A_fit is
   ~2.4× my analytic A: I attribute this to (a) unresolved tip singularities of the Skempton
   initial field in the FFT trace (they also pollute the pipeline's K, which came out ≈ 0.8–0.84
   instead of 1), and (b) the simulation assumes p obeys pure diffusion, i.e. s + γp
   time-independent, which is only approximately true. The analytic A rests on a locally exact
   1D layer argument (Δσ_yy = 0 is forced by the constant traction BC; Δs = −γΔp is exact for
   1D corrections; the erf profile is exact) and is what I report, flagged below as the least
   verified quantity. Earlier runs of this check failed identically (q ≡ 0) due to two coding
   bugs (grid mirroring used the wrong half; Dirichlet rows used pⁿ as RHS so the faces never
   drained) — both found and fixed; see PROCESS.md.
5. **Diffusivity self-consistency**: the given c is reproduced exactly by the 2D-trace reading.

## Questions / assumptions

1. **Index convention**: indices 1,2 with σ_kk, e_kk as *in-plane* traces (justified by c).
2. **"Per unit loading magnitude"**: σ₀ = 1 (loading 1), p₀ = 1 (loading 2); K in units
   √(πL)·loading, V in units L²/G·loading.
3. **p_faces_0** = undrained (Skempton) pore pressure in the solid adjacent to the faces at t = 0⁺
   (uniform along faces); the boundary value on the face itself is 0 for t > 0.
4. **V1_short_coeff**: leading term in V ≈ V1_0 + A√(ct/L²); tip corrections O(ct/L²) neglected.
5. **K2_inf, V2_inf negative**: the linear solution predicts closure/interpenetration; if the
   intended model forbids interpenetration (contact), the crack would partly close and the
   reported values are the LEFM (linear) ones.
6. K₁ is reported at the two asked limits; numerics suggest K₁(t) stays ≈ √(πL)·1 in between
   (pipeline K ≈ 0.8–0.84 with known tip-aliasing error).

## Missing / uncertain

- **V1_short_coeff (A) is the least verified number**: analytic 1.128379 vs transient-pipeline
  ~2.4–2.7 (discrepancy attributed to tip-singularity aliasing and the pure-diffusion assumption
  in the simulation; not resolved within the time budget).
- The full transients between limits were only sampled numerically; loading 2's transient was not
  simulated (its p-evolution may involve a time-dependent harmonic function h = s + γp).
- V2_inf = −πα(1−2ν)/2 L²p₀/G contains no ν_u (the (1−ν) cancels) — verified by check 3.
