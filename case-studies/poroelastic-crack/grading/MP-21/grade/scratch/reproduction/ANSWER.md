# ANSWER — MP-21: suddenly pressurized crack in a fluid-saturated porous solid

Reproduce everything with: **`python3 solve.py`** (writes `results.json`, prints all values; full log in `run_output.txt`).

## Model and method

Biot poroelasticity exactly as given (Rice–Cleary plane-strain form; the two given constitutive laws are mutually consistent: ζ = 0 ⇒ 2G e_ij = σ_ij − ν_u σ_kk δ_ij, i.e. the undrained solid is the drained one with ν → ν_u, G unchanged). Key structure used:

- **t = 0⁺ (undrained):** ζ ≡ 0 everywhere; mechanical response = elasticity with (G, ν_u). For traction BVPs on the cracked plane the *stress* field is the Griffith field, independent of the elastic constants.
- **t → ∞ (drained/steady):** ∇²p = 0 with the face BCs; stress = elasticity with (G, ν) once p is uniform/harmonic-compatible.
- **Small times (loading 1):** diffusion boundary layer of thickness √(ct) ≪ L on each crack face; locally a 1D "oedometer" problem (e_xx kinematically constrained, σ_yy = −1 fixed by the face traction), solved analytically (the layer diffusivity turns out to be exactly c — see check 2).
- Classical Griffith-crack results: K = P₀√(πL); opening Δu_y = (2(1−ν)/G)P₀√(L²−x²); V = π(1−ν)P₀L²/G (plane strain); on the faces σ_xx = σ_yy = −P₀.

## Answers, (c): ν = 0.2, ν_u = 0.4, α = 0.89 (per unit loading magnitude)

| quantity | value | general formula |
|---|---|---|
| K1_0 | **1.0** [√(πL)] | dP·√(πL) |
| K1_inf | **1.0** [√(πL)] | dP·√(πL) |
| V1_0 | **1.88496** [L²/G] | π(1−ν_u)dP·L²/G |
| V1_inf | **2.51327** [L²/G] | π(1−ν)dP·L²/G |
| V1_short_coeff | **1.12838** | A = 8(ν_u−ν)/(√π(1−ν)), V ≈ V1_0 + A√(ct/L²)·L²/G |
| p_faces_0 | **0.74906** [dP] | 2(ν_u−ν)/(α(1−2ν))·dP |
| K2_0 | **0.0** | 0 |
| K2_inf | **0.0** | 0 |
| V2_inf | **0.0** | 0 |

Note V1_0 < V1_inf: under sustained face pressure the crack *opens* with time (undrained → drained creep). For the hydraulically pressurized crack (σ_n = −P, p = P = loading 1 + loading 2): K = P√(πL) at both limits; V(0⁺) = π(1−ν_u)PL²/G, V(∞) = π(1−ν)PL²/G.

## Derivation (brief)

**Loading 1 (σ_n = −1, p = 0 on faces).**
- *t = 0⁺:* undrained elastic Griffith solution; K = √(πL) (SIF of the Griffith crack depends only on face tractions); V = π(1−ν_u)L²/G. Induced pressure p_u = −ησ_kk with η = (ν_u−ν)/(α(1−2ν)); on the faces σ_xx = σ_yy = −1 ⇒ σ_kk = −2 ⇒ **p_faces_0 = 2η** (uniform along the faces, away from tips).
- *t → ∞:* steady p ≡ 0 (faces drained, far field zero ⇒ harmonic ⇒ 0); drained elastic Griffith solution; K = √(πL), V = π(1−ν)L²/G.
- *Small times:* on each face a layer of thickness ~√(ct) drains from 2η to 0. Locally 1D (uniform along x): e_xx fixed at its undrained value (2ν_u−1)/2G, σ_yy = −1. The constitutive + mass laws give σ_kk = const − α(1−2ν)p/(1−ν) and hence p_t = c∇²p exactly (the (1−ν_u)/(1−ν) diffusivity factor cancels the σ_kk-coupling). So p = 2η·erf(y/√(4ct)) and δe_yy = (ν_u−ν)(erf−1)/(G(1−ν)) ≤ 0: face displacement δu_y(0,t) = −∫δe_yy dy = 2(ν_u−ν)√(ct)/(G√π(1−ν)) (face moves away from the solid ⇒ opening increases). Both faces over length 2L: ΔV = 4L·δu_y ⇒ **A = 8(ν_u−ν)/(√π(1−ν))**. Tip-region contributions are O(ct/L²), subleading.

**Loading 2 (σ_n = 0, p = 1 on faces).**
- *t = 0⁺:* no mechanical load and no diffusion yet (p = 1 holds on a measure-zero set); σ ≡ 0 ⇒ **K2_0 = 0** (V2_0 = 0 as well).
- *t → ∞:* in 2D, a constant-pressure finite source pressurizes the whole plane pointwise: p → 1 everywhere (logarithmically slowly; check 3b), consistent with "vanishes at infinity" at every finite t. With p uniform, σ ≡ 0 satisfies equilibrium, face tractions, decay; strain e_ij = α(1−2ν)/(2G)·δ_ij is compatible. Hence the final state is uniform swelling u = β_sw(x,y), β_sw = α(1−2ν)/(2G), with **zero relative displacement across the crack: K2_inf = 0, V2_inf = 0**. Transiently the infiltration *closes* the crack (V2(t) < 0, K2(t) < 0): the swelled near-face matrix pushes the faces toward each other, δu_y ~ −α(1−2ν)√(ct)/(G√π(1−ν)) for ct ≪ L²; this is a transient only.

## Checks (what was compared, output, how it could have failed)

All from `python3 solve.py` (full log in `run_output.txt`).

1. **Elastic Griffith crack, two independent solutions.** Compared (i) Westergaard closed form Z = P₀(z/√(z²−L²)−1) [equilibrium residual by finite differences; σ_yy, σ_xx on the faces; K from σ_yy·√(2πr) ahead of the tip; face displacement from −∫e_yy dy] and (ii) an independent distributed-dislocation solve (climb-dislocation kernel σ_yy = σ_xx = (G/2π(1−ν))PV∫b/(x−ξ)dξ, b = φ/√(L²−ξ²), φ expanded in odd Chebyshev modes via PV∫w T_n/(t−a)dt = πU_{n−1}(a)), against the analytic K = P₀√(πL), V = π(1−ν)P₀L²/G and the weight-function K for a quadratic face load. Pasted output:
   ```
   equilibrium d(sigma_xx)/dx + d(tau_xy)/dy at (0.37L,0.53L) = -9.853e-11
   sigma_xx(face) = [-1.002415 -1.000308 -0.9998 -0.999692 -0.997585]  (expect -1)
   sigma_yy(face) = [-1. -1. -1. -1. -1.]  (expect -1)
   K from sigma_yy(L+1e-04,0): raw = 1.747520, bias-corrected = 1.77258678  (expect 1.77245385)
   u_y(x=0.0L) = 0.80000110  expect (1-nu)*sqrt(L^2-x^2)/G = 0.80000000
   DD solve [uniform P0=1, drained nu]:    K = 1.77245385, V = 2.51315776  (expect 2.51327412)
   DD solve [uniform P0=1, undrained nu_u]: K = 1.77245385, V = 1.88486832  (expect 1.88495559)
   DD solve [quadratic, drained nu]:        K = 0.88622693  (= 0.5*sqrt(pi L), weight function)
   PV-identity test n=1: dense = 3.139740001  pi*U_(n-1)(a) = 3.141592654
   ```
   Could have failed if: kernel sign wrong (V would come out negative), K–b relation wrong (K off), σ_xx on faces ≠ −1 (p_faces_0 wrong), or the Westergaard formulas wrong (equilibrium/face-BC residuals nonzero). The DD V is low by 4.6e−5 relative (clamped endpoint of the √-singular opening integral) — consistent truncation error, and K is exact to 8 digits.
2. **Local 1D poroelastic drainage layer (coefficient A).** Compared an implicit-FD solution of the layer problem (p(y,0)=2η, p(0,t)=0) against the analytic p = 2η·erf(y/2√(ct)) and against the surface-displacement formula 2(ν_u−ν)√(ct)/(G√π(1−ν)):
   ```
   oedometer identity: 1-eta*alpha(1-2nu)/(1-nu) = 0.750000000000 ; (1-nu_u)/(1-nu) = 0.750000000000
    t=1.0e-04: delta u numeric 2.84931401e-03  analytic 2.82094792e-03  ratio 1.010056
    t=1.0e-03: ... ratio 1.001002 ;  t=1.0e-02: ... ratio 1.000100
   profile check at t=1.00e-02: max|p_num - 2*eta*erf(...)| = 1.626e-05  (p scale 0.749)
   ```
   Ratios < 1 at t ≤ 1e−5 are grid under-resolution (dy = 0.004 > √(4ct)) — the check *should* fail there and does. Could have failed if the effective layer diffusivity were not exactly c, or ∫erfc ≠ √(4ct)/√π, or a factor-2 error in the two-faces × 2L bookkeeping. (The identity line verifies the algebraic cancellation that makes the diffusivity c.)
3. **Loading-2 long-time state.** (a) *Negative check:* the decaying harmonic field g = 1−Re[z/√(z²−L²)] (=1 on the slit, →0 at ∞) — which would give V2_inf = −πα(1−2ν)/2 — against the FD Laplace box solution:
   ```
   box X= 6.0:  p_FD(0.0,2.0) =  0.43488   g_harmonic =  0.10557
   box X=12.0:  p_FD(0.0,2.0) =  0.55396   g_harmonic =  0.10557
   min p_FD = 0.0000 (FD stays in [0,1]; g_harmonic is unbounded <0 at the tips) -> different
   ```
   The harmonic field is unbounded (−∞) at the tips and is *not* the parabolic limit (the box solution grows with box size). (b) Parabolic FD (slit held at p = 1):
   ```
   ct =  10.0:  p(0,2)=0.3957  p(0,5)=0.1102  p(3,3)=0.1572   (-> 1)
   ct =  50.0:  p(0,2)=0.5271  p(0,5)=0.2731  p(3,3)=0.3208
   ct = 200.0:  p(0,2)=0.6034  p(0,5)=0.3858  p(3,3)=0.4276
   ```
   p rises toward 1 at every fixed point (logarithmically slow; the finite box eventually caps it). This pins down the final state p = 1, σ = 0, uniform swelling ⇒ K2_inf = V2_inf = 0. Could have failed if a decaying steady state existed (p would stall below 1) or if my kernel/BC assembly were wrong (the earlier assembly bug produced exactly that: divergent column integrals).

## Questions (ambiguities and assumptions used)

1. **"t → ∞" for loading 2.** A pressure-maintained slit in an infinite 2D matrix has no steady state with p → 0 at ∞: p → 1 pointwise (2D recurrence). I used that limit ("vanish at infinity" holds at every finite t). The alternative "harmonic steady state" reading gives V2_inf = −πα(1−2ν)/2 = −0.8387 — rejected by check 3(a) (that field is unbounded at the tips and is not the diffusion limit). A third reading (finite reservoir) is not in the problem.
2. **K is the SIF of the total stress** σ_yy (not of Biot effective stress); the problem's σ_n = σ_yy convention suggests total stress.
3. **Loading 2 transiently gives negative opening/K** (faces interpenetrate). Linear theory allows it; no contact/constraint model is in the problem. Only the endpoint values are reported.
4. Plane-strain constants are used exactly as given (the two constitutive laws are self-consistent as a 2D theory; no 3D→plane-strain reduction was attempted).

## What is missing / uncertain

- The short-time coefficient A comes from a leading-order local (1D) analysis; verified against the FD layer solution (ratio → 1.0001), but not against a full 2D transient poroelastic simulation (tip-region and higher-order corrections O(ct/L²) neglected). Likewise K1(t), V1(t) at intermediate times were not computed (not requested).
- The claim V2(t) → 0 (return of the transiently-closed crack) follows from the pointwise limit state σ → 0, e = β_sw·δ (compatible); no full 2D poroelastic transient simulation was run to see the late-time approach directly.
- Numeric in check 1(b) has 4.6e−5 relative truncation error in V (endpoint clamp of the singular opening integral); K values are exact to the printed digits.
