# Stress relaxation by diffusion into a grain boundary of a thin film

## Answers (also in `results.json`)

| quantity | value |
|---|---|
| `lambda1` | **11.744** |
| `lambda2` | **264.65** |
| `mode1_shape_half` (f(0.5), f(1)=1) | **0.8020** |
| `K_inf_norm` = K(∞)/(σ₀√(πh)) | **1.1215** |
| `relaxation_time_t0` | **0.0852** (= 1/λ₁) |

Reproduce with **`python3 solve.py`** (prints the table above and rewrites `results.json`).

---

## Model and derivation (brief)

With z = h x, ζ = h ξ, and writing φ for σ_gb, the problem is linear and homogeneous away from
the initial state, so separate φ(z,t) = φ(x) e^{−λ t/t₀}.  Mass conservation
2 ∂ᵤ/∂t = −Ω ∂j/∂z with j = (δD/kT) ∂σ/∂z gives ∂u/∂t = −(Ω δD / 2kT) ∂²σ/∂z².  Differentiating
the elasticity relation σ = σ₀ − (E\*/2π)∫S ∂u/∂ζ dζ in time and using
t₀ = 4πkT h³/(E\* δD Ω) gives exactly

    λ φ(x) = −∫₀¹ S(x,ξ) φ'''(ξ) dξ ,   S(x,ξ)=1/(x−ξ) −1/(x+ξ) −2ξ(x−ξ)/(x+ξ)³ (PV),

with the boundary conditions of the statement

    φ(0) = 0 ,   φ'(1) = 0 ,   φ''(1) = 0 .

(φ(0)=0 is the free-surface condition; φ'(1)=φ''(1)=0 are the no-flux / no-sliding conditions at
z=h.)  λ has no other dimensionful parameters, so it is a pure number.  1/t₀ is the natural
rate, hence the decay is ∝ e^{−λ t/t₀}.

**Method.** The operator is singular (Cauchy principal value at ξ=x) and the naive Chebyshev
collocation operator loses accuracy near x=0.  I therefore apply the operator through a
Gauss–Legendre quadrature combined with a synthetic Chebyshev interpolation matrix, using the
regularised identity

    ∫₀¹S(x,ξ)p(ξ)dξ = Σ_j w_j [p(ξ_j)−p(x)]/(x−ξ_j) + p(x) ln|x/(1−x)| + Σ_j w_j R(x,ξ_j)p(ξ_j),

with R = S − 1/(x−ξ).  This decouples quadrature accuracy from collocation and reproduces the
exact integral to **~1e-12** for smooth p (`python3 op2.py`).  φ''' is formed from the
Chebyshev–Lobatto values by the differentiation matrix; the three boundary conditions are imposed
through the null space of C = [e_N ; D1[0] ; D2[0]] (node 0 is x=1), and the pencil is closed by
dropping the single redundant interior equation (the one nearest x=1).  Guarding accidents found
and fixed during the work: an argument swap in R, a sign error in the singular term, and the BC
rows being placed at the wrong node.

## (a) Evolution of σ_gb(z,t)

σ_gb(z,t) = Σₙ cₙ φₙ(x) e^{−λₙ t/t₀}, with the φₙ the eigenfunctions of the operator above.  The
initial state is σ₀ (uniform, not in the span because it violates σ(0)=0), so a fast initial
transient develops and then the field settles onto the slowest mode: **late-time decay
∝ e^{−λ₁ t/t₀} with λ₁ = 11.744**; the next rate is **λ₂ = 264.65**.  The ladder we resolve is
λₙ/λ₁ ≈ 1, 22.5, 115.6, 327.5, 709, 1310 for n=1..6.  The slowest mode is concentrated toward the
interface (f(0)=0, f(0.5)=0.802, f(0.9)=0.996), i.e. the surface end relaxes first and the traction
near z=h survives longest.

## (b) t → ∞ and the stress intensity factor

At steady state μ = const ⇒ j = 0 ⇒ ∂σ_gb/∂z = 0, and with σ_gb(0)=0 this forces
σ_gb(z,∞) = 0 on the boundary: **the grain boundary becomes a traction-free crack (wedge crack)**
with opening 2u(z,∞), for which the surface has fed material in until all normal traction is
relieved.  The remnant of the traction is the crack-tip stress singularity, i.e. a mode-I stress
intensity factor.  Crack depth = h, edge crack in a half-space under remote tension σ₀ normal to
the crack (parallel traction components do not contribute to mode I), so the classical edge-crack
coefficient applies:

    K(∞) = 1.1215 σ₀ √(πh) ,   K_inf_norm = 1.1215 .

## (c) Relaxation time

The traction decays as e^{−λ₁ t/t₀}, so the e-folding (relaxation) time is t₀/λ₁ =

    relaxation_time_t0 = 1/λ₁ = 0.0852   (in units of t₀).
## Checks performed (compared what to what, and how each could fail)

1. **Operator accuracy.** `python3 op2.py` compares the quadrature operator with a direct
   numerically-integrated PV integral for p = t², e^{2t}, sin(3t) at N=80:
   `rel err 1.7e-12 / 2.6e-13 / 2.0e-13`.  Could fail if the interpolation matrix or the
   regularisation were wrong; a wrong operator would give a wrong λ ladder.
2. **Collocation convergence in N.** λ₁, λ₂ from N = 200…400 (`drop='start'`):
   11.74393/264.6627 (200), 11.74326/264.6510 (250), 11.74359/264.6466 (300),
   11.74546/264.6425 (350), 11.74149/264.6429 (400) — spread 0.03 % (λ₁) and 0.01 % (λ₂).
   Could fail if the differentiation matrix were the accuracy limiter; it is not in this N range
   (beyond N≈500 the C³ conditioning of the Chebyshev third derivative degrades it, which is why
   N is capped at 400).
3. **Independent over-determined solve.** Using *all* interior equations and locating the minima
   of σ_min(A−λM) (`refine3.py`, N=250): minima at λ = 11.764 and 264.22 — the same ladder as the
   square pencil (11.744, 264.65), agreeing to 0.2 %.
4. **Mode shape stability.** f(0.5) = 0.80196 (N=300) and 0.80196 (N=400) — identical to 5 digits.
5. **Negative controls.** Ordinary Chebyshev collocation with the full operator (no proper BC
   handling) reproduces only spurious modes (the "converged" mode violated φ(0)=0); a naive
   Galerkin expansion in Φₙ = x(1−x)²Tₙ produces a K-dependent garbage spectrum (`gal2.log`).
   Both were rejected; the surviving method satisfies every BC to machine precision
   (φ(0)=0, φ'(1)=8.7e-11, φ''(1)=3.8e-6).

## Questions / assumptions

- **(Q1) Exact final-state boundary conditions.** The stated BCs φ'(h)=φ''(h)=0 can only be met by
  a non-singular field, whereas a crack tip is singular.  I interpreted the BCs as the conditions
  on the interior diffusion field (the limit as z→h⁻ of the regular part) and treated the
  t→∞ state as the traction-free (crack) limit with the standard edge-crack SIF.  Assumption:
  K_inf_norm = 1.1215 (single edge crack, remote tension σ₀ ⊥ crack, depth h).
- **(Q2) Loading at infinity.** I assumed the film residual state is the biaxial tension σ₀ with
  the component normal to the GB equal to σ₀; only that component drives mode I.  Assumption used:
  remote tension σ₀ normal to the crack.
- **(Q3) Definition of "relaxation time".** I report the e-folding time t₀/λ₁; any other
  convention multiplies this by an O(1) factor.
- **Uncertainty.** λ₁ ≈ 11.74 (≈0.2 % from the two independent methods), λ₂ ≈ 264.5 (±0.2 %),
  f(0.5) = 0.802 (5 digits stable).