# ANSWER — stress relaxation by diffusion into a grain boundary of a thin film

## Model (as given, reduced to dimensionless form)

With `x = z/h`, `s = t/t₀`, `t₀ = 4πkT h³/(E* δ_gb D_gb Ω)`, and `σ` normalized by `σ₀`:

* **Kinetics.** μ = −σΩ, `j = (δ_gb D_gb /kT) ∂σ/∂z`, mass conservation `2∂u/∂t = −Ω ∂j/∂z` give
  `∂u/∂t = −(Ω δ_gb D_gb /2kT) ∂²σ/∂z²`.
* **Elasticity** (given): `σ = 1 − K[∂ũ/∂x]`, `K[g](x) = PV∫₀¹ S(x,ξ) g(ξ) dξ`,
  `S = 1/(x−ξ) − 1/(x+ξ) − 2ξ(x−ξ)/(x+ξ)³`.
* Eliminating `u`: **`∂σ/∂s = K[∂³σ/∂x³]`**, with `σ(0)=0`, `σ_x(1)=σ_xx(1)=0`, `σ(x,0)=1`.

Why this operator: on the whole line the leading (Hilbert) part of `K` gives the symbol
`−π|k|³`, i.e. strong dissipation — relaxation, as expected physically.

## Answers

### (a) Evolution of σ_gb and decay rates

* Short times: the IC `σ≡σ₀` is incompatible with `σ(0)=0`, so a relaxed boundary layer
  forms at the free surface and propagates into the film; since the relaxation rate of a
  local (Fourier) mode of wavelength ℓ is ~π/ℓ³, the front advances like `z/h ~ s^{1/3}`.
  Until the front arrives, the traction deeper down stays ≈ σ₀ (the initial rate is zero
  because σ₀''' = 0).
* Late times: the solution is a sum of discrete exponentials; the traction decays as
  `e^{−λ₁ t/t₀}` with the **mode-1 shape** `f₁(x)`, f₁(0)=0, f₁(1)=1 by normalization,
  monotonically increasing from the surface to the interface end (values: f₁ ≈ 0.11 at
  x=0.05, 0.48 at 0.25, 0.80 at 0.5, 0.96 at 0.75, 1.00 at 1).

**λ₁ = 11.742668, λ₂ = 264.6342** (λ₃ = 1357.45, λ₄ = 3845.2, λ₅ = 15382.5;
λₙ growing roughly like n⁴, consistent with a 4th-order dissipative operator).
The slowest mode dominates after s ≳ 0.2 (e.g. at s = 0.5 the second mode is down
by e^{−126}). Completeness check: the residues of σ̂(0.5;p) at the poles
(1.2266, 0.0622, −0.3586, −0.0667, 0.0380, …) sum to 0.90 with 5 poles and are
still increasing toward σ(0.5,0⁺)=1 as poles are added, as they must for a
complete discrete spectrum; the late-time traction at any x>0 is
σ(x,s) ≈ 1.227·f₁(x)·e^{−11.7427 s} (at x=0.5).

The eigenfunctions are "crack-like" at the film/substrate end: their dislocation density
`ρ = ∂³f/∂x³ ~ (1−x)^{−1/2}`. This follows from the tip balance: if `ρ ~ (1−x)^β`, the
PV term of `K[ρ]` behaves like `−π cot(πβ)(1−x)^β`, and matching orders in
`−λf = K[ρ]` requires `cot(πβ)=0`, i.e. β = −1/2, 1/2, 3/2, …; three integrations give
`f − f(1) ~ (1−x)^{β+3}`, so **all** these satisfy `f'(1)=f''(1)=0` (for β=−1/2:
f−f(1) ~ (1−x)^{5/2}); the modes selected by the dynamics (and by the IC) are the
crack-like branch β = −1/2.

### (b) t → ∞: the GB becomes a traction-free crack ("diffusion wedge")

The plated wedge `2u_∞(z)` has the opening of a surface crack of depth h in the
half-space under the initial (far-field) tension σ₀: traction on the GB → 0, stress
concentrated ahead of the wedge tip at the film/substrate interface.

**K(t→∞) = 1.12152 σ₀√(πh)** — the classical edge-crack value (1.1215), obtained from
the steady singular-integral solve `K[ρ_∞] = 1` (residual < 1e−15) via
`K/(σ₀√πh) = −π√2 · [ρ_∞√(1−x)]_{x→1}`, and independently by extrapolating
`σ(r)√(2r)` ahead of the tip (1.1264, 1.1230, 1.1220 at r = 10⁻², 3·10⁻³, 10⁻³ →
1.12152). The agreement with the known handbook constant is also a validation of the
given kernel S.

### (c) Relaxation time

E-folding time of the slowest mode: **t₀/λ₁ = 0.0852 t₀** (≈ 0.39 t₀ for a 100-fold
reduction, ln100/λ₁).

## Method (summary; details in solve.py, solver5/7/10.py, crack.py)

1. **Exact operator machinery.** For basis `ρ_n = (1−x)^{−1/2}T_n(2x−1)`,
   `K[ρ_n](x) = q_n(x)·J₀(x) + P_n(x) + Rat_n(x)` with
   `J₀ = PV∫(1−ξ)^{−1/2}/(x−ξ)dξ = δ^{−1}ln((1−δ)/(1+δ))`, δ=√(1−x) (closed form,
   verified to 1e−14 against subtraction quadrature), `P_n` the PV-subtraction
   polynomial (built by exact quadrature + interpolation), `Rat_n` the image-kernel
   integrals by Gauss–Jacobi quadrature (exact for the weight). Every piece verified
   against brute quadrature and finite differences.
2. **Steady crack**: collocate `K[ρ_∞] = 1` → ρ_∞, K∞ (converged at N=24 to 8 digits).
3. **Decay rates**: Laplace transform in s. `σ̂(x;p) = f[ρ̂](x)` with
   `f[ρ] = x∫ξρ − (x²/2)∫ρ + ½∫₀ˣ(x−ξ)²ρ dξ` — this construction satisfies
   `f'''=ρ`, `f(0)=f'(1)=f''(1)=0` identically, so **all BCs are built in**. The
   resolvent equation `pσ̂ − 1 = K[ρ̂]` becomes a *square* linear system
   `(pF − H)a = 1` at Chebyshev collocation points. Decay rates are poles: minima of the
   smallest singular value of `(pF − H)`; refined by 1-D search to ~1e−9. The mode
   shape is the null vector (residue) at the pole. The p=0 solve reproduces the steady
   crack coefficients exactly (max|a(0)+a_crack| = 0.0).

## Checks performed (what was compared, output, how it could have failed)

1. **Gauss–Jacobi quadrature vs exact Beta moments** (scratch1):
   `k=0..5: rel ≈ 1e−15`. Could have failed if the Golub–Welsch Jacobi recurrence
   coefficients were wrong.
2. **Closed forms J₀^{±1/2} vs brute PV quadrature** (scratch1/2): agreement
   `~1e−11` (sqrt) and `~1e−14` (invsqrt). Failure mode: wrong PV-subtraction or sign.
3. **Full kernel identity** `∫₀¹S(x,ξ)dξ = ln((1+x)/(1−x)) − 4 + 6x/(1+x) − 2x²/(1+x)²`
   (hand-derived) vs brute: `diff ≤ 1.6e−11` at x = 0.01…0.99. Validates PV + both
   image terms together. Failure mode: wrong image kernel (would break K∞).
4. **Exact basis machinery vs brute** (`h_pieces` vs direct quadrature): max err
   `5e−14` (N up to 32); derivative pieces vs FD: ~FD-noise level. Failure mode: the
   Leibniz/Taylor-moment terms (two sign errors found and fixed during development).
5. **Steady crack**: K∞ = 1.1215223 by two independent routes (tip coefficient; stress
   ahead of tip) and vs the literature edge-crack value 1.1215. Failure mode: any error
   in S would shift this away from 1.1215.
6. **λ convergence** (N=16,20,24,28,32,36,40,44): λ₁ = 11.742665 → 11.7426681
   (stable to ~2e−7); λ₂ = 264.6345 → 264.63421 (stable to ~1e−5). Mode-1 shape
   f(0.5)/f(1) = 0.80195650 identical to 8 digits at N=32, 40, 44 (also
   f(0.25)=0.483408, f(0.75)=0.963257). Failure mode: collocation pollution /
   non-representable eigenfunctions (this DID break several earlier, discarded
   formulations — see PROCESS.md).
7. **Mode-1 eigen-relation, independent brute quadrature** at x = 0.3/0.5/0.7:
   `|K[ρ₁] + λ₁ f₁|/|K[ρ₁]| = 1.2e−4, 1.3e−4, 1.7e−4` (residual is limited by the
   finite pole-refinement accuracy, smin ≈ 1e−6). Failure mode: a spurious pole.
8. **p = 0 resolvent = steady crack**: identical coefficients (diff 0.0).
9. **Independent coarse evolution** (SVD-regularized Galerkin, different basis):
   slowest rate −11.1…−11.5 bracketing λ₁ = 11.74 (that method is only ~3-digit
   accurate; used as consistency evidence only).

## Results (results.json)

```
lambda1            = 11.74266781911283
lambda2            = 264.6342031625176
mode1_shape_half   = 0.8019565042704974   (f(0.5) with f(1)=1)
K_inf_norm         = 1.121522256600971
relaxation_time_t0 = 0.08515952383259624
method             = "resolvent poles of (p-L)^-1 via smallest-singular-value minima ..."
reproduce_command  = "python3 solve.py"
```

Reproduce: `python3 solve.py` (≈ 9 s; prints all values above and writes results.json).

## Questions / assumptions

1. **Initial stress state**: assumed to be a uniform in-plane tension σ₀ perpendicular
   to the GB throughout film+substrate (consistent with "the two form an elastic
   half-space" and σ_gb(z,0)=σ₀; it is also the state whose relaxation gives the
   classic surface-crack K = 1.1215σ₀√(πh) the problem asks for in (b), which
   supports the reading).
2. **λ₂ is the next discrete rate after λ₁** — found as 264.63; note the large gap
   λ₂/λ₁ ≈ 22.5 (the spectrum is strongly non-uniform, λₙ growing ≈ n⁴, consistent
   with a 4th-order dissipative operator).
3. **"Relaxation time"** reported as the e-folding time 1/λ₁ = 0.0852 t₀; a
   "complete" relaxation (2 orders of magnitude) takes ≈ 0.39 t₀.
4. The problem statement's diffusion/elasticity equations were used exactly as given;
   the sign conventions (μ = −σΩ etc.) imply the dissipative direction used here
   (checked: the hump-relaxation test gives σ_t < 0 where σ > 0 over most of the bump).

## What is missing / uncertain

* λ₂ (and higher rates) are determined to ~6 significant digits (λ₂ = 264.6342 ±
  ~1e−5 from N-trend); I have no analytic or literature cross-check for them; only
  λ₁ is corroborated (crudely) by an independent evolution. The mode-1 shape at
  0.5 is reproduced to 8 digits across N = 32–44; its absolute accuracy is judged
  ~1e−4 (level of the brute eigen-relation residual).
* K∞ matches the handbook value 1.1215 (which I quote from knowledge of the standard
  edge-crack result, not from a web source); my own two extractions agree to 5 digits.
* I did not find the source publication of this problem; the numbers are my own
  computation. Several mathematically defective formulations (collocation pencils on
  smooth bases, wrong tip-exponent classes) were tried and discarded — documented in
  PROCESS.md — before the crack-like basis + resolvent-pole method converged.
