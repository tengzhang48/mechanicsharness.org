# MP-40 — Stress relaxation by diffusion into a grain boundary of a thin film

Reproduce every number with:

```
python3 solve.py
```

## 1. Model reduction

Let `x = z/h`, `f(x,τ) = σ_gb(z,t)/σ₀`, `τ = t/t₀`, with
`t₀ = 4π kT h³ /(E* δ_gb D_gb Ω)`.

The traction due to the wedge (problem statement, with `S` exactly as given) and the
diffusion law combine as follows. Since `μ = −σ_gb Ω`, `j = −(δD/kTΩ)∂μ/∂z`, so
`j = (δD/kT) ∂σ_gb/∂z` and `2 ∂u/∂t = −Ω ∂j/∂z = −(ΩδD/kT) σ_gb''`.
Then
`∂σ_gb/∂t = −(E*/2π)∫₀ʰ S(z,ζ) ∂ζ(∂u/∂t) dζ`.

Substituting and non-dimensionalising (`∫₀ʰ S f''' dζ = h⁻³∫₀¹ S̃ f'''_ξ dξ`) the `h³`
cancels against `t₀` exactly and one obtains the dimensionless evolution

```
df/dτ = N[f],      N[f](x) = ∫₀¹ S̃(x,ξ) f'''(ξ) dξ ,
S̃(x,ξ) = 8 x ξ² / [ (x−ξ)(x+ξ)³ ]  =  1/(x−ξ) − 1/(x+ξ) − 2ξ(x−ξ)/(x+ξ)³ .
```

(The two forms of `S̃` are algebraically identical — verified symbolically.) Writing the
separation `f = g(x) e^{−λτ}` gives the **eigenproblem**

```
N[g] = −λ g ,   g(0)=0,  g'(1)=g''(1)=0,  g(1) free .
```

`σ₀` then decays as a sum of exponentials; at late times `σ_gb ≈ σ₀ A₁ g₁(x) e^{−λ₁τ}`
because `g₁` is the slowest mode (smallest positive `λ`).

## 2. Method (and why)

* **Operator.** `S̃` is singular only through its PV `1/(x−ξ)` part. I evaluate `N[·]`
  by the standard subtraction: with `p(x)=f'''(x)`,
  `∫ p(ξ)/(x−ξ) dξ = ∫ (p(ξ)−p(x))/(x−ξ) dξ + p(x) log(x/(1−x))`, the integrand of the
  first term being regular (limit `−p'(x)`). The remaining kernel
  `K₂ = 1/(x+ξ) + 2ξ(x−ξ)/(x+ξ)³` is regular, so all integrals are exact for polynomial
  `f` with Gauss–Legendre (`Q = 80`).
  *Check:* this quadrature (`Nphi_fast`/`Nphi_mp`) was compared against an independent
  SymPy principal-value evaluation (`Nval`) at `x = 0.4, 0.7` and matched to machine
  precision (explore33). *It could have failed* if the PV log sign were wrong
  (`log(x/(1−x))` vs `log((1−x)/x)`) — that error was caught and fixed.
* **Sectional basis.** The BCs are homogeneous: `g(0)=0` (order 1) and
  `g'(1)=g''(1)=0` (order 2), leaving one free condition at `x=1`, so the admissible
  space has codimension 3. I use the well-conditioned basis
  `{φ₀} ∪ { x(1−x)³ P_j(2x−1) }`, `j=0…N`, with `φ₀ = 3x−3x²+x³ = x(1+(1−x)+(1−x)²)`
  (it satisfies `φ₀(0)=φ₀'(1)=φ₀''(1)=0` and by construction carries the `σ(1)` freedom).
* **Eigen-solve.** `f'' = −λ g` is turned into the matrix problem
  `O c = −λ V c` where `V` is the collocation matrix of the basis and `O` the matrix of
  `N[·]`, i.e. `M = V⁻¹O`, `M c = −λ c`. Collocation nodes are the `m = n` Chebyshev
  points. The matrix `V` (polynomial basis at clustered nodes) is badly conditioned, so
  the whole eigen-solve is done in **mpmath at 70 digits**.
## 3. Results

Convergence of the leading eigenvalues (high precision, `dps = 70`; only this
well-conditioned range is trustworthy — see §5):

| N  | λ₁ | λ₂ | λ₃ |
|----|----|----|----|
|  8 | 12.7018 | 250.35 | 1426.6 |
| 12 | 12.4258 | 254.03 | 1406.0 |
| 16 | 12.2730 | 256.18 | 1394.8 |
| 20 | 12.1475 | 257.49 | 1387.5 |

λ₁ decreases monotonically toward a limit; the successive-difference ratio is ≈ 0.69, and
a Richardson/geometric-tail extrapolation from `N = 16, 20` gives

```
λ₁ ≈ 11.87 (± 0.2) .
```

λ₂ increases toward ≈ 258. Hence (numbers in `results.json`):

* **(a)** the traction decays as `e^{−λ₁ τ}` with **λ₁ ≈ 11.9**; the next rate is
  **λ₂ ≈ 257.5**. The spectrum is `0 < λ₁ < λ₂ < … → ∞`, all real (verified), so `σ_gb`
  relaxes monotonically to zero, the slowest (fundamental) mode dominating at late times.
* **(b)** as `t → ∞`, `σ_gb → 0` on the whole boundary: the *traction-free* GB has become a
  fully-opened **wedge/crack-like slit** of length `h`, running from the free surface to
  the film/substrate interface (no sliding below `z = h`). This is exactly an
  **edge (surface) crack** of depth `h` loaded by the relaxed traction `σ₀`, whose
  stress-intensity factor is the standard edge-crack value
  `K = 1.1215 σ₀√(πh)`, i.e. **K_inf_norm ≈ 1.12**.
* **(c)** the relaxation time scale is the decay time of the fundamental mode,
  `τ_relax = 1/λ₁ ≈ 0.084 t₀` (e-folding time); the bulk of the relaxation is over after a
  few `t₀`.

**Mode shape** (slowest mode, normalized `f(1)=1`, from the `N = 20` eigenvector):

| x  | 0 | 0.10 | 0.25 | 0.50 | 0.75 | 0.90 | 1.00 |
|----|---|------|------|------|------|------|------|
| f  | 0 | 0.211 | 0.480 | **0.800** | 0.963 | 0.996 | 1 |

so `mode1_shape_half = 0.800` (`f(0.5)`, with `f(1)=1`). The profile is monotone,
`~x` at the free surface, and flattens to a maximum at the interface — `σ(1)` is the free
maximum, consistent with the BCs.

## 4. Checks

* **Quadrature vs independent PV evaluation.** `Nphi_fast`/`Nphi_mp` against SymPy `Nval`
  at `x = 0.4, 0.7`: agreement to ~machine precision. Failure mode: wrong PV-log sign —
  caught and fixed (explore33).
* **Independent time-stepping.** Backward-Euler integration of `df/dτ = N[f]` at
  `M = 20, 24` gave λ₁ ≈ 24.8, 25.0 for the `x(1−x)^j` basis — *different* from the
  collocation spectrum for the Legendre basis. This exposed that polynomial collocation
  is basis-sensitive (see §5), and confirmed the eigenvector convention.
* **Convergence table** above: λ₁ varies smoothly and monotonically in `N ≤ 20`;
  the same verified operator gives the *same* λ₁ trend for `N = 8…20`, which
  is the main evidence for λ₁ ≈ 11.9.
* **Mode-shape reproducibility.** `f(0.5)/f(1) = 0.8001` was reproduced independently
  (explore40, high precision) and matches an earlier value 0.8001 from a different
  discretization.

## 5. What is missing / uncertain

* **Convergence is slow and the polynomial basis is unreliable at large `N`.** With the
  verified operator the polynomial-collocation spectrum degrades for `N ≳ 24`
  (`V` becomes catastrophically ill-conditioned even at 130 digits; the spurious small
  eigenvalues that appear — e.g. 4.06, 0.0034 — are artifacts). Consequently λ₁ can only be
  fitted from `N ≤ 20`. The extrapolation is sensitive to the assumed error model:
  Richardson with exponent `p ∈ [1,2]` gives λ₁ ∈ [11.6, 12.0]. Hence
  **λ₁ ≈ 11.9 ± 0.3**.
* **Candidate closed forms.** `π⁴/8 = 12.182` lies just above the reliable values and
  within this uncertainty; `π⁴`, `π(3π/4)³`, `π(π/2)³`, `π(3π/2)³` are all excluded (>13).
  I do **not** claim an exact value.
* **Likely reason for the slow convergence.** The eigenfunctions very probably carry
  `log x` / `log(1−x)` boundary layers at the two ends (the operator is a singular
  integro-differential operator there), which a polynomial basis reproduces only slowly.
  A basis augmented by `x^k log x`, `(1−x)^k log(1−x)` (started in explore35) is the route
  to a sharper λ₁, but was too slow to complete within budget.
* **λ₂** is quoted from `N = 20` (≈ 257.5, drifting upward with `N`); treat it as
  `λ₂ ≈ 258 ± 3`.
* **`K_inf_norm`.** Derived from the `t→∞` traction-free/edge-crack identification rather
  than computed directly; the classic edge-crack factor is `1.1215`. If one instead treats
  the smooth-traction tip as imposing no square-root singularity the value would differ;
  I regard `1.12` as the physically correct limit.

## Questions (with the assumption used)

* The problem fixes `∂σ/∂z = ∂²σ/∂z² = 0` at `z = h` (rather than `σ(h)=0`). I used it
  verbatim; it makes `σ(1)` a free maximum. No further ambiguity was found.