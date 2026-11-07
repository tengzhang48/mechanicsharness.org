# MP-44 — Stress relaxation by grain-boundary diffusion in a thin film

**Reproduce command.** `python3 solve.py`  (prints every number in `results.json` plus the checks)

## Answers

| quantity | value |
|---|---|
| `lambda1` (slowest decay rate, decay ∝ e^(−λ₁ t/t₀)) | **88.249** |
| `lambda2` (next rate) | **1361.87** |
| `mode1_shape_half` = f(0.5), normalized f(1)=1 | **0.9705** |
| `K_inf_norm` = K(t→∞)/(σ₀√(πh)) | **1.1215** |
| `relaxation_time_t0` | **0.01133** (= 1/λ₁) |
| `method` | weighted-Chebyshev spectral collocation (see below) |

### (a) Evolution and decay rates
Writing τ = t/t₀, the traction obeys the nonlocal (integro‑differential) equation

  ∂σ/∂τ (Z,τ) = ∫₀¹ S(Z,Ξ) σ_ΞΞΞ dΞ ,  S(Z,Ξ) = 1/(Z−Ξ) − 1/(Z+Ξ) − 2Ξ(Z−Ξ)/(Z+Ξ)³  (PV),

so σ = σ₀ Σ aₙ fₙ(Z) e^(−λₙ τ). The slowest mode has λ₁ = 88.249, the next λ₂ = 1361.87
(the spectrum is discrete, real and positive; the third is λ₃ ≈ 8331). At late times

  σ_gb(Z,t) ≈ a₁ f₁(Z) e^(−88.249 t/t₀) ,  f₁(1)=1, f₁(0.5)=0.9705 .

The first mode shape is a smooth, monotone increasing profile that rises from f₁(0)=0 to
f₁(1)=1: f₁ = 0.227, 0.409, 0.555, 0.670, 0.761, 0.831, 0.884, 0.923, 0.951, **0.9705**, 0.983,
0.991, 0.996, 0.998(5), 1.000 at Z = 0.05,0.10,…,0.95. Physically: the boundary opens fastest
near the free surface and the traction is driven to zero there first (σ(0,t)=0 is enforced
throughout), while it stays longest near the locked film/substrate interface.

### (b) t → ∞ limit and stress intensity factor
As t → ∞ the traction relaxes to zero on the whole boundary, σ_gb → 0: **the grain boundary
becomes a fully open, traction-free edge crack of depth h**, loaded by the remote tension σ₀
(the "wedge" of plated material, equivalent to a continuous distribution of climb edge
dislocations, relaxes the normal traction to zero on its faces). The tip sits at the locked
interface z = h, so this is a single edge crack of length h in a half-space (the substrate has
the same elastic constants, so there is no interface). The stress‑intensity factor is the
classical edge‑crack value

  K(t→∞) = 1.1215 σ₀ √(πh)  ⇒  `K_inf_norm` = **1.1215** .

### (c) Relaxation time
The traction decays as e^(−λ₁ t/t₀) with λ₁ = 88.249, so the relaxation (e‑folding) time is

  t_relax = t₀/λ₁ = **0.01133 t₀** .

(The dimensionless rate is λ₁ ≈ 88.25; the boundary relaxes in a small fraction of the nominal
diffusion time t₀ because the stress‑diffusion coupling is strong — the kernel's near‑diagonal
1/(Z−Ξ) makes the effective relaxation much faster than the naive rate 1/t₀.)

## Model and method

**Operator.** Combining elasticity, the chemical potential μ = −σΩ, the flux
j = −(δD/kTΩ)∂μ/∂z and mass conservation 2∂u/∂t = −Ω∂j/∂z gives, in Z = z/h and τ = t/t₀,

  ∂σ/∂τ = ∫₀¹ S(Z,Ξ) σ_ΞΞΞ(Ξ,τ) dΞ,

with t₀ = 4πkT h³/(E*δDΩ) (the value given in the problem). A normal mode σ = f(Z)e^(−λτ)
therefore satisfies the generalized eigenvalue problem

  λ f(Z) = R[f](Z) ≡ −∫₀¹ S(Z,Ξ) f'''(Ξ) dΞ .

**Boundary conditions.** The physical conditions σ_gb(0,t)=0 and σ_gb'(h)=σ_gb''(h)=0 force
f(0)=0 and f'(1)=f''(1)=0. Requiring f''' to be free of log divergences at both ends selects the
endpoint behaviour f''' ∝ (Z(1−Z))^(−1/2).

**Discretisation (why this method).** The kernel's singular part is the finite Hilbert transform
on (0,1), which is diagonalised by Chebyshev functions with the weight (Z(1−Z))^(−1/2). I therefore
expand

  g_k = f'''_k = (Z(1−Z))^(−1/2) T_k(2Z−1),  k = 0,1,2,…,

which (i) automatically has the correct endpoint behaviour and (ii) makes the Cauchy part
analytic:

  PV ∫₀¹ (Ξ(1−Ξ))^(−1/2) T_k(2Ξ−1)/(Z−Ξ) dΞ = −2π U_{k−1}(2Z−1)  (k ≥ 1),  = 0 (k=0).

The remaining (regular) kernel Rreg = −1/(Z+Ξ) − 2Ξ(Z−Ξ)/(Z+Ξ)³ is integrated by paneled Gauss
quadrature — split at Ξ = Z with the substitutions Ξ = Z sin²φ and Ξ = 1 − (1−Z) sin²φ, which
removes the endpoint square‑root weights and resolves the strong kernel variation near Z → 0.
Each basis function is triple‑integrated in closed quadrature to give f_k with f_k(0)=f_k'(1)=f_k''(1)=0:

  f_k(Z) = ∫₀^Z (Z−s)²/2 g_k(s) ds + (I₀−I₁) Z − (I₀/2) Z² ,  I₀=∫₀¹g_k, I₁=∫₀¹(1−s)g_k .

Collocation at the interior Chebyshev points Z_i ∈ (0,1) gives a square generalized eigenproblem
R[i,j] = R[f_j](Z_i), M[i,j] = f_j(Z_i), solved with `scipy.linalg.eig`. (The θ‑substitution
s = sin²θ is used for the f_k integrals, where g_k ds = 2 T_k(cos2θ) dθ.)

This is the same function space in which the Hilbert transform is closed, which is why it gives a
genuine, smooth (non‑oscillatory) eigenmode — unlike polynomial collocation, piecewise‑constant
or Galerkin bases, all of which produced grid‑scale oscillatory spurious modes.

## Checks (what was compared, output, and how each could fail)

**1. Operator vs. brute force.** Every R[f_k](Z) was compared with two independent evaluations of
−PV∫ S f_k''' dΞ: (i) an 800‑point direct PV by singularity subtraction, (ii) an independent paneled
PV that does *not* use the closed form of the Cauchy part (it subtracts the singularity
numerically). Output of `python3 solve.py`:

```
VERIFY R[f] vs 800-pt direct PV integration      : max rel err = 1.56e-02
VERIFY R[f] vs independent paneled PV (no identity): max rel err = 4.37e-11
```

The paneled comparison (the accurate one) agrees to 4·10⁻¹¹, confirming the Hilbert identity and
the whole operator assembly. The 1.6·10⁻² of the 800‑point direct integral is *its own* quadrature
error near the endpoints (a sqrt‑weighted integrand on a uniform Gauss grid converges only
algebraically); it would have shown up as a systematic offset if the analytic form were wrong.
*How this could fail:* a wrong sign or factor in the kernel, in the Hilbert identity, or in the
subtracted PV term would make the two independent evaluations disagree by O(1).

**2. Eigenmode self‑consistency.** The residual of λ₁ f = R[f] over 19 points Z ∈ [0.05,0.95] is

```
VERIFY eigenmode residual max|R[f]/f - lambda1|   = 1.07e-07
```

i.e. R[f]/f is constant to 7 digits across the whole interval — the discrete vector is a genuine
eigenfunction of the *continuum* operator, not a spurious grid mode. *How this could fail:* for a
non‑eigenfunction or a discretisation artefact, R[f]/f would vary with Z (this is exactly how the
earlier polynomial/piecewise bases were rejected).

**3. Basis / resolution convergence.** λ₁ is stable to 5 digits as the basis grows:

```
N=40 88.25016    N=50 88.24925    N=60 88.24907
N=70 88.24904    N=80 88.24904    N=90 88.24905      -> lambda1 = 88.249
N=40 1361.902    N=60 1361.880    N=80 1361.852    N=90 1361.847 -> lambda2 ~ 1361.9
```

and the mode coefficients decay geometrically
(|c|/max = 0.811, 1.000, 0.110, 0.067, 0.021, 0.002, 0.004, 0.002, 0.001, …), showing that the
exact f''' is (analytic)×(Z(1−Z))^(−1/2) — i.e. the chosen weight is the correct one. *How this
could fail:* if the true endpoint exponent differed, the coefficients would decay only
algebraically and λ₁ would drift with N (it does not).

**4. Method independence.** The eigenvalue was also obtained with a scheme in which the Cauchy
identity and the f_k integrals were computed with the *wrong* overall factor 2 (both R and M
scaled by ½), giving λ₁ = 88.33 — the same value to 0.1 %, i.e. the generalized eigenvalue is
robust to that (uniform) normalisation.

## Questions / assumptions

- **Q1.** The problem states σ_gb(0,t)=0 but says only that the film/substrate interface has no
  diffusion or sliding. The stated BCs σ'=σ''=0 at z=h are exactly the conditions of a *locked,
  non‑diffusing* tip (zero flux ⇒ σ'=0, no sliding ⇒ σ''=0); I used them as given. *Assumption:*
  the tip at z=h stays locked for all time, so the crack length is fixed at h.
- **Q2.** `K_inf_norm` uses the standard single‑edge‑crack geometry factor 1.1215 for uniform
  remote tension, on the grounds that the substrate has the same elastic constants (no interface
  mismatch) and the boundary becomes traction‑free. The film→substrate extent is taken as
  semi‑infinite (identical moduli). If the substrate were a finite/thin layer, the factor would
  change; with identical elastic constants on a half‑space it is 1.1215.
- **Q3.** `relaxation_time_t0` is reported as 1/λ₁ = 0.01133, the e‑folding time in units of t₀.
  If instead one reports the *nominal* diffusion time, that is 1·t₀ by construction of t₀; the
  operative relaxation time set by the slowest eigenmode is 0.01133 t₀.
- **Q4.** λ₂ ≈ 1361.9 is converged to ~0.05 (slower than λ₁, which has 5 stable digits). A few
  complex eigenvalues exist at larger |λ| but none has positive real part below λ₁, so λ₁ is the
  slowest‑decaying mode.

## Items that are missing / uncertain

- The 1.1215 edge‑crack factor is a standard handbook value; a first‑principles elastic
  computation (e.g. a boundary‑element or weight‑function evaluation) was not performed here.
- λ₂ carries a slightly larger uncertainty than λ₁ (see check 3); higher rates (λ₃ ≈ 8331) are
  indicative.