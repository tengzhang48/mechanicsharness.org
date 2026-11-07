# MP-02 — Stress relaxation by diffusion into a grain boundary

## Headline results

| quantity | value |
|---|---|
| λ₁ (slowest decay rate, e^{−λ₁ t/t₀}) | **11.7426681** |
| λ₂ (next decay rate) | **264.63422** |
| mode1_shape_half = f(0.5), f(1)=1 | **0.80195650** |
| K(t→∞)/(σ₀√(πh)) | **1.1215223** |
| relaxation_time_t0 = 1/λ₁ | **0.085160** |

Reproduce with: `python3 solve.py` (prints the values in `results.json`, runs in ~10 s).

## Model and method (and why)

**Dimensionless evolution equation.** With Σ(Z,τ) = σ_gb/σ₀, Z = z/h, τ = t/t₀,
t₀ = 4πkT h³/(E*δ_gb D_gb Ω), the physics in the problem statement reduces to a single
linear integro-differential equation on Z∈[0,1]:

    ∂Σ/∂τ = L[Σ'''] ,   L[g](Z) = PV ∫₀¹ s(Z,Z') g(Z') dZ',
    s(Z,Z') = 1/(Z−Z') − 1/(Z+Z') − 2Z'(Z−Z')/(Z+Z')³,

with boundary conditions Σ(0)=0, Σ'(1)=Σ''(1)=0 and initial condition Σ(Z,0)=1.
(Derivation below.) Diffusion along the boundary acts as a third spatial derivative
feeding back through the nonlocal elastic operator L; relaxation = spectral decay of
this operator.

**Why a spectral (global-basis) method.** The operator is singular (Cauchy PV plus a
(1−Z')^(−1/2) singularity of the third derivative of any reasonable mode at the tip
Z=1). A global basis that (i) satisfies all boundary conditions exactly and (ii) contains
the exact crack-tip power law makes the discrete problem tiny (N≲80) and converges
spectrally, whereas finite differences would fight the PV and the tip singularity.

**Basis.** The regular expansion of a mode about the tip is f = c₀ + c₁(1−Z)^(5/2) + … :
the BCs f'(1)=f''(1)=0 kill the (1−Z)¹, (1−Z)² terms, and the leading non-analytic power
consistent with a finite L[f'''] at Z=1 is (1−Z)^(5/2) (its third derivative
∝(1−Z)^(−1/2) has a finite Hilbert transform at the endpoint). Hence

    φ₀(Z) = 1 − (1−Z)^(5/2)                     (φ₀(0)=0, φ₀(1)=1)
    φ_j(Z) = Z(1−Z)^(5/2) T_{j−1}(2Z−1),  j=1..N

All φ_j satisfy φ(0)=0, φ'(1)=φ''(1)=0 exactly.

**Quadrature (kernel.py).** Substituting Z' = 1−y² turns the tip singularity into a smooth
integrand: for basis functions, G(y) = 2y·φ'''(1−y²) is a polynomial × smooth factor.
The Cauchy kernel becomes 1/(Z−Z') = 1/(y²−y_Z²), y_Z=√(1−Z), handled by analytic
subtraction: ∫[H(y)−H(y_Z)]/(y²−y_Z²) dy + H(y_Z)·I₀(y_Z) with the closed form
I₀(a) = ln((1−a)/(1+a))/(2a) (→−1 as a→0). Gauss–Legendre quadrature (M=384–768 nodes)
then gives machine-precision values of L[g] (unit tests below).

**Eigenproblem (eig_spectral.py).** Collocate at Chebyshev–Lobatto points
Z_i = (1−cos(πi/(N+1)))/2, i=1..N+1 (includes Z=1, where L[φ_j'''] is finite for this
basis): A_ij = φ_j(Z_i), B_ij = −L[φ_j'''] (Z_i), solve the generalized problem
B c = λ A c with scipy.linalg.eig. Positive real λ = decay rates; the eigenvector gives
the mode shape.

## Derivation in brief

- Chemical potential μ = −σ_gb Ω ⇒ flux j = −(δ_gb D_gb/kTΩ)∂_z μ = (δ_gb D_gb/kT) ∂_z σ_gb.
- Mass conservation 2∂_t u = −Ω ∂_z j ⇒ ∂_t u = −(δ_gb D_gb Ω/2kT) ∂²_z σ_gb.
- Elasticity: σ_gb = σ₀ − (E*/2π) L[u'] (kernel S given in the problem). Differentiate in t
  and insert ∂_t u' = −(δ_gb D_gb Ω/2kT) σ_gb''':

      ∂_t σ_gb = (E*δ_gb D_gb Ω/4πkT) L[σ_gb''']  =  (h³/t₀) L[σ_gb''']  ⇒  ∂_τ Σ = L[Σ'''].

- Modes Σ = f(Z) e^{−λτ}: **λ f = −L[f''']**, f(0)=0, f'(1)=f''(1)=0.
- Steady state (τ→∞): Σ→0 requires L[Ã]=1 for the rescaled dislocation density
  Ã = (E*/2πσ₀) ∂_Z U. With the crack Ansatz Ã = R(Z)/√(1−Z) and the crack-opening
  relation u ∝ K√r, K/(σ₀√(πh)) = −π√2 R(1).

## (a) How σ_gb(z,t) evolves; λ₁, λ₂

The uniform initial traction first relaxes near the surface (where Σ(0)=0 is enforced and
plating starts); higher modes die out extremely fast (λ₂/λ₁ ≈ 22.5, λ₃ = 1357.4,
λ₄ = 3845.2 — the spectrum is real and rapidly growing, roughly the ((2n−1)/2)³ ladder).
By τ ≳ 0.02 the shape has locked onto the slowest mode and the whole boundary decays as a
single exponential:

    σ_gb(z,t) ≈ σ₀ · f₁(z/h) · exp(−λ₁ t/t₀),   λ₁ = 11.7426681,  λ₂ = 264.63422.

Slowest mode shape (normalized f₁(1)=1): zero at the surface, rises steeply over the outer
third of the boundary, maximal and flat (f₁'(1)=f₁''(1)=0) at the tip:

    z/h : 0.1    0.2    0.3    0.4    0.5    0.6    0.7    0.8    0.9    1.0
    f₁  : 0.2132 0.3999 0.5603 0.6942 0.8020 0.8843 0.9426 0.9788 0.9962 1.0000

So the last traction to disappear is the deep part near the substrate.

## (b) The t→∞ limit: a traction-free edge crack

The wedge grows monotonically and σ_gb → 0 everywhere: the grain boundary becomes a fully
open, traction-free, crack-like slit of depth h, i.e. the classic edge crack in a
semi-infinite plate. Solving L[Ã] = 1 with the crack-singular basis gives

    K(t→∞) = 1.1215223 · σ₀√(πh),

which agrees with the Wigglesworth/Tada edge-crack constant 1.1215 to all quoted digits — a
non-trivial consistency check (our normalization and kernel reproduce a classical,
independently known number). Note that at every finite time the tip traction is bounded
(flat profile enforced by the BCs); the √-singular crack field is attained only as τ→∞.

## (c) Relaxation time scale

Asymptotic e-folding time **1/λ₁ = 0.0852 t₀**. Operationally (time-domain simulation from
the uniform initial state): the spatial mean traction halves by τ ≈ 0.085, falls to 1/e of
its early value by τ ≈ 0.11, and the tip traction is below 1% of σ₀ by τ ≈ 0.4. So the
traction relaxes on a scale of ~0.1 t₀ (reported value: relaxation_time_t0 = 1/λ₁).

## Checks

All run in this folder; `python3 test_kernel.py`, `python3 eig_poly.py`,
`python3 steady.py`, `python3 sim.py`, `python3 solve.py`.

1. **Quadrature unit tests** (`test_kernel.py`). Compared the Gauss–Legendre PV machinery
   against (i) the closed form S₀(Z) = L[1], (ii) an independent subtraction-based reference
   for g = (1−Z')^(−1/2) including the finite-part limit at Z=1 (= −2 + regular part), and
   (iii) a polynomial g. Output (pasted):

        test g=1 vs S0 analytic:      diffs between -1.4e-07 (Z=0.001) and 1e-16
        test g=(1-Z')^{-1/2}:          diffs ~1e-11 … 3e-10 ; Z=1: diff=-9.19e-14
        test g=Z'^2:                   diffs ~1e-16

   Could have failed if: PV subtraction or image terms wrong, endpoint finite-part mishandled
   → diffs ≫ 1e-10.

2. **Convergence in basis size N** (fractional basis, M=768). Pasted:

        N= 40  lam1=11.7426686106  lam2=264.6343022204  f(0.5)=0.8019565005
        N= 48  lam1=11.7426682058  lam2=264.6342579604  f(0.5)=0.8019565024
        N= 56  lam1=11.7426681012  lam2=264.6342374644  f(0.5)=0.8019565029
        N= 64  lam1=11.7426680757  lam2=264.6342283101  f(0.5)=0.8019565030
        N= 72  lam1=11.7426680713  lam2=264.6342241574  f(0.5)=0.8019565030

   λ₁ stable to ~9 digits, λ₂ to ~7 (limit ≈ 264.63422), f(0.5) to ~9. Could have failed if
   values drifted or oscillated.

3. **Independent discretization** (`eig_poly.py`): purely polynomial basis
   φ₀=1−(1−Z)³, φ_j = Z(1−Z)³T_{j−1}(2Z−1), interior Chebyshev collocation. Converges
   algebraically (the (1−Z)^(5/2) tip behavior is not representable) but unambiguously
   toward the same numbers; power-law extrapolation gives λ₁ → 11.7420 and λ₂ → 264.63
   vs 11.74266807 / 264.63422. Pasted:

        N= 10   lam1=12.57770369   lam2=245.03346022
        N= 24   lam1=12.14015828   lam2=256.62664801
        N= 44   lam1=11.96895727   lam2=260.30980560
        N= 80   lam1=11.87174242   lam2=262.24074281   (still approaching 11.7427 / 264.63)

   Could have failed if it converged elsewhere → shared-formula or assembly bug.

4. **Spectrum reality / residuals.** At N=30: all 31 generalized eigenvalues finite; the
   resolved spectrum is real (first 8: 11.7427, 264.6344, 1357.447, 3845.184, 8327.4,
   15382.3, 25596.8, 39554.7); the only genuinely complex eigenvalues are discretization-tail
   artifacts with |λ| ≥ 4.49×10⁵, three orders of magnitude above λ₂, so no slow complex mode
   is missed; mode-1 generalized residual 8.2×10⁻¹². Could have failed: a complex pair or a
   smaller real eigenvalue below λ₂, or a large residual.

5. **Time-domain simulation** (`sim.py`): exact exponential integration of
   dc/dτ = −A⁻¹Bc from the uniform initial state (least-squares projection). Pasted:

        late-time decay rate from tip traction: 11.742669  (lambda1 spectral = 11.742668)
        tau=0.3: sigma(0.5)/sigma(1) = 0.80195650   (mode1 spectral = 0.80195650)

   Dynamics reproduces both λ₁ and the mode shape. Could have failed: wrong sign/operator in
   the evolution, or an eigenpair that is not actually a decay mode.

6. **Steady state vs classical fracture** (`steady.py`): K_inf_norm converges
   (Nc=32→96: 1.1215221959 → 1.1215222524 → 1.1215222555 → 1.1215222554) to the
   Wigglesworth/Tada edge-crack value 1.1215223. Could have failed: any error in the E*/2π
   normalization or kernel would shift this number; agreement to 7 digits validates both.

## Questions / assumptions

- The kernel S(z,ζ) is taken as stated (it already contains the free-surface image terms);
  no additional image construction was needed.
- `mode1_shape_half`: the slowest mode's traction is maximal at the tip, so the requested
  normalization f(1)=1 is natural; value at z/h=0.5 reported (0.80195650).
- `relaxation_time_t0` is reported as the asymptotic e-folding time 1/λ₁ = 0.08516; the bulk
  relaxation is effectively complete by ~0.4–0.5 t₀ (see (c)).
- The initial uniform traction is incompatible with Σ(0,t)=0 at t=0 (a jump); the BCs are
  enforced for all t>0 and the initial condition is projected modally. This affects only the
  very-early transient, none of the reported numbers.
- Uncertainty: λ₁ and f(0.5) are good to ~9 shown digits; λ₂ to ~7 (264.63422 ± 2e-5);
  K_inf_norm to ~1e-9.

## Files

- `kernel.py` — L-operator (PV quadrature, smoothing substitution), analytic pieces.
- `eig_spectral.py` — fractional-power spectral eigensolver (main method).
- `eig_poly.py` — independent polynomial-basis eigensolver (check).
- `steady.py` — t→∞ crack problem, K_inf.
- `sim.py` — time-domain simulation (check).
- `test_kernel.py` — quadrature unit tests.
- `solve.py` — one-command reproduction; writes and prints `results.json`.

