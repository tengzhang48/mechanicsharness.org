# Method

## Formulation

The crack plane is y = 0, the front is x = a(z), and the loading is fixed at infinity
(the far field of a straight front with uniform K's K₁⁰ > 0, K₂⁰ = 0.6 K₁⁰, K₃⁰ = 0.8 K₁⁰).
Everything is linear in the small amplitudes A/λ (Part A) and ε (Part B), so the response is
computed from the **first-order (in front-perturbation) operator**, evaluated per z-Fourier mode.

### 1. Crack-plane stiffness matrix (exact, symbolic-linear algebra)

For a displacement-discontinuity (DD) distribution B e^{i(px+qz)} on the plane y = 0 the
traction on that plane is t = M(p,q)B, with (μ=1)

- M_yy = −r/(2(1−ν)),  r = √(p²+q²)
- in-plane (x,z) block: −[(1−ν) r I + ν k kᵀ/r] / (2(1−ν)),  k = (p,q)

Derivation (`m_stiffness.py`): half-space solutions u = (a + yb)e^{−ry}e^{i(px+qz)}
(b = c(ip,−r,iq), c fixed by Navier), glued across y = 0 by the y-/xz-reflections, which
impose the sector mixed conditions (opening sector: t_x = t_z = 0; shear sector: t_y = 0).
Verified: closed forms match the exact 3×3 sector solves to 1e-16; symmetry; isotropy
(M(Rk) = R M Rᵀ); Navier residual ~1e-16; q→0 reproduces the standard 2D crack-plane
constants −μ|p|/(2(1−ν)) (modes I, II) and −μ|p|/2 (mode III).

### 2. Per-mode first-order problem (scale-invariant; solved at |k| = 1)

For a front perturbation mode δâ e^{ikz}, the correction DD b(X) on the crack (X = depth
behind the front, in units 1/|k|) satisfies T[b] = 0 on X > 0 (Hadamard finite part),
b → 0 as X → ∞, with edge behaviour b ~ c₋ X^{−1/2} + c₊ X^{1/2}, where the input edge is
c₋ = D_i K_i⁰ δâ/2 (D₁ = D₂ = 2(1−ν)√(2/π), D₃ = 2√(2/π), plane-strain/anti-plane pairs).

Wiener–Hopf analysis (exact, `first_order.py`):

- **Mode I (y-sector)**: M_yy = cc√(p+i)√(p−i) factorises; the unique solution is
  b = c₋ e^{−X}X^{−1/2}, which is *exactly* traction-free on the crack, with ahead traction
  [1/(4(1−ν))] e^{−|x|}|x|^{−3/2}. Hence A₃₂ = 1/(√(2π)D₁) (exact), A₁₂ = −A₃₂, and
  **δK₁/K₁⁰ = −(|k|/2) δâ** — the Gao–Rice (1989) result, reproduced analytically.
- **Shear sector**: the b₁-columns (X^{−1/2}e^{−X} in each component) have closed-form
  Abel tractions: on the crack e^{−X}·(√π cc ν/√2)[[−1,i],[i,1]] (pole of M·Φ̂ at p = i,
  residue computation) and ahead (erfc/Γ closed forms), giving A₃₂ = −Q̂/2 — i.e.
  A₃₂ = 1/(√(2π)D₂), 1/(√(2π)D₃) *exactly*, off-diagonal A₃₂ ≡ 0. The remaining √X-edge
  part is solved by weighted least squares (collocation on X ∈ [10⁻⁴, 25], Laguerre basis
  e^{−X}X^{1/2}L_m, m ≤ 8; LS residual 1.7e-4, convergent in basis size).

Kernel evaluation (`kernels_x.py`): closed-form Bessel kernels
𝒦_yy = K₁(|t|)/(2π(1−ν)|t|), 𝒦_xx = [K₁/|t| + νK₀]/(2π(1−ν)), 𝒦_zz = [K₁/|t| − νK₀/(1−ν)]/2π,
𝒦_xz = −iν sgn(t)K₁/(2π(1−ν)), applied with fixed-panel Gauss–Legendre quadrature and
Hadamard finite parts (subtraction model matched to value/slope/curvature at the collocation
point; closed forms for the model integrals via Ei).

### 3. Resulting first-order operator (global basis, per mode n, σ = sign(n))

With γ = ν/(2−ν):

- δK₁^g = |k|δâ (−1/2) K₁⁰
- δK₂^g = |k|δâ [(−1/2+γ)K₂⁰ + iσγ K₃⁰]
- δK₃^g = |k|δâ [(−1/2−γ)K₃⁰ + iσγ K₂⁰]

The problem's K₂, K₃ are defined in the local (n, t) basis: K₂^prob = K₂^g − a'K₃⁰,
K₃^prob = K₃^g + a'K₃⁰-related, with a' = ∂a/∂z = iσ|k|δâ (mode-wise). This rotation was
verified against the exactly solvable inclined-front anti-plane problem
(K₂^prob = −ψK₃⁰ for a front inclined by ψ in a pure mode-III field).

## Checks run (all in this folder)

1. `m_stiffness.py`: exactness, symmetry, isotropy, 2D limits, Navier residual (all ≤ 1e-15).
2. `check_J.py`, `check_I2.py`: finite-part closed forms (Ei-based) to 12 digits.
3. `check_hadamard.py`, `check_arbiter.py`, `check_plus.py`, `check_b1ahead.py`:
   X-domain kernel application vs independent contour (branch-cut) quadrature — agree to
   ~1e-8 (ahead); b₁-column closed forms vs cut quadrature to 1e-12.
4. `check_b1on.py`: numerical on-crack traction of the shear b₁-columns reproduces the
   residue closed forms e^{−X}(√πccν/√2)[[−1,i],[i,1]].
5. Exact internal identities: A₁₂[X^{1/2}e^{−X}] = −A₁₂[b₁] for xx, zz (holds to 1e-5);
   D-consistency A₃₂ = 1/(√(2π)D_i) exact; off-diagonal A₃₂ = 0.
6. Extraction cross-checks: √x-series fit on the total ahead traction, algebraic large-p
   tail of MΦ̂ (Watson pairing, calibrated on b₁ and φ₀⁺), and a slip-route (c₊) extraction.
   The first two agree to 1e-5. The slip route was **rejected**: the exact inclined-front
   anti-plane test shows that displacement-jump asymptotics do not pair with the
   traction-defined K's at first order in the coupled (II+III) 3D problem.
7. `check_nu.py`: the closed form R = −½I + (ν/(2−ν))[[1,i],[i,−1]] is reproduced by the
   full pipeline at ν = 0.2, 0.25, 0.35 to ~5e-6.
8. Basis-size and collocation convergence (ML = 4…8, Ncol, X-range): R stable to ~1e-5.
9. Part B: equilibrium G(z) = Γ(z) verified pointwise to machine precision; f real,
   zero-mean; Fourier series converged (q^n, n ≤ 60).

## Part A

Single mode δa = A cos(2πz/λ): C₁ = −π, S₁ = 0; C₂ = 2πRe T₂, S₂ = −2πIm T₂,
C₃ = 2πRe T₃, S₃ = −2πIm T₃ with T₂ = (−1/2+γ)(0.6) + iγ(0.8) − i(0.8),
T₃ = (−1/2−γ)(0.8) + iγ(0.6) + i(0.6) (the ∓i0.8, ±i0.6 are the local-basis rotation).
C_G, S_G from δG via Irwin with the problem-basis K's.

## Part B

G⁰ = Γ₀ from the mean of Γ (g has zero mean), giving K₁⁰√((1−ν²)/(EΓ₀)) =
[(1−ν²)(1+0.6²) + (1+ν)0.8²]^{−1/2} = 0.663098. Per mode n:
κ₁ = −π|n|f̂_n, κ₂ = 2π|n|T₂(σ)f̂_n, κ₃ = 2π|n|T₃(σ)f̂_n, and the equilibrium
2[(1−ν²)(κ₁ + 0.6κ₂) + (1+ν)0.8κ₃]/S₀ = ĝ_n = q^{|n|} is solved algebraically for f̂_n
(S₀ = (1−ν²)(1+0.36)+(1.3)(0.64) = 2.0696). f, κ₁, κ₂, κ₃ evaluated at s = k/16;
κ₃ max found on a fine grid (peak is smooth; value/location stable to <1e-5).

## Assumptions

- First-order accuracy in A/λ and ε (as required by the task's O(A/λ), O(ε) definitions);
  the operator is exact within first-order theory.
- The far field is the straight-front K-field (the loading definition of §2); no T-stress
  or higher-order far-field terms (none are specified, and K-perturbations are independent
  of them at first order for this class of loadings).
- Crack faces traction-free and non-touching (small amplitudes); linear elasticity, ν = 0.3.

## Files

- `m_stiffness.py` — exact crack-plane stiffness matrix + verifications
- `kernels_x.py` — closed-form 2.5D kernels, finite-part application, ahead evaluation
- `first_order.py` — exact WH results, b₁ closed forms, shear-sector LS solve, extractions
- `solve_AB.py` — final assembly of Part A and Part B, writes `results.json`
- `check_*.py`, `diag_*.py` — verification scripts described above
