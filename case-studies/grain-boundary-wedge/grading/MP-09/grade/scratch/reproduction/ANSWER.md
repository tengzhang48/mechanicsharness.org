# MP-09: Stress relaxation by diffusion into a grain boundary — ANSWER

**Reproduce all numbers with one command:**

```
python3 solve.py
```

(prints every value in `results.json`; runtime ≈ 45 s; files: `solve.py`,
`sblib.py`, `gblib.py`; log of the run used for this answer: `run.log`)

---

## Model and reduction

With x = z/h, τ = t/t₀, t₀ = 4πkTh³/(E*δ_gbD_gbΩ), and σ in units of σ₀,
combining elasticity, μ = −σΩ, Fick's law and mass conservation gives the
linear integro-differential equation

    σ_τ = A σ,   (Aσ)(x) = PV ∫₀¹ S(x,ξ) σ_ξξξ(ξ) dξ,
    S(x,ξ) = 1/(x−ξ) − 1/(x+ξ) − 2ξ(x−ξ)/(x+ξ)³,

with BCs σ(0)=0 (source/sink at the free surface), σ'(1)=σ''(1)=0 (no
diffusion/sliding at the substrate), σ(x,0)=1. Modes f(x)e^(−λτ) satisfy
A f = −λ f. The operator is **not** self-adjoint and has a genuine
eigenvalue spectrum; ker(A) = span{1} is excluded by the BCs.

**Key structural fact:** eigenfunctions behave as (1−x)^(7/2) at the tip
(f''' ~ (1−x)^(1/2), a weak singularity). All methods must resolve this.

## (a) Evolution, eigenvalues, slowest mode

The uniform initial traction first goes through a short **multiscale
transient**: the incompatibility at the mouth (σ(0,0⁺)=1 vs BC σ(0)=0) is
resolved instantly, the near-tip region adjusts on the fast scales 1/λ₂,
1/λ₃, … (the traction near the tip briefly *rises* above σ₀, to ≈1.3σ₀ at
τλ₁ ≈ 0.1), and then the slowest mode takes over:

    σ_gb(z,t) ≈ c₁ σ₀ f₁(z/h) e^(−λ₁ t/t₀),   c₁ ≈ 1.5 (measured 1.53–1.55),

with f₁(0)=0, f₁(1)=1, f₁ ≥ 0, single-signed (no interior zero). The
traction relaxes monotonically to **zero everywhere** — at late times the
grain boundary carries no load.

Numerical values (converged to ~5e-7 across method, resolution and
eigensolver tolerances):

    λ₁ = 11.742668
    λ₂ = 264.6342        (λ₂/λ₁ ≈ 22.5: the transient is over quickly)
    f₁(0.5) = 0.8019565  (normalization f₁(1) = 1)

## (b) Steady state: the boundary becomes an edge crack; K

Since σ_gb → 0, the final state has a **traction-free grain boundary**: the
plated wedge of material (which stops growing because σ_gb → 0 shuts off
the chemical-potential gradient) holds the boundary open exactly like a
**classical edge crack of length h in a half-space under remote tension
σ₀**. Therefore

    K(t→∞) = 1.1215 σ₀ √(πh)   — K_inf_norm = 1.1215

(more precisely 1.12152; the accepted handbook value 1.1215 for a single
edge crack in a semi-infinite plate). Verified by solving the crack
integral equation PV∫S v̂ = 1 directly and extracting K by two independent
routes, which also reproduce the classical mouth opening
u(0)E*/(σ₀h) = 2.9086 (our COD_mouth = 2u(0) = 5.8172) and satisfy the
universal near-tip COD–K relation exactly.

## (c) Relaxation time

    t_relax = 1/λ₁ = 0.085160 t₀

i.e. the traction decays by e every 0.085 t₀ (half-life of the exponential
tail ln2/λ₁ = 0.059 t₀; the tail dominates after τ ≳ 2/λ₁ ≈ 0.17 t₀).
Time-domain integration confirms the fitted decay rate matches λ₁ and that
max|σ| → 0.

---

## Methods

**Spectral eigensolver (primary).** Tip-resolving coordinate s = √(1−x);
basis φ_k(s) = s^k − 1, k = 5,…,N+4, which enforces all three BCs (the
k≥5 start also removes the spurious ker-mode s⁴, whose third-derivative
coefficient c₄ = 0). Collocation at Chebyshev points in s. The kernel
action is *exact*: with dξ = 2r dr (ξ = 1−r²),

    (Aφ_k)(s) = 2 c_k R_{k−5}(s),  c_k = −k(k−2)(k−4)/8,
    R_m(s) = PV∫₀¹ S̃ r^m dr = M_m + G_m,

M_m by the forward-stable recurrence M_m = 1/(m−1) + s²M_{m−2} (analytic
PV seeds), G_m (smooth image terms) by 4200-point Gauss–Legendre. The
generalized pencil (Θ, F) is numerically singular for N ≳ 32 (cond(F)
~ 1e19: columns s^k−1 become indistinguishable), which destroys the QZ
eigenvalues if used directly. **Fix: truncated-SVD projection** — SVD
F = UΣVᵀ, keep singular values > 10⁻¹²σ_max, QZ the reduced pencil
(UᵀΘV, Σ). λ₁ then varies < 5e-7 over N = 60…110 and tol = 1e-11…1e-13.

**Finite-difference cross-checks.** Two independent panel discretizations
(A: cell-centered with quadratic mouth ghost; B: node-based), kernel
integrated analytically over panels (exact PV). Their λ₁ approach from
opposite sides and bracket the spectral value; Richardson limits agree with
the spectral value to ~2e-5; f₁(0.5) is likewise bracketed.

**Crack solve (part b).** v̂(ξ) = w(ξ)/√(1−ξ), w(ξ) = Σ c_m(1−ξ)^m,
moment collocation (T v̂)(s) = 2R_{2m}(s) with the same exact moments —
well-conditioned (c = O(1)). K from the tip coefficient,
K/(σ₀√(πh)) = −π√2·c₀ (identity from the near-tip fields), and
independently from the outside stress σ_yy(x) = 1 − ∫S v̂,
K = lim_{ε→0} √(2ε) σ_yy(1+ε) (evaluated with exact external moments,
Richardson in √ε). An earlier version of this check used the
monomial-in-ξ basis for post-processing, which mislabeled c₀ as Σc_m and
gave a spurious 1.1107 — found and fixed (see PROCESS.md); the two routes
now agree to 1.2e-5.

**Time domain.** Implicit Euler on the FD-A operator (N = 512, dt = 2e-4,
single LU factorization), σ(0) = 1: confirms σ → 0, gives c₁ and the
decay-rate check.

## Checks (what was compared, output, how it could have failed)

1. **FD operator vs exact A f.** f = x²(1−x)³ satisfies all BCs and
   A f = −18T₀ + 72T₁ − 60T₂ exactly (T_m = moments of the kernel action on
   monomials, from gblib). Output:
   `FD-A N=1024: max|L f - Af| = 4.908e+00 (interior 1.908e-01)`,
   `FD-B N=1024: max|L f - Af| = 2.041e+01 (interior 2.368e-01)`
   (errors concentrated in 2–4 cells at the ends where the FD closure
   smears the tip singularity; the interior error is O(dx)). Could have
   failed: any error in the kernel, panel integration, or BC closures would
   show O(1) interior error.
2. **Spectral convergence + off-grid residual** of the *continuous*
   eigen-equation at points never collocated (strongest test of the
   eigensolver; would expose spurious eigenvalues from the singular
   pencil):
   `N= 60: lam1=11.7426678465 lam2=264.6341789549 f(1/2)=0.8019564986 resid 1.7e-06`
   `N=110: lam1=11.7426682103 lam2=264.6342256450 f(1/2)=0.8019565032 resid 1.0e-06`
3. **FD bracketing** (independent discretizations, opposite approach
   directions — would fail if the spectral answer were contaminated by the
   pencil singularity or a wrong BC):
   `N= 1024: FD-A lam1=11.73407269 | FD-B lam1=11.75127120`
   `Richardson limits: FD-A ~ 11.74268899 ; FD-B ~ 11.74267930 ;
    spectral 11.74266821` (bracket width ~1e-5; the extrapolations carry
   ~2e-5 residual error, the N=1024 values strictly bracket the answer).
   f(1/2): FD-A 0.80172 / FD-B 0.80219 bracket 0.80195650.
4. **Time evolution vs eigenvalues.** Implicit Euler from σ=1: late-time
   fitted decay rate 11.7127 vs FD-A's own discrete λ₁(N=512) = 11.7255 and
   Euler rate factor (1−λdt/2) = 0.99883 → 11.7118 expected — agreement to
   1e-3. A wrong spectrum (complex pairs, positive rates) or missing ker
   removal would show here as non-decay or a wrong slope. Also shows
   σ → 0 (`final max|sigma| at tau*lam1=6.0: 3.853e-03`, still decaying
   ∝ e^(−λ₁τ)) — a nonzero steady traction would have contradicted the
   edge-crack conclusion.
5. **K, two independent routes + classical value:**
   `crack M=24: c0=-0.252431220520 K_tip=1.1215223627 K_out=1.1215340709
    COD_mouth=5.81719168`
   vs classical 1.1215 (Tada) and mouth opening 2×2.9086 = 5.8172.
   M-convergence 10→24: K_tip = 1.1214990 → 1.1215224 (last-digit drift
   ~2e-7). Could have failed: wrong kernel factors/conventions would move
   K off 1.1215; the two extraction routes (tip coefficient vs outside
   stress) use disjoint functionals and only agree if the solve is correct.
   The near-tip COD–K consistency relation COD ≈ (8K/E*)√(r/2π) was also
   verified against the solution (coefficients 6.3435 vs 6.3443, 1e-4).

## Uncertainties / limitations

- λ₁, λ₂, f(0.5): accurate to ~5e-7 (spectral convergence spread); FD
  brackets consistent at 2e-5. Higher modes available but not needed.
- K: computed 1.121522(2); the outside-stress route sits 1.2e-5 high
  (extrapolation residual). Both equal the classical 1.1215 to the
  accuracy with which that constant is known (1.12152).
- c₁ ≈ 1.53σ₀: measured from FD time evolution with the spectral mode
  profile; the slow creep of the fitted c1 (1.532→1.554 over
  τλ₁ = 0.8→6) is the FD operator's 0.15% eigenrate offset times
  e^(λ₁τ) back-extrapolation; the quoted value is the honest
  late-asymptotic amplitude of the tail, defined through the evolution
  (the initial data σ=1 violates the mouth BC, so a purely modal
  projection of the initial condition is not well-defined at t = 0).
- Physics assumed as stated: linear isotropic elasticity, no grooving,
  flat surface, fast surface diffusion, no interface diffusion/sliding.

## Questions / assumptions

- None blocking. Convention notes: the wedge opening is 2u (as given), so
  the crack-opening displacement is COD = 2u; the K normalization
  σ₀√(πh) is the standard one; E* = E/(1−ν²) as stated (plane strain).

## Files

- `gblib.py` — kernel S, antiderivatives, exact moments in x (validation).
- `sblib.py` — s-coordinate spectral pencil, SVD-projected QZ eigensolver,
  off-grid residual, crack moment collocation, external moments, K/COD.
- `solve.py` — runs everything, prints all values, writes `results.json`.
- `run.log` — output of the run used for this answer.
