# ANSWER — stress relaxation by diffusion into a grain boundary of a thin film

## Summary of answers

| quantity | value |
|---|---|
| `lambda1` | **11.7427** |
| `lambda2` | **264.63** |
| `mode1_shape_half` = f₁(z/h=0.5), f₁(1)=1 | **0.8020** |
| `K_inf_norm` = K(t→∞)/(σ₀√(πh)) | **1.1215** |
| `relaxation_time_t0` | **0.0852 t₀** (= 1/λ₁; 1/e at ≈0.10 t₀, 95% relaxed by ≈0.27 t₀) |

Reproduce with **`python3 solve.py`** (runs in ~2 s, prints every number above plus all
validation checks, and rewrites `results.json`).

## (a) How σ_gb(z,t) evolves

With x = z/h, τ = t/t₀, s = σ_gb/σ₀, w = (E\*/2πσ₀)∂u/∂z the problem becomes

- elasticity: s = 1 − A[w],  A[w](x) = PV∫₀¹ Ŝ(x,ξ) w(ξ) dξ,  Ŝ = 1/(x−ξ) − 1/(x+ξ) − 2ξ(x−ξ)/(x+ξ)³
- diffusion:  w_τ = −s_xxx  (from 2u_t = −Ω(δ_gbD_gb/kT)σ_zz)
- BCs: s(0,τ)=0, s_x(1,τ)=s_xx(1,τ)=0; IC s(x,0)=1.

Right after t = 0⁺ the surface reservoir (μ=0) pulls the traction to zero near z = 0
(a boundary layer of the coupled elastic-diffusion problem); the relaxed region grows
and the traction profile evolves through a sequence of shapes. At late times the
evolution is a sum of decaying eigenmodes, s(x,τ) ≈ Σ cₙ fₙ(x) e^{−λₙτ}, with
**λ₁ = 11.7427, λ₂ = 264.63, λ₃ = 1357.4, …** (the spectrum is real, positive,
discrete; λₙ grows ≈ n^{4.5}, the asymptotics expected for this 3rd-order nonlocal
operator with square-root endpoint singularities). The slowest mode shape (normalized
f₁(1)=1) rises smoothly from 0 at the free surface to 1 at the interface:
f₁(0.1)=0.213, f₁(0.5)=0.802, f₁(0.9)=0.996; near the surface f₁ ≈ 2.26 x, near the
tip f₁(x) − 1 ∝ (1−x)^{5/2}.

## (b) The t → ∞ limit

Diffusion stops only when μ is uniform, and μ(0)=0 forces σ_gb ≡ 0: the grain
boundary becomes a **traction-free, atom-filled slit — a "diffusion wedge" that is
elastically exactly a surface crack of depth h** (the plated material carries no
stress, so the filled crack ≡ traction-free crack). By superposition this is an edge
crack in a half-space (film and substrate have identical moduli) loaded by uniform
pressure σ₀ on its faces, so

  K(t→∞) = 1.1215 σ₀√(πh),

the classical edge-crack value, which my own crack-limit computation reproduces
(1.12152 from the wedge-tip dislocation-density coefficient, 1.12152 from the
extrapolated traction ahead of the tip). During the transient the tip field is milder
than a crack's; the crack-like 1/√ singularity (and K) grows as 1 − O(e^{−λ₁t/t₀})
toward this limit.

## (c) Relaxation time

The traction relaxes on the time scale **t/t₀ ≈ 1/λ₁ = 0.085**, i.e.
t_relax ≈ 0.0852 t₀ (e-folding of the dominant mode). Including the early transient,
s at mid-boundary falls to 1/e by t ≈ 0.10 t₀ and to 5% by t ≈ 0.27 t₀. In physical
units, t_relax = 0.0852 · 4πkTh³/(E\* δ_gb D_gb Ω).

## Method

1. **Elasticity kernel as an operator**: discretized A by product integration with the
   PV singularity subtracted,
   PV∫g(ξ)/(x−ξ)dξ = ∫[g(ξ)−g(x)]/(x−ξ)dξ + g(x)ln(x/(1−x)), Gauss–Legendre quadrature.
2. **Final state**: solved 1 − A[w] = 0 with w = v/√(1−x) (the crack's square-root
   dislocation density at the tip), least-squares Chebyshev collocation; K extracted
   from the tip density coefficient (K = π√2 E\*σ₀^{-1}√h-normalized coefficient) and,
   independently, from the traction ahead of the tip.
3. **Eigenproblem**: separable modes obey A[f‴] = −λf with f(0)=0, f′(1)=f″(1)=0. I
   solved it in terms of g = f‴ in the basis g = x(1−x)^{−1/2}·Chebyshev, with the
   three traction BCs built in analytically (f is reconstructed from g by exact
   integration: f = (∫tg)x − (∫g/2)x² + ∫₀ˣ(x−ξ)²g/2). All quadratures use the
   substitution ξ = 1−t² which makes every integrand smooth. Generalized eigenvalue
   problem eig(A, −F).
4. **Relaxation time**: from λ₁; cross-checked by a 10-mode reconstruction of the
   evolution.

## Checks performed (all pasted from `python3 solve.py` output)

- **Kernel quadrature vs exact formula** A[1](x) = ln(x²/(1−x²)) + 6x/(1+x) +
  2ln(1+x) − 2x²/(1+x)² − 4 − 2ln x (derived by hand, A[1](x→0⁺) = −4):
  ```
  x=0.01: numeric -3.920789451921  exact -3.920789451909
  x=0.5:  numeric -1.123609933554  exact -1.123609933554
  x=0.9:  numeric  1.337790779720  exact  1.337790779720
  ```
  Could have failed if the PV split or image terms were wrong; it pins the kernel to 12 digits.
- **Crack limit** (validates kernel + solver + my K normalization simultaneously):
  ```
  K_inf/(sigma0 sqrt(pi h)) from wedge-tip density coefficient : 1.121522
  K_inf/(sigma0 sqrt(pi h)) from traction ahead of tip (extrap): 1.121516
  handbook edge-crack in half-space, uniform face pressure      : 1.1215
  emergent w(0) = -0.249998 (exact -0.250000)
  ```
  The emergent w(0) = −1/4 is an exact analytic consequence of A[c](0⁺) = −4c
  (i.e. of σ(0)=0). A wrong kernel sign/normalization would have broken both agreements.
- **Eigenvalue convergence** (g-formulation; could have failed if the assumed
  endpoint singularities x·(1−x)^{−1/2} were wrong — then values would drift with Nv):
  ```
  Nv= 16:  11.7426636   264.63391
  Nv= 24:  11.7426678   264.63413
  Nv= 48:  11.7426681   264.63422
  Nv= 64:  11.7426681   264.63422
  ```
  (also unchanged for quadrature sizes Mt = 200…1600, and with uniform instead of
  Chebyshev collocation points: 11.7427, 264.6 — could have failed if collocation were
  under-resolved.)
- **Independent discretization** (differentiation-based s-space pencil; converges ~N⁻²
  because of the (1−x)^{5/2} tip behavior of f, extrapolated):
  ```
  N=40: lam1=11.77468 lam2=265.3428  f(0.5)/f(1)=0.80238
  N=60: lam1=11.75682 lam2=264.9493  f(0.5)/f(1)=0.80215
  N=96: lam1=11.74818 lam2=264.7574  f(0.5)/f(1)=0.80203
  N^-2 extrapolated: lam1=11.7426 (g-form 11.7427); lam2=264.634 (g-form 264.634)
  ```
- **Pointwise eigenpair residual** |A[g] + λ₁f| on an independent 2000-node grid:
  ≤ 7×10⁻⁷ everywhere (typical 1×10⁻⁷), vs terms of size ~6.
- **Late-time decay rate from the mode expansion** (could have failed if λ₁ were a
  spurious mode or had zero amplitude):
  ```
  fitted d/dtau ln s at x=0.25: -11.7427   (-lambda1 = -11.7427)
  fitted d/dtau ln s at x=0.5 : -11.7427   (-lambda1 = -11.7427)
  fitted d/dtau ln s at x=0.9 : -11.7427   (-lambda1 = -11.7427)
  s(0.5) falls to 1/e: tau = 0.1005 t0 ;  to 5%: tau = 0.2705 t0
  ```
- **Spectral asymptotics**: successive eigenvalue ratios match λₙ ∝ n^{4.5}
  (11.74, 264.6, 1357, 3845, 8327, 15382 → ratios 22.5, 5.13, 2.83, 2.17, 1.85 vs
  theoretical 22.6, 5.14, 2.83, 2.17, 1.85), as expected for a 3rd-order operator
  with ½-power endpoint singularities — an internal consistency test of the spectrum.

## Questions / assumptions

1. **σ_gb(0,t)=0 vs the elasticity relation**: S(0,ζ)=0 pointwise, so the BC is
   meaningful only as a limit; analytically A[c](0⁺) = −4c makes σ(0)=0 ⟺ w(0)=−1/4
   (verified to emerge from the crack solve). For eigenmodes it becomes f(0)=0.
2. **Final-state K**: I assume K is defined in the t→∞ limit at the tip z=h from the
   1/√(z−h) field just ahead of the tip (the GB itself becomes traction-free);
   equals the handbook edge-crack value 1.1215σ₀√(πh), which I also recompute.
3. **Relaxation time** is reported as 1/λ₁ (e-folding); "95% relaxed" would be
   ≈ 3.2/λ₁ ≈ 0.27 t₀. The problem does not define the measure, so both are given.
4. The IC σ(z,0)=σ₀ is incompatible with σ(0)=0 at the single point z=0 (instantaneous
   surface layer); this affects only the short-time transient, not the λ's.
5. Everything assumes the problem's own model (isotropic, identical film/substrate
   moduli ⇒ half-space; flat surface; no sliding/diffusion at z=h). Surface grooving
   and the exact atomic structure at the GB/surface junction are neglected as stated.

## Files

- `solve.py` — single-command reproducer (prints everything, writes results.json).
- `results.json` — the reported quantities.
- `gbcore.py`, `crack.py`, `eig.py`, `eig_g.py`, `dbg*.py`, `check1.py` — development
  and validation scripts referenced above.
- `PROCESS.md` — step log.
