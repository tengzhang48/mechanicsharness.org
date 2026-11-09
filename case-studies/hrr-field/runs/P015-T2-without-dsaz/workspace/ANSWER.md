# Answer — tensile sharp notch, plane strain, power-law (n = 4)

## Requested numbers

| quantity | value |
|---|---|
| **p**  (σθθ(r,0) ∝ r⁻ᵖ) | **0.08764** (0.0876) |
| **R**  = σrr(r,0)/σθθ(r,0) on θ = 0 | **0.3330** |

The opening amplitude is taken positive, i.e. σθθ(r,0) > 0 (tensile) as r → 0.

## Method

**Constitutive model.** Near a sharp notch tip the plastic term of the
Ramberg–Osgood law dominates the elastic term (because σ_e → ∞, and
ε^p ∝ σ_eⁿ with n = 4 > 1), so the leading small-strain field obeys the
incompressible power-law law

  ε_ij = (3/2) B σ_e^{n−1} s_ij ,   n = 4 ,   (ε_kk = 0).

Plane strain (ε_zz = 0) ⇒ s_zz = 0 ⇒ σ_zz = ½(σ_rr + σ_θθ); the inplane
field is deviatoric with active stress ∇·σ = 0.

**Self-similar (HRR-type) ansatz.** For a power-law material the dominant
field is separable: σ_ij = r^{−s} Σ_ij(θ), equivalently an Airy function
Φ = r^{2−s} f(θ), and displacement u_i = r^{1−sn} U_i(θ).

**Geometry / BCs (this task).** Material occupies −2π/3 < θ < 2π/3, so the
half-angle is α = 2π/3 (missing wedge 120°). Mode-I opening under remote
tension in y:
* symmetry line θ = 0:  u_θ = 0,  σ_rθ = 0;
* traction-free flanks θ = ±α:  σ_θθ = 0,  σ_rθ = 0.

**Eigenvalue problem.** The constitutive law + incompressibility + the two
equilibrium equations reduce (unknowns U(θ), V(θ), pressure P(θ), exponent s)
to a nonlinear BVP that we solve directly in the **displacement
formulation**:

  inc :  (2 − s n) U + V′ = 0
  dev :  d_rr = (1−sn)U ,  d_θθ = −(1−sn)U ,  d_rθ = ½(U′ − sn V)
         K = e2^((1−n)/2n) ,  s_ij = (2/3) K d_ij ,  e2 = (2/3)(Σ d² + 2 d_rθ²)
  S_ij = s_ij + P δ_ij
  eq_r :  −sn S_rr + (S_rr − S_θθ) + S_rθ′ = 0
  eq_θ :  (2 − sn) S_rθ + S_θθ′ = 0

solved with `scipy.integrate.solve_bvp` (unknown parameter s), amplitude
gauged by U(0) = 1.

## Checks

1. **n = 1 validation (exact).** At n = 1 the model reduces to the
   incompressible linear problem; the solver returns **s = 0.38427** for
   α = 2π/3, which equals **1 − λ** of the Williams mode-I root
   λ = 0.61573 (sin 2λα + λ sin 2α = 0). Continuing in α to the crack
   (α → π) gives **s → 0.5000**, the classical HRR/Williams crack value
   1/2 at n = 1. This pins the physical mode-I branch.

2. **Branch continuation in n.** Following that mode-I branch from n = 1
   to n = 4 (steps Δn = 0.25, residual ≤ 5·10⁻⁹ throughout) gives a smooth
   monotone sequence
   n: 1 → s 0.38427;  2 → 0.18481;  3 → 0.11970;  **4 → 0.087639**.

3. **Finite energy selects the branch.** Two real mode-I eigenvalues exist
   at n = 4: s = 0.08764 and s = 0.41236 (they satisfy s₁ + s₂ = 0.5). The
   strain-energy density scales as W ∝ σ_e^{n+1} ∝ r^{−s(n+1)}; the energy
   in a neighbourhood of the tip is finite iff s(n+1) < 2, i.e. s < 0.4.
   The mode with s = 0.41236 has infinite energy and is excluded; the
   finite-energy leading symmetric opening field is **s = 0.087639**.

4. **Direct evaluation of R.** From the converged field at n = 4,
   α = 2π/3 (grid M = 800, collocation residual 6·10⁻⁹):

       σ_rr(0) = 0.5365 ,  σ_θθ(0) = 1.6110 ,  R = 0.3330 .

## Conclusion

    p = 0.08763859   (tolerance 0.0003)  → report 0.0876
    R = 0.33301436   (tolerance 0.002 )  → report 0.3330

**p = 0.08764**,  **R = 0.3330**.

### Grid / tolerance convergence of the final run (conv.py)

| back grid Mb | fine grid Mf | solve_bvp tol | p | R |
|---|---|---|---|---|
| 200 | 400  | 1e-10 | 0.08763859 | 0.33301436 |
| 400 | 800  | 1e-11 | 0.08763859 | 0.33301436 |
| 500 | 1000 | 1e-12 | 0.08763859 | 0.33301436 |
| 600 | 1500 | 1e-12 | 0.08763859 | 0.33301436 |

Identical to 8 decimals → both requested numbers are resolved far inside the
stated tolerances.