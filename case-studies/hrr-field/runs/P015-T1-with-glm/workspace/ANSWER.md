# Answer — HRR field, plane strain, n = 6.5

## Results

| quantity | value |
|---|---|
| **p** in σθθ(r,0) ∝ r⁻ᵖ | **0.133333** (= 2/15 = 1/(n+1); to 10 digits 0.1333333333) |
| **H** in σθθ(r,0) = H σ₀ [J/(α ε₀ σ₀ r)]ᵖ | **1.8936** (computed 1.89357) |
| **C** = U(R)/(J R) | **0.4573** (computed 0.45731) |

Supporting values: HRR integral I(n=6.5, plane strain) = 4.818379, hoop-stress angular
function σ̃θθ(0) = 2.335257 (normalization max_θ σ̃e = 1), energy coefficient
E_n = 1.436668. α, σ₀, ε₀, E set only the scales (ε₀ = σ₀/E = 2×10⁻³) and cancel
from all three dimensionless numbers; J and R likewise cancel — the answers are
independent of the load level and of R within the range where the leading term
dominates (checked explicitly at J = 5×10⁻⁵ J/m², R = 10⁻⁴ m).

## Method

The material law is a pure power law, so the near-tip field is the Hutchinson (1968)
plane-strain HRR field with Airy function Φ = A rˢ f(θ), s = (2n+1)/(n+1):

- σrr = A r⁻ᵖ (s f + f″), σθθ = A r⁻ᵖ s(s−1) f, σrθ = A r⁻ᵖ(1−s) f′, p = 2−s = 1/(n+1);
  σzz = (σrr+σθθ)/2 and σe² = (3/4)(σrr−σθθ)² + 3σrθ² (plane strain, εzz = 0 ⇒ szz = 0);
- f solves the 4th-order angular compatibility ODE (equilibrium + power-law
  compatibility) with f′(0) = f‴(0) = 0 (mode-I symmetry) and f(π) = f′(π) = 0
  (traction-free faces), obtained by shooting on f″(0);
- W = (n/(n+1)) αε₀σ₀ (σe/σ₀)^{n+1} ∝ r⁻¹; displacements u = αε₀(A/σ₀)ⁿ rᵖ(ũr, ũθ)
  follow from the strain–displacement relations (ũr from εrr, ũθ from εθθ).

Using the knowledge-unit solver (`knowledge/hrr.py`, rtol 10⁻¹¹, angular grid up to
32001 points) at n = 6.5:

1. **p** — the eigenvalue condition f(π) = f′(π) = 0 admits a non-trivial solution only
   at s = (2n+1)/(n+1) = 13/7.5 = 1.86667, so stresses ∝ r^{s−2} = r^{−1/(n+1)}
   ⇒ p = 1/7.5 = 2/15.
2. **H** — evaluating the full-circle J of the similarity field fixes the amplitude:
   J = α ε₀ σ₀ (A/σ₀)^{n+1} I, I = ∫_{−π}^{π}{...}dθ ≡ 4.818379 (Hutchinson's eq. (24)
   integral over −π..π). Hence σθθ(r,0) = σ₀[J/(αε₀σ₀I r)]ᵖ σ̃θθ(0) and
   **H = σ̃θθ(0)/I^{1/(n+1)} = 2.335257/4.818379^{2/15} = 1.89357** (normalization-free:
   rescaling σ̃ leaves H invariant since (n+1)p = 1).
3. **C** — since W ∝ r⁻¹ exactly, U(R) = R·(n/(n+1))·(J/I)·∫_{−π}^{π}σ̃e^{n+1}dθ, so
   **C = (n/(n+1)) ∫_{−π}^{π} σ̃e^{n+1} dθ / I = 0.45731** (= E_n/π, E_n = 1.436668),
   independent of R.

## Checks (verifications)

1. **J from first principles (strongest).** I reconstructed the *dimensional* field
   (σ₀ = 400 MPa, ε₀ = 2×10⁻³, α = 0.02, J = 5×10⁻⁵ J/m²), built the displacements by
   my own route (ũr from εrr, ũθ by quadrature of the εθθ relation — different from the
   code's shear-strain route), and evaluated the task's J integral over the full circle
   −π..π at R = 10⁻⁴ m with Cartesian displacement gradients ∂u/∂x computed (i) from my
   own closed-form polar→Cartesian formulas and (ii) by central finite differences of
   u(x,y). Results: J(direct)/J(given) = 0.9999999999 (analytic gradients) and
   1.0000000006 (FD gradients); analytic vs FD gradients agree to 4×10⁻¹⁷ absolute
   (2×10⁻⁸ relative). This validates the factor-2/full-circle convention of I, the
   displacement algebra (including rotation of the polar basis in ∂u/∂x) and the
   amplitude relation J = αε₀σ₀(A/σ₀)^{n+1}I. The same H (1.89356973) and C
   (0.45730578) come out of the reconstructed field directly.
2. **Eigenvalue re-derived.** Two-parameter shooting on (s, f″(0)) enforcing both
   f(π)=0 and f′(π)=0 at n = 6.5 returns s* = 1.8666666667 vs exact (2n+1)/(n+1)
   (difference 5×10⁻¹⁵), i.e. p = 2−s* = 0.1333333333; the log–log slope of the
   reconstructed σθθ(r,0) is −0.1333333333.
3. **Independent compatibility check.** The two constructions of ũθ (εθθ quadrature vs
   shear-strain relation) agree pointwise to 1.5×10⁻⁹ (max |ũθ| ≈ 2.2).
4. **External anchors.** The same pipeline at n = 3 (plane strain) reproduces
   Hutchinson's Table 1: I = 5.50732 (paper 5.51), E_n = 1.22716 (paper 1.23); the
   paper's plane-strain I row interpolated (cubic, log–log) to n = 6.5 gives 4.810 vs
   our direct 4.818 (0.2%, table-rounding level).
5. **Grid convergence.** n_θ = 2001→32001: H = 1.89356973 constant to 8 digits, C
   constant to 8 digits, I stable to 7 digits.

Files: `compute.py` (results), `verify3.py` (checks 1–3), `debug_ut.py`,
`verify.py` (superseded intermediates).

## Physical sanity

H grows with n toward the perfectly plastic Prandtl value σθθ(0) → (2+π)σ₀/√3 ≈ 2.97σ₀
(H = 1.267 at n = 3, 1.894 at n = 6.5), as expected since [J/(αε₀σ₀r)]ᵖ → 1; the
high-triaxiality, low-σe character of the plane-strain crack front (σ̃e(0) = 0.55 vs
max 1 at θ ≈ 97°) is reproduced.
