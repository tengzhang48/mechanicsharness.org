# Topic: crack-tip and notch-tip asymptotic fields

Prototype, 2026-10-06 (log P240). Part of the mechanics library: a **topic** holds what several
papers share; each paper is a **knowledge entry** under `entries/`. See
`pilots/KNOWLEDGE_ENTRY_DESIGN.md` for the format. Nothing here is a recipe: use any method,
and use the checks in `checks/` to test your own solution.

## 1. The question this topic answers

What is the leading stress and strain field near the tip of a sharp crack or notch, how is its
strength measured, and when does that description hold? The answer is a separable field

  σᵢⱼ(r, θ) = A r^(−p) Σᵢⱼ(θ),

where p is the stress exponent, Σᵢⱼ(θ) the angular functions, and A an amplitude that only the
remote problem (geometry and load) can fix.

## 2. Governing equations and assumptions

Assumptions, and where each enters (2D, quasi-static, no body force, stationary tip):

| Assumption | Where it enters |
|---|---|
| Small strain | Strain–displacement relations below; finite strain changes the problem (§5) |
| Faces traction-free at θ = ±α (α = π: crack) | Boundary conditions Σθθ = Σrθ = 0 at θ = ±α |
| Plane strain (εzz = 0) or plane stress (σzz = 0) | The out-of-plane stress in the constitutive law |
| A strain energy exists (linear elasticity, or J2 deformation theory) | Path independence of J; the energy argument for the crack exponent |
| Deformation theory ≈ flow theory | Only under proportional loading without unloading (§5) |
| The leading term dominates | Only for r small compared with the geometry; the size of that region is not given by the asymptotics |

Equations:

- **Equilibrium** with σ = A r^(−p) Σ(θ): (1 − p)Σrr − Σθθ + Σ′rθ = 0 and (2 − p)Σrθ + Σ′θθ = 0.
  An Airy function φ = A r^(2−p) F(θ) satisfies these identically, with
  Σrr = F″ + (2 − p)F, Σθθ = (2 − p)(1 − p)F, Σrθ = −(1 − p)F′.
- **Constitutive laws.** Linear: Hooke's law. Power law (J2 deformation theory):
  εᵢⱼ = (3/2) α ε₀ (σₑ/σ₀)^(n−1) sᵢⱼ/σ₀, with sᵢⱼ the deviator and σₑ = √(3sᵢⱼsᵢⱼ/2).
  Plane strain with the pure power law is incompressible, so σzz = (σrr + σθθ)/2.
  For an elastic–power-law material (Ramberg–Osgood type) with n > 1, the power-law term
  dominates as r → 0 because the stress is singular; elastic strains are of lower order.
- **Compatibility.** Strains scale as r^(−np) (power law) or r^(−p) (linear); the strains must
  derive from a displacement field. In polar components:
  εrr = ∂ur/∂r, εθθ = (ur + ∂uθ/∂θ)/r, 2εrθ = (1/r)∂ur/∂θ + ∂uθ/∂r − uθ/r.

The result is an **eigenvalue problem**: nontrivial Σ(θ) exist only for special p.

## 3. How the exponent is fixed

- **Admissibility.** The strain energy in any finite neighbourhood of the tip must be finite:
  p(n + 1) < 2 for the power law, p < 1 for linear elasticity. Among admissible eigenvalues,
  the most singular one gives the leading term as r → 0.
- **Crack (α = π), flanks parallel to x.** J is path-independent (see §4). For a separable
  field J ∝ r^(1 − p(n+1)), so a finite, nonzero J requires **p = 1/(n + 1)**: p = 1/2 in the
  linear case, the HRR exponent for the power law.
- **Notch (α < π).** The flanks are not parallel to x, so J taken along the flanks does not
  vanish, J is not path-independent, and J does not fix p. The exponent must be found from
  the eigenvalue problem itself (for example by shooting on p and one free initial value, with
  the face conditions as targets). Linear case (Williams), mode I: λ = 1 − p is the smallest root
  in (0, 1] of sin(2λα) + λ sin(2α) = 0.
- **Limits to check any new solution against:** n = 1 (the power law becomes incompressible
  linear elasticity, and traction-free eigenvalues do not depend on ν); α → π (a crack);
  for the power law with a crack, p = 1/(n+1).

## 4. Amplitude, J and conventions

- **J** = ∫Γ (W n₁ − tᵢ ∂uᵢ/∂x₁) ds, with W the strain energy density (W = n/(n+1) σᵢⱼεᵢⱼ for
  the pure power law, ½σᵢⱼεᵢⱼ for linear elasticity).
  It is path-independent when a strain energy exists, the material is homogeneous in x₁, there
  are no body forces, and the flanks are traction-free and parallel to x₁.
- **Compute ∂uᵢ/∂x₁ in Cartesian components**, or include the derivatives of the rotating polar
  basis vectors. Dropping them gives a plausible wrong number: I ≈ 5.65 instead of 3.86 for
  n = 3 plane stress, an error made and later traced in this library.
- **Agreement of J on a few contours is a necessary check, not a proof.** For a separable field,
  J is r-independent as soon as p = 1/(n+1), so contour agreement tests only the exponent.
- **Normalization.** Σ(θ) is defined up to a factor. Common conventions: max σ̃ₑ = 1
  (Hutchinson), F(0) = 1, Σθθ(0) = 1. The crack field is then
  σᵢⱼ = σ₀ [J/(α ε₀ σ₀ Iₙ r)]^(1/(n+1)) σ̃ᵢⱼ(θ), and Iₙ depends on the convention. **State
  the convention whenever you quote Iₙ or an amplitude**; I, K and the field must share it.
- **Plane strain J and K (linear):** J = (1 − ν²)K²/E. The factor (1 − ν²) belongs to the linear
  elastic relation only; it does not enter J = α ε₀ σ₀ Iₙ K^(n+1) for the power-law field.

## 5. Leaving the scope

| If this changes | What happens | What still holds |
|---|---|---|
| Notch instead of crack | p is an eigenvalue; J does not fix it; J on an arc depends on r | Equilibrium, compatibility, face conditions, admissibility; n = 1 and α → π limits |
| Mode II or mixed mode | Antisymmetric (or mixed) angular functions; different eigenvalues; mixity enters | The same equations with other symmetry conditions |
| Flow theory, unloading or non-proportional loading | No strain energy; deformation-theory results become approximations, valid only if loading is nearly proportional (this is usually assumed, not shown) | Equilibrium and compatibility |
| Finite strain (e.g. rubber) | Different asymptotic structure; the small-strain fields do not apply (Knowles & Sternberg 1973) | J-type energy arguments, in the appropriate (finite-strain) form |
| Crack between dissimilar materials | Linear interface cracks have oscillatory singularities (Williams 1959; England 1965) | Energy release rate as the robust measure |
| Growing crack | Elastic unloading behind the tip; the stationary field does not describe steady growth | Equilibrium; energy balance |
| Distance from the tip comparable to geometry | Non-singular terms (T-stress, higher-order terms) matter (Larsson & Carlsson 1973) | The leading term as the inner limit |
| Very large n (perfect plasticity) | The smooth separable structure changes (sectors, possible discontinuities) | Equilibrium and yield condition |
| Scales comparable to microstructure | The continuum field itself fails | — |

Passing the checks in §6 shows that a field is a correct solution **of the stated model**. Whether
the model fits the problem is a separate question, answered with this table and with physical
arguments, and it should be justified explicitly whenever an assumption is changed.

## 6. Checks any solution must satisfy (`checks/crack_tip_checks.py`)

`check_field(theta, S_rr, S_tt, S_rt, p, n, alpha, plane, law, nu, mode, reported_J)` takes the
angular stress functions of a candidate solution (from any code) and reports, each with its
conditions and tolerance:

1. traction-free faces; 2. the declared symmetry (mode I or II);
3. equilibrium in angular form, and 4. in Cartesian components by finite differences;
5. compatibility of the strains from the declared law (Saint-Venant, Cartesian);
6. finite energy near the tip; 7. the crack exponent fixed by J (crack only) or, for a linear
   notch, the Williams exponent;
8. compatibility through displacements integrated from the strains (mode I);
9. J computed from the field in Cartesian components (report), compared with a reported J or I
   if one is given; and the normalization of the field (report).

The checks import no reference solver. `checks/test_checks.py` shows them passing on correct
fields (Williams crack and notch, HRR from the Hutchinson entry's reference code, a power-law
notch) and failing on planted errors, each for the intended reason; `checks/mutant_detection.txt`
lists which checks catch which error. Tolerances assume angular functions resolved to about
10⁻⁶; if a check fails marginally on a coarse grid, refine before blaming the model.

## 7. Entries in this topic

| Entry | Status |
|---|---|
| Hutchinson 1968, power-law crack-tip field (P015) | `entries/hutchinson_1968/` (prototype) |
| Rice & Rosengren 1968 (P016); Williams 1957 (P006); Rice 1968, the J integral (P009) | planned |
| Knowles & Sternberg 1973, finite strain (P086); Larsson & Carlsson 1973, non-singular terms (P013) | planned |
| Williams 1959 (P052), England 1965 (P053): interface cracks | planned |
