# Topic: crack-tip and notch-tip asymptotic fields

Version 0.2, frozen 2026-10-07 (logs P244, P249; v0.1 at P240). Part of the mechanics library: a **topic**
holds what several papers share; each paper is a **knowledge entry** under `entries/`. See
`pilots/KNOWLEDGE_ENTRY_DESIGN.md` for the format. Nothing here is a recipe: use any method,
and use the checks in `checks/` to test your own solution. Section 6 says where the solution
holds; read it before applying a separable field to a real problem.

Labels used below: **[derived]** shown here or computed with this topic's checks;
**[standard]** textbook mechanics; **[from memory]** a literature result quoted from memory and
not yet checked against its source.

## 1. The question this topic answers

What is the leading stress and strain field near the tip of a sharp crack or notch, how is its
strength measured, and **when does that description hold**? The candidate answer is a
separable field

  σᵢⱼ(r, θ) = A r^(−p) Σᵢⱼ(θ),

where p is the stress exponent, Σᵢⱼ(θ) the angular functions, and A an amplitude that only the
remote problem (geometry and load) can fix.

## 2. Governing equations and assumptions

Notation: the material occupies −β < θ < β; the faces are at θ = ±β (β = π: a crack; β < π: a
notch). α is reserved for the hardening constant below.

| Assumption | Where it enters |
|---|---|
| Small strain | Strain–displacement relations; fails as r → 0 (§6, item 1) |
| Faces traction-free at θ = ±β, and open | Boundary conditions Σθθ = Σrθ = 0 at θ = ±β; admissible only if the faces separate (§6, item 5) |
| Plane strain (εzz = 0) or plane stress (σzz = 0) | The out-of-plane stress in the constitutive law; each is a limit of a 3D body (§6, item 8) |
| A strain energy exists (linear elasticity, or J2 deformation theory) | Path independence of J; the energy argument for the crack exponent |
| Deformation theory = flow theory | Only under proportional loading (§6, item 6) |
| One term dominates | Only in an annulus around the tip (§6, items 1–2) |

Equations [standard]:

- **Equilibrium** with σ = A r^(−p) Σ(θ): (1 − p)Σrr − Σθθ + Σ′rθ = 0 and (2 − p)Σrθ + Σ′θθ = 0.
  An Airy function φ = A r^(2−p) F(θ) satisfies these identically, with
  Σrr = F″ + (2 − p)F, Σθθ = (2 − p)(1 − p)F, Σrθ = −(1 − p)F′.
- **Constitutive laws.** Linear: Hooke's law. Power law (J2 deformation theory):
  εᵢⱼ = (3/2) α ε₀ (σₑ/σ₀)^(n−1) sᵢⱼ/σ₀, with sᵢⱼ the deviator and σₑ = √(3sᵢⱼsᵢⱼ/2).
  Plane strain with the pure power law is incompressible, so σzz = (σrr + σθθ)/2.
  For an elastic–power-law material with n > 1, the power-law term dominates as r → 0 because
  the stress is singular; elastic strains are of lower order.
- **Compatibility.** Strains scale as r^(−np) (power law) or r^(−p) (linear) and must derive
  from a displacement field: εrr = ∂ur/∂r, εθθ = (ur + ∂uθ/∂θ)/r,
  2εrθ = (1/r)∂ur/∂θ + ∂uθ/∂r − uθ/r, with εrθ the **tensor** shear (half the engineering shear).

The result is an **eigenvalue problem**: nontrivial Σ(θ) exist only for special p.

## 3. How the exponent is fixed

- **Admissibility.** The strain energy in any finite neighbourhood of the tip must be finite:
  p(n + 1) < 2 for the power law, p < 1 for linear elasticity. Among admissible eigenvalues,
  the most singular one gives the leading term as r → 0. [standard]
- **Crack (β = π), flanks parallel to x.** J is path-independent (§4). For a separable field
  J ∝ r^(1 − p(n+1)), so a finite, nonzero J requires **p = 1/(n + 1)**: p = 1/2 in the linear
  case, the HRR exponent for the power law. [derived]
- **Notch (β < π).** The flanks are not parallel to x, so J taken along them does not vanish, J is
  not path-independent, and J does not fix p. The exponent comes from the eigenvalue problem
  itself (e.g. shooting on p and one free initial value, with the face conditions as targets).
  For a power-law notch J on a circle grows with r, which the checks confirm. [derived]
- **Linear (Williams) exponents**, λ = 1 − p, smallest real root of [standard; checked here]
  mode I: sin(2λβ) + λ sin 2β = 0;  mode II: sin(2λβ) − λ sin 2β = 0. In mode II, λ = 1 always
  solves the equation and is the rigid rotation, not a stress field; it counts only where it is a
  double root (tan 2β = 2β, the critical angle below), which gives p = 0.
- **When a singularity exists at all** [derived, `checks/test_checks.py::test_singularity_thresholds`]:

  | total wedge angle 2β | 180° | 240° | 257.45° | 270° | 300° | 360° (crack) |
  |---|---|---|---|---|---|---|
  | mode I p | 0 | 0.384 | 0.431 | 0.456 | 0.488 | 0.5 |
  | mode II p | none singular | none singular | 0 (double root) | 0.092 | 0.269 | 0.5 |

  Mode I is singular only at a re-entrant corner (2β > 180°); mode II only beyond 257.45°,
  where λ = 1 becomes a double root of the mode II equation (tan 2β = 2β). A shear singularity
  needs a sharper corner than a tensile one.
- **Limits to check any new solution against:** n = 1 (the power law becomes incompressible
  linear elasticity, and traction-free eigenvalues do not depend on ν); β → π (a crack); for
  the power law with a crack, p = 1/(n+1).

## 4. Amplitude, J and conventions

- **J** = ∫Γ (W n₁ − tᵢ ∂uᵢ/∂x₁) ds, with W the strain energy density (W = n/(n+1) σᵢⱼεᵢⱼ for
  the pure power law, ½σᵢⱼεᵢⱼ for linear elasticity). Path-independent when a strain energy
  exists, the material is homogeneous in x₁, there are no body forces, and the flanks are
  traction-free and parallel to x₁. [standard]
- **Compute ∂uᵢ/∂x₁ in Cartesian components**, or include the derivatives of the rotating polar
  basis vectors. Dropping them gives a plausible wrong number: I ≈ 5.65 instead of 3.86 for
  n = 3 plane stress, an error made and later traced in this library.
- **Agreement of J on a few contours is a necessary check, not a proof.** For a separable field
  J is r-independent as soon as p = 1/(n+1), so contour agreement tests only the exponent.
- **Normalization.** Σ(θ) is defined up to a factor. Common conventions: max σ̃ₑ = 1
  (Hutchinson), F(0) = 1, Σθθ(0) = 1. The crack field is then
  σᵢⱼ = σ₀ [J/(α ε₀ σ₀ Iₙ r)]^(1/(n+1)) σ̃ᵢⱼ(θ), and Iₙ depends on the convention. **State the
  convention whenever you quote Iₙ or an amplitude**; I, K and the field must share it.
- **Plane strain J and K (linear):** J = (1 − ν²)K²/E. The factor (1 − ν²) belongs to the linear
  elastic relation only; it does not enter J = α ε₀ σ₀ Iₙ K^(n+1) for the power-law field.
- The checks set all material constants to 1; a reported J or I is compared in those units.

## 5. Leaving the scope

| If this changes | What happens | What still holds |
|---|---|---|
| Notch instead of crack | p is an eigenvalue; J does not fix it; J on an arc depends on r | Equilibrium, compatibility, face conditions, admissibility; n = 1 and β → π limits |
| Mode II or mixed mode | Antisymmetric (or mixed) angular functions; different eigenvalues; for nonlinear laws no superposition (§6, item 7) | The same equations with other symmetry conditions |
| Flow theory, unloading or non-proportional loading | No strain energy; deformation-theory results become approximations (§6, item 6) | Equilibrium and compatibility |
| Finite strain (e.g. rubber) | Different asymptotic structure; small-strain fields do not apply (Knowles & Sternberg 1973) | Energy arguments in the finite-strain form |
| Crack between dissimilar materials | Linear interface cracks have oscillatory singularities (Williams 1959; England 1965) and predict face interpenetration very near the tip (§6, item 5) | Energy release rate as the robust measure |
| Growing crack | Elastic unloading behind the tip; the stationary field does not describe steady growth | Equilibrium; energy balance |
| Distance from the tip comparable to geometry | Non-singular terms (T-stress, higher-order terms) matter (Larsson & Carlsson 1973) | The leading term as the inner limit |
| Very large n (perfect plasticity) | The smooth separable structure changes (sectors, possible discontinuities) | Equilibrium and yield condition |

## 6. Where the separable solution holds (applicability conditions)

The separable field is an **intermediate asymptotic**: it can hold only in an annulus
r_in ≪ r ≪ r_out around the tip, and only under conditions it does not itself contain. Each
row is physics the solution leaves out; when a condition fails, something must be added (another
model, another term, or another criterion). Before applying a separable field, say which rows
hold for your problem and what you add where one does not. The checks (§7) test whether a field
solves the model; this section is about whether the model fits.

| # | Condition | Physics the solution leaves out | When it fails | Basis |
|---|---|---|---|---|
| 1 | **Inner limit**: r large compared with the zone where small strain, continuum or the material law fail | The singular strains violate small strain as r → 0, so every singular solution is self-inconsistent at the tip | Crack-tip blunting and large strains in metals over a few crack-opening displacements δ_t, with δ_t of order J/σ₀; fracture process zone; microstructure | [standard]; the size "a few δ_t" and δ_t ≈ d_n J/σ₀ with d_n ≈ 0.3–0.8 [from memory: McMeeking 1977; Shih 1981] |
| 2 | **Outer limit**: one term dominates the next ones where it matters | The leading term ignores the remote geometry; the next terms (T-stress for linear fields, Q or A₂ for power-law fields) measure constraint | One-parameter characterization (K or J alone) fails; use two parameters (J–T, J–Q) or the full solution. For J dominance, ligament and crack length must be large compared with J/σ₀ | [standard] (Larsson & Carlsson 1973 for T); size factors ≈ 25 J/σ₀ in bending, ≈ 200 J/σ₀ in tension [from memory: McMeeking & Parks 1979] |
| 3 | **A singularity exists** | — | Below the thresholds in §3 (mode I: 2β ≤ 180°; mode II: 2β ≤ 257.45°) stresses stay bounded: there is no tip field to characterize, and failure is a strength (stress) question | [derived] |
| 4 | **Eigenvalues are distinct** | The pure power form assumes no two exponents coincide | Where exponents coincide (e.g. λ = 1 is a double root in mode II at 2β = 257.45°), or a forced higher-order term meets an eigenvalue, the separable form is incomplete and r^λ ln r terms appear; near such angles the "leading" term is ill-conditioned | double root [derived]; its connection to the Sternberg–Koiter wedge paradox at the same angle [from memory] |
| 5 | **The faces open** | Traction-free faces assume separation | If the field would make crack faces interpenetrate (closing load, U_θ(π) > 0), contact (and friction) is the missing physics. Linear interface cracks predict interpenetration in a very small zone at the tip, so the oscillatory field is itself an approximation there | opening condition [derived, check "crack faces open"]; interface contact zones [from memory: Comninou 1977] |
| 6 | **Loading is proportional** | Deformation theory replaces flow theory | Inside the separable zone, under monotonic J, the stress at a fixed material point is σ ∝ J^(1/(n+1)) × a fixed tensor, so it is exactly proportional there. Proportionality fails in the elastic–plastic transition, under crack growth (elastic unloading behind the tip, non-proportional ahead) and under cyclic loading | proportionality in the separable zone [derived from the separable form]; the rest [standard] |
| 7 | **Superposition** (linear fields only) | The power-law problem is nonlinear | Mode I and mode II power-law fields cannot be added; mixed mode needs its own solution and a mixity parameter, and the near-tip mixity need not equal the far-field one | [standard]; mixity parameterization [from memory: Shih 1974] |
| 8 | **2D idealization** | A 3D body near its crack front | Plane strain holds near the front where r is small compared with the thickness; plane stress only where r is large compared with it. A plane-stress tip field therefore needs its zone to exceed the thickness | [standard] |
| 9 | **From field to fracture** | The amplitude is a field parameter, not a fracture criterion | For a crack, K or J is a fracture parameter through the energy argument (under items 1–2). A notch amplitude has units of stress × length^p and is not by itself one: crack onset at a notch needs both a stress and an energy condition (Leguillon 2002) | crack part [standard]; the coupled notch criterion [from memory: Leguillon 2002, in the corpus as P040, not yet read here] |

## 7. Checks any solution must satisfy (`checks/crack_tip_checks.py`)

`check_field(theta, S_rr, S_tt, S_rt, p, n, beta, plane, law, nu, mode, leading, tol, tol_eq, reported_J, h, h2)`
takes the angular stress functions of a candidate solution (from any code) and reports, each
with its conditions and tolerance:

1. traction-free faces; 2. the declared symmetry (mode I or II);
3. equilibrium in angular form (tight tolerance, so it also tests the exponent), and
4. in Cartesian components by finite differences;
5. compatibility of the strains from the declared law (Saint-Venant, Cartesian);
6. finite energy near the tip; for a leading term with p ≤ 0, a report that the field is not
   singular (§6, item 3);
7. the crack exponent fixed by J (leading term of a crack) or, for a linear notch, the Williams
   exponent (modes I and II); `leading=False` skips these for higher-order or non-singular terms.
   A crack declared with a non-singular leading term (p ≤ 0) fails, because a crack's leading
   term is non-singular only if the singular amplitude vanishes: check a T-stress with `leading=False`;
8. **mode I only**: crack faces do not interpenetrate (U_θ(π) ≤ tol × max|U|; §6, item 5; for a
   notch the face displacement is reported only); compatibility through displacements
   integrated from the strains; J computed from the field in Cartesian components (report),
   compared with a reported J or I if given (crack only);
9. the normalization of the field (report) and a resolution warning on coarse grids.

The checks import no reference solver. They were validated by tests kept in the library
(`checks/test_checks.py`): passing on correct fields (Williams crack, notches in modes I and II,
a 90° wedge, a T-stress field, HRR, a power-law notch) and failing on planted errors, each for
the intended reason; `checks/mutant_detection.txt` lists which checks catch which error.
Tolerances assume angular functions resolved to about 10⁻⁶; if a check fails marginally on a
coarse grid, refine before blaming the model.

## 8. Entries in this topic

| Entry | Status |
|---|---|
| Hutchinson 1968, power-law crack-tip field (P015) | `entries/hutchinson_1968/` (prototype) |
| Rice & Rosengren 1968 (P016); Williams 1957 (P006); Rice 1968, the J integral (P009) | planned |
| Knowles & Sternberg 1973, finite strain (P086); Larsson & Carlsson 1973, non-singular terms (P013) | planned |
| Williams 1959 (P052), England 1965 (P053): interface cracks; Leguillon 2002 (P040): notch onset | planned |

## Change log

- **v0.2 frozen (P249)** after Grok's re-check (`corpus-runs/library/P015-recheck-grok/RECHECK.md`,
  verdict "freeze after the listed small fixes"; all v0.1 high-severity findings resolved; §6 judged
  physically right and honestly labelled). Fixes: the mode II rigid-rotation root λ = 1 is divided
  out of the equation, so the double root at 257.45° gives p = 0 and nearby roots are no longer
  missed (it crashed there); a crack with a non-singular leading term keeps failing, with the
  physical reason stated; the face check tests for interpenetration with a tolerance (a T-stress
  passed or failed on round-off); an all-zero field is rejected. Tests: 27.
- **v0.2 (P244)** after Grok's review (`corpus-runs/library/P015-review-grok/REVIEW.md`, verdict
  "usable after the listed fixes") and Teng's direction to check the physics behind each
  failure, not only the solution. The checks: residuals normalized by the field scale (correct
  T-stress and mode II fields were rejected); the Williams exponent in modes I and II, returning
  non-singular roots (a 90° wedge crashed); no crash at 1 + q = 0; angular-equilibrium tolerance
  10⁻⁵ (an exponent error of 10⁻⁴ on a power-law notch passed); `leading` switch; wedge
  half-angle renamed β; the crack-opening check (a closing crack passed); resolution warning.
  The topic: §6 added (each review finding traced to the physics it reveals); existence
  thresholds computed (§3); tensor-shear convention stated; mode II equation added.
- **v0.1 (P240)** first prototype.
