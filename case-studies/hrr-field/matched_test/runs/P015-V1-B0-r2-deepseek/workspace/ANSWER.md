# Answer: slip and strain contours near a stationary slit (mode-III power law)

## Requested numbers

All in mm, N, MPa. The leading field is used for all `r <= 0.30 mm`.

**1. Relative face slip Δw(r) = w(r,+π) − w(r,−π)**

| r (mm) | Δw (mm) |
|---|---|
| 0.0200 | **1.31813 × 10⁻³** |
| 0.0800 | **1.79369 × 10⁻³** |

**2. Signed stresses at r = 0.0200 mm, θ = 60°** (e_r = (cosθ,sinθ), e_θ = (−sinθ,cosθ))

| component | value (MPa) |
|---|---|
| τ_rz = τ·e_r | **+38.5866** |
| τ_θz = τ·e_θ | **+138.011** |

**3. Radii where |γ| = 0.0100**

| ray | r (mm) |
|---|---|
| θ = 0° | **0.046894** |
| θ = 90° | **0.025066** |
| θ = 180° (upper face, limit from inside) | **0.013398** |

## Governing model and how the loading fixes its magnitude

The problem is **antiplane shear (mode III)**: the only displacement is `w(x,y)`, the
engineering shear is `γ = ∇w`, and the constitutive law couples only the shear vectors,

    γ = (γ0/τ0)(τ/τ0)^(n-1) τ   ⇔   τ = τ0 (|∇w|/γ0)^(1/n) (∇w/|∇w|).

Substituting into equilibrium (no body forces, ∂/∂z = 0) gives the scalar **p-Laplacian**

    ∇·( |∇w|^(1/n − 1) ∇w ) = 0,

with odd symmetry `w(x,−y) = −w(x,y)` (mode III) and traction-free faces
`τ_θz(r,±π) = 0`, i.e. `∂w/∂θ = 0` at θ = ±π.

The leading field is separable, `w = A r^s f(θ)`, `f` odd. Equilibrium fixes the
angular ODE for `f`; the exponent is fixed by **J**. With `γ ∼ r^(s−1)` and
`τ ∼ γ^(1/n)` one has `W ∼ τγ ∼ r^((s−1)(1+1/n))`, so

    J ∼ ∫ W r dθ ∼ r^(1 + (s−1)(1+1/n)).

A finite, non-zero J therefore requires `s = 1/(n+1) = 0.222222`, and the stresses
carry the HRR exponent `τ ∼ r^(−1/(n+1))`. Solving the angular ODE

    (g² + q f′²) f″ = −q s² f f′² − (1−p) s f g²,
    g² = s²f² + f′²,  p = s = 1/(n+1),  q = 1/n − 1,

on 0 ≤ θ ≤ π with `f(0)=0`, `f′(0)=1`, `f′(π)=0` (odd extension to θ<0) gives
`f(π)=1.698433`. The **amplitude** A is then fixed by the given flux:

    J = τ0 γ0^(−1/n) A^((n+1)/n) · I,   I = ∫_{−π}^{π}[ (n/(n+1)) g^(1+1/n) cosθ
        − g^(1/n−1) s f (s f cosθ − f′ sinθ) ] dθ,

with `I = 2.019595`, so `A = 9.25609 × 10⁻⁴ mm^(1−s)`. The requested slip is
`Δw(r) = 2 A r^s f(π)`; the stresses are

    τ_rz = τ0 γ0^(−1/n) A^(1/n) r^((s−1)/n) g^(1/n−1) (s f),
    τ_θz = τ0 γ0^(−1/n) A^(1/n) r^((s−1)/n) g^(1/n−1) f′,

and `|γ| = A r^(s−1) g(θ)`, so `r = (A g/0.01)^(1/(1−s))`.

## Assumptions used / changed

- **Model choice.** The specified law is a pure power law with no additive elastic term,
  so there is no elastic–plastic transition: the p-Laplacian above is the *exact* governing
  equation, and the single power-law exponent `1/(n+1)` is the leading one. No superposition
  or elastic term was added.
- The leading term was extended to `r = 0.30 mm` as the task instructs ("geometry outside
  r = 0.30 mm introduces negligible corrections"); small strain and continuum are taken as
  valid in the requested range, which is an assumption of the task, not of the solution.
- The amplitude sign was fixed by the stated loading (upper face moves +z ⇒ A > 0, Δw > 0);
  the field has no additive rigid translation (task statement).

## Verification (quantitative)

1. **Independent J check (energy flux / path independence).** The J contour integral was
   evaluated directly in Cartesian components from the field at four radii:
   J = 0.120000000 N/mm at R = 0.02, 0.05, 0.10 and 0.30 mm — exactly the prescribed
   0.120 N/mm and independent of R, as required for a crack with `s = 1/(n+1)`.
2. **Known limiting solution (n = 1).** For `n = 1` the same code gives `s = 1/2`,
   `I = π`, and reproduces the classical linear mode-III crack: with μ = τ0/γ0 = 40000 MPa,
   K_III = 2Aμ√(π/2), G = K_III²/(2μ) = 0.120000 N/mm. The exponent and normalization are
   therefore anchored to the exact linear result.
3. **Boundary residual.** The free-face condition is satisfied to machine precision,
   `f′(π) = 1.4 × 10⁻¹⁶`. A Cartesian finite-difference divergence of τ over interior angles
   (θ = 10°…175°, r = 0.05…0.20 mm) is at the few × 10⁻⁶ level relative to `|τ|/r`,
   consistent with discretization error of the singular field.

## What is missing / uncertain

- The linear (in-plane) `crack_tip_checks.py` provided in `./knowledge` does not cover
  antiplane mode III, so its `check_field` was not applied; the mode III verification above
  (J path-independence and the `n = 1` limit) replaces it.
- The angular shape `f(θ)` has no closed form used here; it is a high-accuracy numerical
  solution of the ODE, stable to ~10⁻¹² across independent integrators (DOP853, Radau) and
  grids, far inside the 2% tolerance.
- At `r` of the requested points the strains are `~10⁻²` and stresses tens of MPa, so
  small strain is self-consistent there; no claim is made nearer the tip.

## Does this determine the load at which the slit advances?

**No.** J = 0.120 N/mm is only the amplitude (state) of the current near-tip field; onset of
advance requires comparing it with a separately measured critical value (toughness / cohesive
criterion), which these results do not contain.
