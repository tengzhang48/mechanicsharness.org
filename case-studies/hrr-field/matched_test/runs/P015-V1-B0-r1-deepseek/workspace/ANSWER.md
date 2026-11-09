# Answer: slip and strain contours near a stationary slit (antiplane power-law crack)

## Requested numbers

**1. Relative face slip** Δw(r) = w(r,+π) − w(r,−π)

| r (mm) | Δw (mm) |
|---|---|
| 0.0200 | **1.31813 × 10⁻³** |
| 0.0800 | **1.79369 × 10⁻³** |

**2. Signed shear stresses at r = 0.0200 mm, θ = 60°**

- τ_rz = **+38.587 MPa**
- τ_θz = **+138.011 MPa**

(both positive; the positive z direction is the direction the upper face moves)

**3. Radii where |γ| = 0.0100**

| ray | r (mm) |
|---|---|
| θ = 0° | **4.68944 × 10⁻²** |
| θ = 90° | **2.50661 × 10⁻²** |
| θ = 180° (upper face, limit from material) | **1.33984 × 10⁻²** |

## Model and how the load fixes its magnitude

The only displacement is out of plane, u_z = w(x,y), so this is **antiplane shear (mode III)**.
With γ = (w_,x, w_,y) = ∇w and τ = (σ_xz, σ_yz), equilibrium is ∇·τ = 0 and the given
engineering-shear law inverts to

    τ = τ0 (|γ|/γ0)^(1/n) γ/|γ|,        τ0 = 80 MPa, γ0 = 0.002, n = 3.5.

Seek the leading separable tip field w = A r^s f(θ). Then |γ| = A r^(s−1) g(θ) with
g = √((s f)² + f'²), and the stress exponent is p = (1−s)/n. Because the crack flanks are
parallel to x, J is path-independent; for a separable field J ∝ r^(1−p(n+1)). A finite,
nonzero J forces p = 1/(n+1), hence **s = p = 1/(n+1) = 2/9**. The angular ODE (equilibrium,
(1−p)T_r + T_θ' = 0 with T_r = g^(1/n−1) s f, T_θ = g^(1/n−1) f') closes with f(0) = 0
(antisymmetry w(x,−y) = −w(x,y)) and f'(π) = 0 (traction-free face). It gives

    f(π) = 1.6984332340,  f(π/2) = 1.3032395152,  f'(π/2) = 0.5418084176.

The loading fixes only the amplitude A. On a circle J = A^(1+1/n) I_J with
I_J = 953.8687, so A = (J/I_J)^(n/(n+1)) = **9.25609 × 10⁻⁴** (mm units) at J = 0.120 N/mm.

## Assumption changed and why

The provided topic/HRR entry is the in-plane J2 (mode I) deformation-theory problem. The
problem statement forbids that: the only displacement is u_z, so I used the **antiplane (mode III)**
version with the law written directly in engineering shear. I kept the law exactly as given
(pure power law, no elastic term), so the shear strain is γ = (w_,x, w_,y) = ∇w.

## Verification

1. **Energy flux / path independence** (strongest check): integrating
   J = ∫_Γ (W N_x − t_z w_,x) ds numerically from the solved field gives
   J = 0.1200000000 N/mm on circles r = 0.02, 0.08, 0.20, 0.30 mm — independent of r and equal
   to the prescribed value (A was chosen precisely this way).
2. **Known limiting solution**: with n = 1 the ODE and eigenvalue reduce to Laplace, and the
   solver returns w = 2A r^(1/2) sin(θ/2) exactly (f(π) = 2.0000000000).
3. **Boundary/constitutive residuals**: max Cartesian equilibrium residual / term magnitude
   ≈ 3 × 10⁻¹⁰; constitutive-law residual ≈ 3 × 10⁻¹⁶; face traction τ_θz(π) ≈ 6 × 10⁻¹³ MPa;
   |f'(π)| ≈ 3 × 10⁻¹⁵.
4. **Independent integrations**: DOP853 from θ = 0 forward and backward to ±π agree to 10⁻¹⁴, and
   f(−θ) = −f(θ).

## Does this determine the load at which the slit starts to advance?

**No.** These results give the stress/strain field and its amplitude J at the current load, but
not a fracture criterion: deciding onset needs a material toughness (a critical J, or a critical
strain/process-zone criterion), which is not supplied, and the singular field is anyway
inconsistent in the inner process zone where small strain fails.

## What is missing or uncertain

- The 2% tolerance is met with large margin; the main numerical uncertainty is quadrature /
  ODE tolerance, at the 10⁻¹⁰ level.
- The result assumes the leading (separable) field dominates for r ≤ 0.30 mm, as the problem
  states; no non-singular (T-stress-like) corrections were added.
- Units: lengths in mm, stresses in MPa; A carries units of mm^(1−s) and is reported in that
  system.
- Scripts: `final.py` (results), `verify2.py` (residuals and J), `solve2.py` (first pass).
