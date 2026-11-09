# Answer: slip and strain contours near a stationary slit (Mode III, power law)

## 1. Requested numbers

Loading: J = 0.120 N/mm; material τ0 = 80 MPa, γ0 = 0.002, n = 3.5 (pure power law).

**(1) Relative face slip** Δw(r) = w(r,+π) − w(r,−π):

| r (mm) | Δw (mm) |
|---|---|
| 0.0200 | **1.3181 × 10⁻³** |
| 0.0800 | **1.7937 × 10⁻³** |

**(2) Signed stresses at r = 0.0200 mm, θ = 60°** (e_r = (cosθ, sinθ), e_θ = (−sinθ, cosθ)):

- τ_rz = **+38.59 MPa**
- τ_θz = **+138.01 MPa**

(both positive; |τ| = 143.30 MPa, |γ| = 1.539 × 10⁻² at this point).

**(3) Radii where |γ| = 0.0100** intersects the rays:

| ray | r (mm) |
|---|---|
| θ = 0° | **4.6894 × 10⁻²** |
| θ = 90° | **2.5066 × 10⁻²** |
| θ = 180° (upper face, limit from within) | **1.3398 × 10⁻²** |

All seven numbers come from one leading-order field (below); numerical errors are
≲ 10⁻⁸ relative, far inside the 2% tolerance. All requested locations lie inside
r = 0.30 mm, where the problem statement says the leading field suffices.

## 2. Governing model and why

This is anti-plane shear (Mode III): the only displacement is w(x,y), the engineering
shear vector is γ = ∇w, equilibrium is div τ = 0, and the constitutive law makes γ
parallel to τ with |γ| = γ0(τ/τ0)ⁿ. Because the loading is quasistatic, monotonic and
proportional (fixed pattern), the specified total (deformation-theory) power law is also
the correct incremental description at every fixed material point, so a nonlinear-elastic
(strain-energy) formulation is exact and Rice's J integral is path independent.

I use the standard **HRR-type separable near-tip field** for anti-plane shear,

  w(r,θ) = A r^q f(θ),  γ_r = A q r^(q−1) f,  γ_θ = A r^(q−1) f',
  τ = τ0(A/γ0)^{1/n} r^s G^{1/n},  τ_rz = τ·(qf/G),  τ_θz = τ·(f'/G),

with G(θ) = √(q²f² + f'²) and s = −1/(n+1). The exponent is not free: on a circle of
radius r, W ~ τγ ~ r^{(q−1)(n+1)/n} and J ~ r·W, so a finite, non-zero, path-independent J
forces **q = 1/(n+1) = 1/4.5 = 0.222222** (stress exponent s = −1/(n+1)); this is the crack
case of the general argument in `knowledge/TOPIC.md` §3.

Equilibrium, (1/r)(rτ_rz)_,r + (1/r)τ_θz,θ = 0, reduces to one angular ODE (p = f'):

  p' = −q² f · (n² − (n−1)η) / (n − (n−1)η),  η = p²/(q²f² + p²),

with f(0) = 0 (antisymmetry w(x,−y) = −w(x,y)) and f'(π) = 0 (traction-free face,
τ_θz(r,π) = 0). The ODE is homogeneous in f, so I integrate it with the normalization
f'(0) = 1 (DOP853, rtol 10⁻¹²). With q = 1/(n+1) the solution satisfies f'(π) = 0 to
3 × 10⁻¹⁴ without any shooting — a strong internal consistency check of the field.

## 3. How the loading fixes the magnitude

J taken on a circle of radius r (CCW from θ = −π to +π, N radial, ds = r dθ) is

  J = ∫ (W cosθ − τ_rz w_,x) r dθ = τ0γ0 (A/γ0)^{(n+1)/n} I,

with w_,x = γ_r cosθ − γ_θ sinθ, W = (n/(n+1))τ0γ0(γ/γ0)^{(n+1)/n} (proportional loading),
and the dimensionless angular integral (in the f'(0) = 1 normalization)

  I = ∫_{-π}^{π} [(n/(n+1))G^{(n+1)/n}cosθ − G^{1/n−1}qf(qf cosθ − f' sinθ)] dθ = 2.019595.

Equating to the applied energy flux J = 0.120 N/mm gives the single free constant:

  A = γ0 [J/(τ0γ0 I)]^{n/(n+1)} = 9.25609 × 10⁻⁴ mm^{n/(n+1)}.

Everything requested follows: Δw(r) = 2A f(π) r^q with f(π) = 1.698433; the stresses from
the formulas above; and the strain contour from |γ| = A G(θ) r^{−n/(n+1)} = 0.01, i.e.
r(θ) = (A G(θ)/0.01)^{(n+1)/n}, with G(0) = 1, G(90°) = 0.614353, G(180°) = qf(π) = 0.377430.

## 4. Assumptions

- No assumption of the problem statement was changed. In particular I kept the pure power
  law with no elastic term (as specified; the singular term dominates as r → 0 anyway) and
  used the leading separable term everywhere, with no rigid translation, per the statement.
- Deformation theory is used as the constitutive model; this is exact here because the
  separable field under monotonically increasing J loads every material point
  proportionally (stress direction fixed, magnitude ∝ J^{1/(n+1)}).
- Units mm, N, MPa throughout; angles in radians.

## 5. Verification (quantitative)

- **Known limiting solution (n = 1):** the identical code at n = 1 reproduces the exact
  linear-elastic Mode III crack field: max|f − 2 sin(θ/2)| = 2 × 10⁻¹²; the angular
  integral I = π to 1.5 × 10⁻¹³; and the amplitude A = √(γ0J/(πτ0)) matching
  w = (K_III/μ)√(2r/π) sin(θ/2) with J = K_III²/2μ to 7 × 10⁻¹⁴. This validates the ODE,
  the J integral, and the amplitude formula end-to-end.
- **Energy-flux (J) check:** J recomputed from the full n = 3.5 field on circles of radius
  0.005 / 0.020 / 0.100 / 0.300 mm equals 0.120 N/mm to machine precision (max deviation
  1 × 10⁻¹⁶): path independence holds and the amplitude matches the applied flux.
- **Boundary residual:** f'(π) = −3 × 10⁻¹⁴ at ODE level; evaluating the Cartesian stress
  just inside the upper face gives |t_z|/|τ| ≤ 1 × 10⁻⁶ (limited by the finite approach
  distance y = 10⁻⁹ mm), i.e. the faces are traction free.
- **Equilibrium:** central finite differences of the Cartesian stress field at 300 random
  points (10⁻⁴ < r < 0.29 mm) give max |div τ|/(|τ|/r) = 7 × 10⁻¹⁰.
- **Constitutive residual:** |γ| from the field vs γ0(τ/τ0)ⁿ from the stress agrees to
  2 × 10⁻¹⁵ (by construction of τ(γ)).
- **Arithmetic spot checks:** Δw(0.08)/Δw(0.02) = 4^q = 1.36079, and the contour radii and
  stress magnitudes were recomputed by hand from A, f(π), G(θ) — all match the script.

## 6. Missing or uncertain

- The angular functions are numerical (ODE + quadrature); their error (~10⁻¹²) is
  negligible against the 2% tolerance, so the numbers above are not the weak point.
- Real-crack physics outside the model: higher-order terms (ignored per the statement,
  valid for the requested r ≤ 0.3 mm), and crack-tip blunting / finite-strain and
  process-zone effects at distances of order the crack-tip slip itself, which the
  small-strain stationary field cannot describe.
- **Crack advance:** No — these results alone do not determine the load at which the slit
  starts to advance, because the stationary analysis only characterizes the near-tip field
  for a given applied J; initiation requires an additional material fracture criterion
  (e.g. a critical J_IIIc, critical face slip, or critical strain over a microstructural
  distance) that is not specified anywhere in the problem.

## 7. Files

- `solve.py` — solver and all checks (run: `python3 solve.py`); output in `solve_output.txt`.

