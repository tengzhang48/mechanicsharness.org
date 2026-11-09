# Answer

## Requested numbers

| # | Quantity | Value |
|---|---|---|
| 1 | Relative face slip Δw at r = 0.0200 mm | **0.0013181 mm** |
| 1 | Relative face slip Δw at r = 0.0800 mm | **0.0017937 mm** |
| 2 | τ_rz at r = 0.0200 mm, θ = 60° | **+38.5866 MPa** |
| 2 | τ_θz at r = 0.0200 mm, θ = 60° | **+138.0106 MPa** |
| 3 | radius where \|γ\| = 0.0100 on θ = 0° | **0.046894 mm** |
| 3 | radius where \|γ\| = 0.0100 on θ = 90° | **0.025066 mm** |
| 3 | radius where \|γ\| = 0.0100 on θ = 180° (face) | **0.013398 mm** |

Both stresses in part 2 are positive in the stated basis (e_r, e_θ).

## Model chosen and why

The displacement is out-of-plane only (u_z = w(x,y)) with antisymmetric loading, so this is **antiplane shear (mode III)**. The only non-zero stresses are τ = (σ_xz, σ_yz), with γ = (w_,x, w_,y) = ∇w and the pure power-law shear law
γ = (γ₀/τ₀)(τ/τ₀)^(n−1) τ, ν → this is a p-Laplacian problem: ∇·τ = 0 with τ = τ₀γ₀^(−1/n)\|∇w\|^(1/n−1)∇w.

Because the crack flanks are traction-free and parallel to x, a strain-energy density exists (pure power law, monotonic proportional loading), so J is path-independent and an energy argument fixes the exponent. For a separable field τ = A r^(−p) T(θ), γ ∝ τ^n ∝ r^(−np) and W = ∫τ·dγ = n/(n+1) τ·γ ∝ r^(−p(n+1)); then J ∝ r^(1−p(n+1)). A finite non-zero J forces

  p = 1/(n+1) = q = 2/9 = 0.2222222 (n = 3.5).

The displacement then has w = C r^q F(θ), and polar equilibrium
∂τ_r/∂r + (1/r)∂τ_θ/∂θ + τ_r/r = 0
gives the angular ODE

  (F′ H)′ + q(1−q) F H = 0,  H = (q²F² + F′²)^(1/n − 1),

with F(0) = 0 (antisymmetry w(x,−y) = −w(x,y)) and F′(±π) = 0 (traction-free faces, τ_θ ∝ F′H). The stress components are

  τ_r = A r^(−q) q F g^(1/n−1),  τ_θ = A r^(−q) F′ g^(1/n−1),  g = √(q²F²+F′²),
  A = τ₀ γ₀^(−1/n) C^(1/n),  τ₀ = 80 MPa, γ₀ = 0.002.

**How the loading fixes the magnitude.** Substituting the field into the stated J,
J = ∫_Γ (W N_x − t_z w_,x) ds = τ₀ γ₀^(−1/n) C^((n+1)/n) ℐ,
with the dimensionless angular integral ℐ = 2.019595277308 (computed numerically). Setting J = 0.120 N/mm gives

  C = [J / (τ₀ γ₀^(−1/n) ℐ)]^(n/(n+1)) = 9.2560906×10⁻⁴ (mm^(1−q)),  
  A = τ₀ γ₀^(−1/n) C^(1/n) = 64.1932 MPa·mm^q.

The sign is fixed by the upper face moving in +z: F(π) = 1.698433 > 0, so w(r,+π) > 0 and Δw > 0. Then
Δw(r) = 2 C r^q F(π), \|γ\| = C r^(q−1) g(θ) with g(0)=1, g(π/2)=0.614353, g(π)=0.377430.

## Assumption I changed and why

The stated law has no elastic (γ₀/τ₀ is the only compliance), so I used the **pure power law**, not a Ramberg–Osgood form: W = n/(n+1) τ·γ. I also took the **leading local field only**, which the task authorizes for r ≲ 0.30 mm, and used **deformation theory under proportional loading** (the stress at a fixed point is exactly proportional to a fixed tensor here). No superposition was used (the problem is mode III and nonlinear).

## Verification

1. **Eigenvalue.** Shooting the ODE from F(0)=0, F′(0)=1 gives F′(π) = 0 (residual −6×10⁻¹⁷) **exactly at q = 1/(n+1)**; for n = 1 it reproduces F = 2 sin(θ/2), the known mode-III field.
2. **Traction-free faces.** τ_θ(π) = 0 to 1×10⁻¹⁶.
3. **Equilibrium residual** (polar, finite differences) max|R|/scale ≈ 2×10⁻⁸.
4. **Constitutive law.** γ vs (γ₀/τ₀)(τ/τ₀)^(n−1)τ agree in magnitude to 1×10⁻¹⁵ and direction to 4×10⁻¹⁶.
5. **J path-independence.** The independent contour integral J = R∫(W cosθ − τ_r w_,x)dθ gives 0.12000000 on R = 0.01, 0.02, 0.05, 0.20 mm for n = 1, 2, 3.5, 6 (necessary check that the exponent is right).
6. **Known limiting solution.** At n = 1 the code reproduces the linear mode-III relation J = K²/(2μ) (μ = τ₀/γ₀ = 4×10⁴ MPa) and w ∝ √r sin(θ/2); the closed-form Δw(0.02) = 5.5279×10⁻⁴ mm matches the code.

## Does this determine the load at which the slit starts to advance?

**No.** These results give the near-tip field and its amplitude for the *prescribed* energy flux J = 0.120 N/mm, but they contain no fracture criterion, so they cannot by themselves predict the onset load — one still needs a critical value (e.g. J_c) and the outer relation between the applied load and J.

## Missing or uncertain

- The problem supplies J directly; the actual remote load that produces J = 0.120 N/mm is not given, so no onset load can be extracted.
- Finite-strain blunting, the fracture process zone, and any elastic term are omitted (the stated small-strain pure power law and the r ≲ 0.30 mm leading-field instruction).
- The reported stresses assume the stated positive basis and the +z motion of the upper face; reversing the loading reverses the sign of all stresses and of Δw.
- Assumed a straight, stationary slit and plane antiplane conditions, as stated.
