# Reference solution

Run `python V2/solve_reference.py` from the task root, or `python solve_reference.py` inside this folder. The code requires numpy and scipy. It prints two calculation routes plus energy, compatibility, material-law, and limiting checks, and asserts their agreement. A run is saved in `reference_output.json`.

## Results

| Requested quantity | Reference | Unit | Relative tolerance |
|---|---:|---|---:|
| T | 4816.472187 | MPa | 1% |
| C_rr at 60° | +0.003942048238 | mm^1/2 | 2% |
| C_θθ at 60° | +0.0006331781351 | mm^1/2 | 2% |
| C_γrθ at 60° (engineering shear) | +0.005731131133 | mm^1/2 | 2% |
| C_y at 90° | +0.009139018941 | mm^1/2 | 2% |
| Lower radius endpoint on θ = 0 | 0.01000000000 | mm | 2% |
| Upper radius endpoint on θ = 0 | 0.06132522855 | mm | 2% |

The radial interval is nonempty and includes both endpoints under the stated non-strict inequalities.

## Route A: reduce the leading compliance, then calibrate

As stress increases without bound in the formal tip limit, the σ_y/σ_e term contributes bounded strain. It does not contribute to the coefficient multiplying r^−1/2. The remaining leading compliance is

\[
\epsilon_{ij}^{\rm lead}
 =\frac{1+\nu_t}{T}\sigma_{ij}^{\rm lead}
  -\frac{\nu_t}{T}\sigma_{kk}^{\rm lead}\delta_{ij},\qquad
\nu_t=\frac12-\left(\frac12-\nu\right)\frac{T}{E}.
\]

This follows by collecting the deviatoric and spherical parts of the **given tensor law**, not by discarding the elastic volume strain. In particular, the bulk modulus remains E/[3(1−2ν)]. Simply replacing E by T while keeping ν = 0.30 is incorrect.

The resulting leading plane-stress crack problem is linear isotropic elasticity with modulus T and Poisson ratio ν_t. Its traction-free tensile field has r^−1/2 stresses and r^1/2 displacements. Let K_loc denote its stress-intensity coefficient in MPa√mm. Evaluating the local flux gives

\[
K_{\rm loc}=\sqrt{JT},\qquad J=K_{\rm loc}^2/T.
\]

For c = cos(θ/2), s = sin(θ/2), define the stress coefficients

\[
\sqrt r(\sigma_{rr},\sigma_{\theta\theta},\sigma_{r\theta})
 =\frac{K_{\rm loc}}{\sqrt{2\pi}}\big(c(1+s^2),c^3,sc^2\big).
\]

On θ = 0, σ_rr = σ_θθ and σ_rθ = σ_zz = 0, so

\[
C_{zz}(0)=-2\nu_t\sqrt{\frac{J}{2\pi T}}.
\]

The right side is monotonic in 0 < T < E and yields a unique admissible inverse. For an explicit stable calculation, put

\[
L=-C_{zz}/\sqrt{J/(2\pi)},\quad b=(1-2\nu)/E,\quad
\sqrt T=\frac{2}{L+\sqrt{L^2+4b}}.
\]

The solution is T = 4816.472187302739 MPa and ν_t = 0.4862386508934207. Use

\[
C_{rr}=(S_{rr}-\nu_tS_{\theta\theta})/T,\quad
C_{\theta\theta}=(S_{\theta\theta}-\nu_tS_{rr})/T,\quad
C_{\gamma r\theta}=2(1+\nu_t)S_{r\theta}/T,
\]

where S is the stress-coefficient vector above. The conventional compatible displacement field is

\[
\frac{(u_x,u_y)}{\sqrt r}
 =\frac{K_{\rm loc}}{2\mu_t\sqrt{2\pi}}
  (\kappa_t-\cos\theta)(\cos(\theta/2),\sin(\theta/2)),
\]

with μ_t = T/[2(1+ν_t)] and κ_t = (3−ν_t)/(1+ν_t). At θ = 90° this gives

\[
C_y=(3-\nu_t)\sqrt{J/(4\pi T)}.
\]

## The two radius screens

On θ = 0, the three leading principal-strain coefficients are

\[
\sqrt{\frac{J}{2\pi T}}\,(1-\nu_t,\ 1-\nu_t,\ -2\nu_t).
\]

At the calibrated ν_t, the absolute out-of-plane coefficient is the largest. Therefore the strain screen gives

\[
r\ge (|C_{zz}|/0.05)^2=0.0100000\ {\rm mm}.
\]

The effective stress ahead is K_loc/√(2πr), so the high-stress screen gives

\[
r\le\frac{JT}{2\pi(10\sigma_y)^2}
 =0.06132522855\ {\rm mm}.
\]

This is more restrictive than r ≤ 0.200 mm. These are tests of the **leading field** as requested. They do not certify a 2% physical error at finite r; at the upper endpoint the dropped σ_y/σ_e term is still 0.1. A calculation with the complete finite-r solution can give different screened endpoints and is a different target from this problem.

## Route B and independent checks

The second route calibrates T by a bracketing root solve using the original 3D tensor compliance with only the bounded intercept term removed. It rotates Cartesian stresses and strains into the requested polar basis, finds the radial bounds by root finding on the stress invariant and the eigenvalues of the full strain tensor, and obtains ν_t from a separate uniaxial compliance experiment. It does not use the explicit calibration quadratic or the closed radial-bound formulas.

| Check | Result |
|---|---|
| Two routes, all seven requested values | Maximum relative difference 1.40 × 10^−14 |
| Direct J using W = σ:ε/2 and finite-difference u_,x at r = 0.015, 0.040, 0.100 mm | 0.800000000000063, 0.799999999999311, 0.800000000000098 N/mm |
| sym(∇u) versus original leading tensor compliance at r = 0.040 mm, θ = 60° | Maximum absolute strain difference 2.21 × 10^−12 |
| Complete law in uniaxial tension at 30 MPa | 0.00429527393719145, equal to σ_y/E + (30−σ_y)/T to 10^−18 |
| T = E special material limit | Original elastic tensor recovered exactly in the numerical check; ν_t becomes ν |

The direct J check validates the normalization used by both inverse routes. The strain/displacement check independently validates C_y and the angular strains. The complete-law uniaxial check verifies that T is the actual tensile slope, not a secant modulus or an independently assumed shear modulus. All requested reference values have been checked.

## Scope

The interval supplies a possible working range for the specified approximation and kinematic screens. It does not eliminate the need for matching, a continuum cutoff, or finite-strain corrections outside that range. Since no fracture resistance or damage/critical-separation law is given, crack initiation is not determined. A solver should complete the requested inverse calculation before making this limitation clear.
