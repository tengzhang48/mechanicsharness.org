# Reference solution

Run `python V1/solve_reference.py` from the task root, or `python solve_reference.py` inside this folder. The code requires numpy and scipy, prints both independent routes and its checks, and asserts their agreement. `reference_output.json` records a run. No external data or software beyond those packages are needed.

## Results

The upper-face slip and both stresses below are positive under the stated convention.

| Requested quantity | Reference | Unit | Relative tolerance |
|---|---:|---|---:|
| Δw(0.0200 mm) | 0.001318126778 | mm | 2% |
| Δw(0.0800 mm) | 0.001793693738 | mm | 2% |
| τ_rz at 0.0200 mm, 60° | 38.58659686 | MPa | 2% |
| τ_θz at 0.0200 mm, 60° | 138.0106166 | MPa | 2% |
| Radius of γ magnitude = 0.0100 at 0° | 0.04689442305 | mm | 2% |
| Radius at 90° | 0.02506612349 | mm | 2% |
| Radius at the upper face | 0.01339840659 | mm | 2% |

The contour is defined by the magnitude of the engineering-shear vector.

## Route A: equilibrium and an angular initial-value problem

Antiplane equilibrium is ∂_x τ_xz + ∂_y τ_yz = 0, while compatibility is automatic when γ = ∇w. Inverting the given law gives

\[
\boldsymbol\tau=\kappa |\nabla w|^{p-2}\nabla w,
\quad p=1+1/n,\quad \kappa=\tau_0\gamma_0^{-1/n}.
\]

The stored-work density is

\[
W=\frac{\kappa}{p}|\nabla w|^p
 =\frac{n}{n+1}\boldsymbol\tau\cdot\boldsymbol\gamma.
\]

Using the complementary energy instead would incorrectly replace n/(n+1) by 1/(n+1).

Write w = a r^λ f(θ) for the leading symmetric shear-loading field, with w odd in θ. The flux scales as r^[1+p(λ−1)]. Nonzero finite J requires

\[
\lambda=\frac1{n+1}=\frac29.
\]

Thus the energy density is proportional to 1/r and has a finite integral over a disk. Set f(0) = 0 and f′(0) = 1 to fix the arbitrary angular normalization. With q² = λ²f² + f′², equilibrium becomes

\[
\big(q^{1/n-1}f'\big)'
 +\lambda(1-\lambda)q^{1/n-1}f=0,
\qquad f'(\pi)=0.
\]

The last condition is zero face traction. `solve_reference.py` integrates this second-order equation using DOP853; the independently imposed face condition is then a check rather than an initial condition. The computed residual is 2.36 × 10^−16 in the f′(0) = 1 normalization.

The dimensionless flux integral is

\[
I=2\int_0^\pi\left[
 \frac{n}{n+1}q^p\cos\theta
 -q^{p-2}\lambda f\big(\lambda f\cos\theta-f'\sin\theta\big)
 \right]d\theta.
\]

Numerical quadrature gives I = 2.019595277307724. Set

\[
a=\left(\frac{J}{\kappa I}\right)^{1/p}
 =0.0009256090617515132\ {\rm mm}^{1-\lambda}.
\]

Then

\[
\Delta w=2ar^\lambda f(\pi),\qquad
(\gamma_{rz},\gamma_{\theta z})=ar^{\lambda-1}(\lambda f,f'),
\]

and the stress vector follows by the inverse shear law. A contour intersection is

\[
r(\theta)=\left(\frac{a q(\theta)}{0.0100}\right)^{1/(1-\lambda)}.
\]

At r = 0.0200 mm the largest |γ| over angle is 0.01940210, so the requested calculations have modest engineering shear. This does not validate an extrapolation to arbitrarily small r.

## Route B: independent phase parametrization and exact flux

An independent way to solve the angular problem is to introduce 0 ≤ ψ ≤ π/2:

\[
\theta=\psi+\operatorname{atan2}(n\sin\psi,\cos\psi),\quad
D=\cos^2\psi+n^2\sin^2\psi,\quad q=D^{-n/[2(n+1)]},
\]

\[
f=\frac{q\sin\psi}{\lambda},\qquad f'=q\cos\psi.
\]

This route uses a scalar bracketing solve for ψ at a requested θ, with no integration of the angular differential equation. It also permits an exact flux evaluation. On substituting t = tan ψ,

\[
I=2n\int_0^\infty\frac{1+nt^2}{(1+n^2t^2)^2}\,dt
 =\frac{\pi(n+1)}{2n}.
\]

Using this exact I and the phase field reproduces **all seven results** from Route A. Maximum relative disagreement is 4.12 × 10^−16; maximum absolute discrepancy in (f,f′) on a 2,001-point angular grid is 5.39 × 10^−15.

Additional checks:

- Direct dimensional flux integrals at r = 0.005, 0.020, and 0.080 mm give 0.1200000000000002 N/mm to displayed precision.
- In the independently known n = 1 limit, f = 2 sin(θ/2), I = π, and G = τ_0/γ_0. The ODE agrees with f within 1.34 × 10^−15. The standard linear relation J = K_III²/(2G) gives the same face slip as the code; for G = 26,000 MPa, K_III = 35 MPa√mm, r = 0.030 mm, both give 0.000372070622507592 mm.
- The nonzero face value q(π) makes the face contour a regular limit of the field, rather than a stress-free-state assumption. Zero traction on a face does not require every stress component there to vanish.

All requested reference numbers have been independently checked. Floating-point differences in later reruns may affect the last few digits only.

## Scope

The constitutive law and J determine a stationary field. They supply no fracture resistance, cohesive law, critical strain criterion, or damage evolution rule. Therefore they do not determine the initiation load. The stated leading-field assumption supplies outer matching; it does not supply a fracture criterion.
