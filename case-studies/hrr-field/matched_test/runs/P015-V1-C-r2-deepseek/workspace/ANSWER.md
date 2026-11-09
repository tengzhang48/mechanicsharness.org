# Slip and strain contours near a stationary slit — answer

Units: mm, N, MPa. Antiplane (mode III) shear, displacement $u_z=w(x,y)$ only, slit on the
negative $x$-axis, tip at the origin, material $-\pi<\theta<\pi$.

## 1. Relative face slip

$$w = A\,r^{\beta}f(\theta),\qquad \beta=\frac{1}{n+1}=\frac{1}{4.5}=0.222222,\qquad A=1.2874\times10^{-3}$$
(amplitude fixed by $J$, below). With $f$ odd in $\theta$, $w(r,-\pi)=-w(r,\pi)$ and
$\Delta w(r)=w(r,+\pi)-w(r,-\pi)=2A\,r^{\beta}f(\pi)$, where $f(\pi)=1.69843$
(exact: $f(\pi)=(n+1)n^{-n/(n+1)}$).

| r (mm) | Δw (mm) |
|--------|---------|
| 0.0200 | **1.83×10⁻³** |
| 0.0800 | **2.49×10⁻³** |

## 2. Signed stresses at r = 0.0200 mm, θ = 60°

$\tau_{rz}=+42.4$ MPa, $\tau_{\theta z}=+151.7$ MPa (magnitude 157.5 MPa). Both components are
positive because $\gamma_r$ and $\gamma_\theta$ are both positive there.

## 3. Radii where |γ| = 0.0100

$|\gamma|=A\,r^{\beta-1}Q(\theta)$, $Q=\sqrt{\beta^2f^2+f'^2}$, so
$r=\big(A\,Q(\theta)/0.01\big)^{1/(1-\beta)}$.

| ray | r (mm) |
|-----|--------|
| θ = 0°    | **0.07167** |
| θ = 90°   | **0.03831** |
| θ = 180° (face limit) | **0.02048** |

## Governing model and how the loading fixes its magnitude

The loading is antiplane and antisymmetric, so this is a mode-III crack problem. The constitutive
law $\gamma=\gamma_0(\tau/\tau_0)^n\hat\tau$ (magnitude $\gamma=\gamma_0(\tau/\tau_0)^n$) inverts to
$\tau\propto|\nabla w|^{1/n-1}\nabla w$, making equilibrium $\nabla\!\cdot\!\tau=0$ the
$p$-Laplacian $\nabla\!\cdot\!\big(|\nabla w|^{m-1}\nabla w\big)=0$ with $m=1/n$. The self-similar
near-tip field is $w=A r^{\beta}f(\theta)$. Path-independence/finiteness of $J$ for a single power
law forces the HRR exponent $\beta=1/(n+1)$; traction-free faces require $f'(\pm\pi)=0$, and the
loading symmetry gives $f$ odd, fixing the angular problem
$$(R f')'+D\beta R f=0,\quad R=(\beta^2f^2+f'^2)^{(m-1)/2},\quad D=1+(\beta-1)/n=n/(n+1).$$
The single remaining scale $A$ is fixed by the supplied energy flux $J=0.120$ N/mm via
$J=\alpha A^{1+1/n}$ (computed below), since $J$ is the only load measure in the local problem.

## Verification

- **Boundary residual:** with $\beta=1/(n+1)$, integrating the angular ODE from $f(0)=0,f'(0)=1$
  gives $f'(\pi)=1.4\times10^{-16}$, i.e. the traction-free condition $\tau_{\theta z}(\pm\pi)=0$
  is satisfied to machine precision (residual $<5\times10^{-8}$ MPa).
- **Energy flux / path independence:** $J=\int_\Gamma(WN_x-t_zw_{,x})ds$ was evaluated directly at
  $r=0.02,\ 1,\ 2$ mm with $A=1$: all give 624.136, ratio 1.0000000, confirming path independence
  and confirming the exponent. Setting $J(A)=624.136\,A^{1.285714}=0.120$ gives
  $A=1.2874\times10^{-3}$.
- **Closed form:** $f(\pi)=(n+1)n^{-n/(n+1)}=1.69843$ matches the numerical integration, and for
  $n=1$ the model reduces to the classical linear mode-III field ($f\propto\sin(\theta/2)$).

## Assumptions changed

- Treated the problem as the leading (HRR) near-tip field and used the stated statement that
  geometry outside $r=0.30$ mm is negligible; no far-field/loaded-boundary condition was imposed.
- No additive rigid translation (as stated); this sets $w(\text{tip})=0$ and the odd branch.
- $J$ is taken as the only load magnitude (the applied load level itself is not specified).

## Does this determine the load at which the slit starts to advance?

**No.** These results give the local fields for the given $J=0.120$ but include no fracture/toughness
criterion (e.g. a critical $J_c$, crack-tip slip, or critical strain) and no relation between the
remote applied load and $J$, so the onset load cannot be determined from them alone.

## Missing / uncertain

- The far-field loading pattern is not given, so only the amplitude-normalized local field is known;
  absolute stresses/slips depend on the assumption that $J$ alone sets the scale (correct for a
  homogeneous power-law medium).
- Numerical input: angular ODE solved with `scipy` DOP853 (rtol 1e-12); J by trapezoidal
  quadrature. Estimated numerical error $\ll$ the 2% tolerance.
- All seven numbers are quoted to 3+ significant figures, well within the 2% tolerance.
