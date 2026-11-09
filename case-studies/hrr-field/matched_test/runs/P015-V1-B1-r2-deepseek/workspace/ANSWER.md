# ANSWER

## Requested numbers

All seven with the leading antiplane (mode III) power-law field (equations below).

**1. Relative face slip** `Delta w(r) = w(r,+pi) - w(r,-pi) = 2 A r^s f(pi)`:

| r (mm) | Delta w (mm) |
|---|---|
| 0.0200 | **1.31813e-3** |
| 0.0800 | **1.79369e-3** |

**2. Signed polar shear stresses at r = 0.0200 mm, theta = 60 deg** (with e_r = (cos t, sin t), e_theta = (-sin t, cos t)):

- `tau_rz = tau_r = +38.587 MPa`
- `tau_theta z = tau_theta = +138.011 MPa`

(`|tau| = 143.30 MPa`, `|gamma| = 0.015386` there; both components are positive because w > 0 on the upper side and f, f' > 0.)

**3. Radii where |gamma| = 0.0100** (`|gamma| = A r^(s-1) g(theta)`, with `g = sqrt(s^2 f^2 + f'^2)`, face value taken as theta -> pi^-):

| ray | r (mm) |
|---|---|
| theta = 0 deg | **4.68944e-2** |
| theta = 90 deg | **2.50661e-2** |
| theta = 180 deg (upper face) | **1.33984e-2** |

## Governing model and how the loading fixes its magnitude

The only displacement is out of plane, so this is **antiplane shear (mode III)** with a **pure power-law** constitutive relation. With `gamma = grad w` (engineering shear = gradient of w) and `tau` parallel to `gamma`,

    |tau| = tau0 (|gamma|/gamma0)^(1/n),   tau0 = 80 MPa, gamma0 = 0.002, n = 3.5.

Quasistatic equilibrium `div tau = 0` then gives the p-Laplacian

    div( |grad w|^((1-n)/n) grad w ) = 0,   i.e.  div(|grad w|^(m-2) grad w)=0,  m=(n+1)/n = 9/7.

The leading near-tip field is the mode III (HRR) separable field `w = A r^s f(theta)`, `|grad w| = A r^(s-1) g(theta)`, so `tau ~ r^((s-1)/n)`. A finite, nonzero, r-independent J for the straight crack requires the exponent `1 + (s-1)(n+1)/n = 0`, hence

    s = 1/(n+1) = 2/9 = 0.22222,   stress exponent p = (1-s)/n = 1/(n+1) = 0.22222.

The angular function obeys (with `q=(s-1)/n`, `alpha=(1-n)/n`, `g = sqrt(s^2 f^2+f'^2)`)

    (g^alpha f')' + (1+q) s g^alpha f = 0,

subject to the antisymmetry `f(0)=0` and the traction-free crack face `tau_theta(pi)=0  =>  f'(pi)=0`. With `s=1/(n+1)` the face condition is satisfied automatically by the odd solution (numerically `f'(pi)=1e-16`), which is the mode III counterpart of the J argument. Integrating from `theta=0` with `f'(0)=1` gives

    f(pi) = 1.6984332,   f'(pi/2) = 0.7636948 (so g(pi/2)=0.614353),  g(0)=1, g(pi)=s f(pi)=0.377430.

**Amplitude from J.** The stored work is `W = integral tau.dgamma = n/(n+1) tau |gamma|` on the proportional path. Writing the given J-integral over a circle of radius r (normal N=e_r, ds=r dtheta) and using `|gamma| = A r^(s-1) g` gives

    J = tau0 gamma0^(-1/n) A^((n+1)/n) I_J,
    I_J = integral_{-pi}^{pi} [ n/(n+1) g^((n+1)/n) cos t
          - s f g^((1-n)/n)(s f cos t - f' sin t) ] dt = 2.0195953.

With J = 0.120 N/mm this fixes

    A = ( J / (tau0 gamma0^(-1/n) I_J) )^(n/(n+1)) = 9.2560906e-4  (mm^(1-s)).

The load level enters only through this amplitude A (i.e. through J); the exponent and angular shape do not depend on it.

## Assumptions / changes

- The provided checks and the Hutchinson entry are for **in-plane** J2-deformation-theory mode I/II fields. This problem is **antiplane**: only `u_z=w` is nonzero and the material law involves the out-of-plane shear vector only. I therefore used the antiplane specialization (a scalar p-Laplacian with exponent `(n+1)/n`), not the in-plane Airy/eigenvalue problem. This is not an approximation to the in-plane field; it is the exact mode III asymptotic structure for this constitutive law.
- Faces are taken traction-free and open (mode III sliding is unconstrained), and I used only the leading term (geometry corrections negligible for r <= 0.30 mm as stated).

## Verification

1. **Linear limit (n=1).** The solver returns `f(pi)=2.0` and the analytic result `I_J=pi/4` after normalizing `f'(0)=1/2` (`=pi` for `f'(0)=1`), and the direct J-integral of `w=sqrt(r) sin(theta/2)` reproduces `J=K^2/(2 mu)`; this validates the ODE, the ODE integration and the J-integral reduction.
2. **Energy flux / path independence.** Evaluating the given J-integral numerically from the Cartesian field on circles R = 0.5, 1, 2 mm returns exactly **0.120 N/mm** each (leading field has the exponent that makes J r-independent), confirming the amplitude normalization.
3. **Boundary residual.** The Cartesian finite-difference equilibrium residual `max|d sigma_xz/dx + d sigma_yz/dy| ~ 1.4e-8` on the unit circle (field scale ~ 1e2 MPa), and the face shear `|sigma_yz(theta=pi)| ~ 1e-11` (traction-free). ODE solution checked to be grid/tolerance-independent (`f(pi)`, `I_J` stable over max_step 0.05 -> 5e-4 and rtol 1e-10 -> 1e-13).

## Do these results alone determine the load at which the slit starts to advance?

**No.** They give only the near-tip field and its amplitude (already fixed by the given J = 0.120 N/mm); crack advance additionally requires a fracture criterion (e.g. a critical energy-release rate J_c or a process-zone stress/opening criterion) that is not contained in the local field, and the small-strain separation-of-length-scales conditions needed for J-dominance would also have to be checked.

## Missing / uncertain

- No tabulated mode III HRR reference data were provided to cross-check the angular function at n=3.5; verification rests on the linear-limit analytic check, path-independence of J, and the equilibrium/face residuals.
- The field is a small-strain intermediate asymptotic; very near the tip it becomes self-inconsistent (finite strains / process zone), and the outer geometry and the actual fracture criterion are outside this calculation.
