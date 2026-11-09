# Answer

Problem: HRR near-tip field, plane strain, fully plastic total-deformation power
law. n = 6.5, alpha = 0.02, sigma0 = 400 MPa, E = 200 GPa, eps0 = sigma0/E
(= 2e-3). Stationary sharp crack on the negative x-axis, traction-free flanks,
symmetric opening, small strain.

The three requested dimensionless numbers:

1. **p = 0.1333333333**   (exactly 1/(n+1) = 2/15, to well within 1e-5)
2. **H = 1.8935697252**
3. **C = 0.4573057839**

## Method

**Angular (eigenvalue) problem.** The near-tip field is written
sigma_ij = K r^(s-2) sigma~_ij(theta), u_i = (alpha eps0/… ) with the standard
HRR separation; plane-strain effective stress is
sigma_e = sqrt( (3/4)(sr-st)^2 + 3 srt^2 ). Substituting into the Airy/equilibrium
and compatibility equations gives a fourth-order ODE for the angular stress
function f(theta). Shown in terms of f and its derivatives:

    sr = K r^(s-2) (s f + f''),  st = K r^(s-2) s(s-1) f,
    srt = K r^(s-2) (1-s) f'.

Imposing traction-free faces f(pi)=f'(pi)=0 selects the eigenvalue. Rather than a
naive 2-D Newton (which can wander to a spurious branch), I **bracketed f(pi)=0
beginning from the analytic value** `s=(2n+1)/(n+1)` and shot on `f''(0)`. This
returns

    s = 1.866666666667  (matches (2n+1)/(n+1) to 2e-16)
    f''(0) = -0.689256761930,   f(pi)=2.0e-15,  f'(pi)=8.3e-17.

The exponent is therefore **p = 2 - s = 1/(n+1) = 2/15 = 0.1333333333**.

**Displacement angular functions.** With the dimensionless deviatoric
constitutive law e_ij = (3/2) se^(n-1) s_ij, the plane-strain displacement
functions follow from the strain-displacement relations,

    ur = rho*c1/b,  ut = (3 rho srt - ur')/m,
    m = n(s-2),  b = n(s-2)+1,  c1 = (3/4)(sr-st),  rho = se^(n-1),

which is exactly eq. (24) of the knowledge unit's paper.

**J normalisation (eq. 24).** The circular J integral
J = ∫ (W cos(theta) - t_i du_i/dx) r dtheta with W = n/(n+1) sigma_ij eps_ij gives,
for the separated field, J = alpha (eps0/sigma0^n) K^(n+1) I with

    I = 2 ∫_0^pi [ (n/(n+1)) se~^(n+1) cos(theta)
                   - sin(theta)(sr(ut-ur') - srt(ur+ut'))
                   - (n(s-2)+1) cos(theta)(sr ur + srt ut) ] dtheta,

the symmetry factor 2 exploiting the mode-I mirror symmetry. This yields
**I = 4.8183790787** for n=6.5, plane strain.

**H.** From sigma_tt(r,0) = K st0 r^(s-2) and J = alpha eps0/sigma0^n K^(n+1) I
one eliminates K in favour of J. Using s_tt(r,0) = H sigma0 [J/(alpha eps0 sigma0 r)]^p,
the coefficients of sigma0, eps0, alpha all cancel and

    **H = st0 * I^(-p)**,   with st0 = s(s-1) f(0) / max se~ = 2.3352574342.

This gives H = 2.3352574342 * 4.8183790787^(-0.1333333333) = **1.8935697252**,
independent of sigma0, eps0 and alpha magnitudes.

**C.** U(R) = ∫_0^R ∫_{-pi}^{pi} W r dtheta dr = (n/(n+1)) (alpha eps0/sigma0^n)
K^(n+1) R * shape2, with shape2 = ∫_{-pi}^{pi} se~^(n+1) dtheta = 2.5424684099.
Dividing by J R and using J = alpha eps0/sigma0^n K^(n+1) I,

    **C = (n/(n+1)) * shape2 / I** = (6.5/7.5)*2.5424684099/4.8183790787
      = **0.4573057839**,

independent of R and of the amplitude K.

## Verification

**(1) Elastic-limit validation of the whole J machinery.** I fed the linear-elastic
mode-I field (in polar form: sigma_r, sigma_theta, sigma_rtheta, u_r, u_theta as
functions of theta) through the *same* Cartesian code path used for the plastic
field (polar→Cartesian stress rotation, radial traction, finite-difference
Cartesian displacement gradients), and recovered

    J_numeric = 4.550000e-12 = K_I^2 (1-nu^2)/E   (ratio = 1.00000000)

identically at R = 1e-4, 1e-3, 1e-2. This exercises the polar-basis rotation and
the derivative machinery exactly, and pins the factor 2.834 I was chasing to a
**fold bug**: continuing theta in [-pi,0] by taking |theta| without flipping the
sign of the ODD components s_rtheta and u_theta corrupts du_y/dx across the crack
line and inflates/contaminates the numeric J. After fixing the fold the numeric
route and the closed form agree.

**(2) Independent numeric J on the plastic field.** Building the full dimensional
field and computing J by Cartesian finite differences (validated above), with the
fixed symmetry continuation:

    J_numeric = 7.7094e+04  (R-independent for R = 3e-4 … 1e-2)

and the re-derived eq. (24) coefficient agrees, I = 4.81837909.

**(3) Self-consistency of H.** For arbitrary input J (12345.678 J/m^2) I solve for
K, evaluate sigma_tt(r,0), and reconstruct H through the *full dimensional chain*
H = sigma_tt(r,0) / (sigma0 [J/(alpha eps0 sigma0 r)]^p): it returns
1.893569725217 for any r — i.e. the definition sigma_tt = H sigma0 [J/(alpha eps0 sigma0 r)]^p
holds pointwise (ratio 1.00000000 at r = 1e-4, 1e-3, 1e-2), and the measured
exponent from sigma_tt(r,0) at two radii is p = 0.133333333333.

**(4) Invariance checks.** H and C are unchanged under field amplitude A in
{0.5, 1, 2} and contour radius R in {3e-4 … 1e-2}, as required.

**(5) Constitutive consistency.** The angular displacement functions satisfy the
plane-strain flow rule eps_e(u_ang) = se~^n exactly (ratio 50.000000 against
alpha*se~^n with alpha=0.02, i.e. eps_e = alpha eps0 (se/sigma0)^n once the eps0
and the r-exponent are restored), confirming the ũ used in eq. (24) belong to the
same field as the stresses.

### Result (final)

| quantity | value | tolerance |
|----------|-------|-----------|
| p        | 0.1333333333 | ±1e-5 |
| H        | 1.8935697252 | ±0.005 |
| C        | 0.4573057839 | ±0.002 |