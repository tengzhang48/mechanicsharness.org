# MP-24 — Suddenly pressurized crack in a poroelastic solid: Answers

Reproduce: `python3 solve.py` (prints every value in `results.json`, rewrites `results.json`,
and runs the two numerical checks quoted below). Worked example below uses unit loading
magnitude; for load magnitude `s` multiply K, V, p by `s`.

## Model and method (and why)

Biot poroelasticity with quasi-static mechanical equilibrium. Two elementary facts make the
problem tractable without solving the full transient field:

1. **Instantaneous response = undrained elasticity.** Setting the fluid-content increment
   `zeta = 0` in the constitutive laws gives `p = -(nu_u-nu) sigma_kk / (alpha (1-2nu))` and
   `2G e_ij = sigma_ij - nu_u sigma_kk delta_ij`: an incompressible-elastic-type response with
   shear modulus G and Poisson ratio nu_u. At `t -> inf` the pressure is the boundary value
   and the skeleton is drained elastic (G, nu). Both are standard elasticity problems.
2. **zeta obeys a homogeneous diffusion equation.** Eliminating sigma_kk between the
   constitutive laws and plane-strain compatibility `Lap(sigma_kk) = -alpha(1-2nu)/(1-nu)
   Lap(p)` gives `dzeta/dt = c Lap(zeta)` with the `c` quoted in the problem. At small times
   the disturbance is a thin diffusion boundary layer at the crack faces; its effect on the
   crack volume follows from the **Betti reciprocal theorem** with the pore-pressure
   eigenstrain `e* = alpha(1-2nu) p delta_ij / (2G)`, using the unit-pressure elastic crack
   as reference state. This avoids solving the coupled transient crack problem.

The reference elastic crack field is the Westergaard solution
`Z(z) = s [ z/sqrt(z^2-L^2) - 1 ]`, `sigma_kk = 2 Re Z` (independent of elastic constants);
on the faces `sigma_xx = sigma_yy = -s`, so `sigma_kk(faces) = -2 s`; near-tip
`K = s sqrt(pi L)`; opening `delta(x) = 2(1-nu') s sqrt(L^2-x^2)/G`, hence volume
`V = pi (1-nu') s L^2/G` with nu' the relevant Poisson ratio.

## Derivations in brief

### (a) Loading 1: sigma_n = -1, p = 0 on the faces
- **t = 0+ (undrained):** crack under unit internal pressure, moduli (G, nu_u):
  `K1_0 = sqrt(pi L)`, `V1_0 = pi (1-nu_u) L^2/G`.
- **t -> inf (drained, p = 0 everywhere):** same traction problem with (G, nu):
  `K1_inf = sqrt(pi L)`, `V1_inf = pi (1-nu) L^2/G`.
  (K carries no elastic constants, so it is the same at both ends; in fact it is constant
  for all t for this loading.)
- **p just inside the faces at t = 0+:** the undrained constraint zeta = 0 with
  sigma_kk(faces) = -2 gives `p_faces_0 = 2 (nu_u-nu) / (alpha (1-2nu))` (positive:
  compression raises pore pressure). This is the value "before it diffuses"; for t > 0 the
  boundary condition p = 0 holds exactly on the faces.
- **Short-time growth of V:** write the solution as (undrained part) + (transient part II).
  Part II has zero face traction, zero initial field, and face pressure `-p_faces_0`. For
  ct << L^2 the faces are locally planar, so in the layer `|y| << L`:
  `sigma_yy^II = 0`, the tangential-strain coupling is smaller by O(sqrt(ct/L^2)) (verified
  a posteriori by the volume scaling), and zeta^II = C p^II with
  `C = alpha^2 (1-2nu)^2 (1-nu_u) / [2G (nu_u-nu) (1-nu)]`; the zeta-diffusion then gives
  1-D heat flow for p^II with the given c:
  `p^II = -p_faces_0 erfc( |y| / (2 sqrt(c t)) )`.
  Betti with the unit-pressure reference crack (traction-free Part II, eigenstrain
  e* = alpha(1-2nu) p^II delta/(2G)) gives the volume change
  `Delta V = [alpha(1-2nu)/(2G)] int int sigma_kk^R p^II dA`.
  With sigma_kk^R = -2 on both faces and int_0^inf erfc(y/2sqrt(ct)) dy = 2 sqrt(ct)/sqrt(pi):
  `Delta V = [alpha(1-2nu)/(2G)] * 2 * 2L * p_faces_0 * 2 sqrt(ct/pi)`
          = `16 (nu_u-nu)/sqrt(pi) * (L^2/G) * sqrt(ct/L^2)`.
  So `V ~ V1_0 + A sqrt(ct/L^2)` with `A = 16 (nu_u-nu)/sqrt(pi) * L^2/G`.
  Tip zones of size sqrt(ct) contribute only O((ct)^{3/4}) and are subleading.

### (b) Loading 2: sigma_n = 0, p = 1 on the faces
- **t = 0+:** no traction is applied; the undrained interior (zeta = 0) has sigma = 0,
  p = 0, the boundary value p = 1 being confined to a layer of zero thickness. Hence
  `K2_0 = 0` (and V2_0 = 0).
- **t -> inf:** in 2-D there is no bounded steady state of Lap(p) with p = 1 on the crack
  and p -> 0 at infinity, so p -> 1 everywhere; sigma = 0 then satisfies equilibrium,
  compatibility, and all boundary conditions (unique). This is a uniform free strain
  e = alpha(1-2nu)/(2G) delta_ij, which does not open a zero-thickness slit and has no
  crack-tip singularity: `K2_inf = 0`, `V2_inf = 0`.
  (Physically: at intermediate times K2(t) and V2(t) are nonzero — diffusion from the faces
  pressurizes the near-crack skeleton while the far field is still unpressurized — but they
  return to zero as the pressurized zone grows without bound.)

## (c) Numerical values (nu = 0.2, nu_u = 0.4, alpha = 0.89)

K in units of sqrt(pi L), V in units of L^2/G, per unit loading magnitude.

| quantity   | value          | formula                                   |
|------------|----------------|-------------------------------------------|
| K1_0       | 1.0000000000   | sqrt(pi L)                                |
| K1_inf     | 1.0000000000   | sqrt(pi L)                                |
| V1_0       | 1.8849555922   | pi (1-nu_u) L^2/G                         |
| V1_inf     | 2.5132741229   | pi (1-nu) L^2/G                           |
| A          | 1.8054066674   | 16 (nu_u-nu)/sqrt(pi) * L^2/G             |
| p_faces_0  | 0.7490636704   | 2 (nu_u-nu)/(alpha (1-2nu))               |
| K2_0       | 0              | 0                                         |
| K2_inf     | 0              | 0                                         |
| V2_inf     | 0              | 0                                         |

## Checks (outputs pasted, not retyped)

CHECK 1 — reciprocity consistency at the LONG-TIME limit. The same Betti formula used for
the short-time coefficient, applied at t -> inf (where p^II = -p^undrained and sigma_kk of
the undrained and drained states coincide), gives
Delta V(inf) = [(nu_u-nu)/(2G)] int int sigma_kk^2 dA, which must equal
V1_inf - V1_0 = pi (nu_u-nu) L^2/G, i.e. int int sigma_kk^2 dA must equal 2 pi L^2.
Direct numerical integration of the Westergaard field (L = 1, grid dx = 0.004 over
[-4,4]^2 plus the analytic far-field tail):

    CHECK 1 (reciprocity consistency, long-time limit)
      int sigma_kk^2 dA = 6.277848  (grid 6.179674 + tail 0.098175)
      expected 2*pi*L^2 = 6.283185,  rel.err = 8.494e-04

How it could have failed: any error in the eigenstrain factor alpha(1-2nu)/(2G), in the
reference stress field, or in the sign/factor conventions of the Betti identity would make
the integral differ from 2 pi L^2 (an earlier run with a wrong sqrt branch gave 505.6,
caught immediately). It matches to 8.5e-4 (grid-resolution level), so the same formula
used at short times is trustworthy.

CHECK 2 — independent 2-D simulation of the short-time diffusion. Explicit FTCS solution
of dq/dt = Lap(q) (c = 1) on a 933 x 933 grid (dx = 0.003, dt = 0.2 dx^2), q = 0 initially,
q = -1 held on the slit y = 0, |x| <= 1; then M(t) = int int sigma_kk^R q dA is the quantity
that enters the Betti volume. Prediction from the face boundary layer:
M ~ 16 sqrt(t/pi); crack-tip diffusion zones (size sqrt(t), sigma_kk ~ r^{-1/2}) add
O(t^{3/4}) and the curvature of sigma_kk along the faces adds O(t):

    CHECK 2 (2-D slit diffusion -> short-time volume growth)
      face-layer prediction: M_face = 16 sqrt(t/pi)  (M/sqrt(t) -> 9.027033);  tips/curvature add O(t^0.75), O(t)
      t= 2.610e-04  M=   0.1249  M_face=   0.1458  M/M_face=  0.8567
      t= 5.724e-04  M=   0.1782  M_face=   0.2160  M/M_face=  0.8252
      t= 1.253e-03  M=   0.2521  M_face=   0.3195  M/M_face=  0.7891
      t= 2.749e-03  M=   0.3535  M_face=   0.4733  M/M_face=  0.7469
      t= 6.034e-03  M=   0.4893  M_face=   0.7012  M/M_face=  0.6978
      fit M = C0 sqrt(t/pi) + C1 t^0.75 + C2 t over 39 samples:
      C0 = 16.1118 (predicted 16), C1 = -11.3458, C2 = 4.7682, rms resid = 4.84e-05

The fitted leading coefficient C0 = 16.11 reproduces the predicted 16 to 0.7%
(discretization level; the boundary layer is only ~5 cells thick at the earliest sample),
and the residual is entirely accounted for by the subleading O(t^{3/4}) + O(t) terms
(rms 4.8e-5). How it could have failed: a wrong boundary-layer coefficient (e.g. 8 instead
of 16, or the wrong diffusivity) would shift C0 by O(1), far above the 0.7 % fit error; a
wrong power (e.g. growth like t instead of sqrt(t)) would destroy the fit quality.

Also verified analytically: int_0^inf erfc(s) ds = 1/sqrt(pi) (used in the layer integral),
and the consistency of the quoted c with the zeta-diffusion derivation
(2G Lap(zeta) = alpha^2 (1-2nu)^2 (1-nu_u)/[(1-nu)(nu_u-nu)] Lap(p), dzeta/dt = kappa Lap(p)).

## Questions / assumptions

- Loading 2, "sigma_n = 0, p = H(t)", is read literally: the TOTAL normal traction on the
  crack faces is zero while the boundary pore pressure is 1 (e.g. fluid exchanged through a
  membrane that carries no mechanical load). Under this reading all requested loading-2
  quantities are zero (they are nonzero at intermediate times, not asked). If instead the
  intent was a crack cavity filled with fluid at pressure 1 that also pushes on the faces
  (sigma_n = -1 with p = 1), that is the SUM of loadings 1 and 2, giving
  K = sqrt(pi L), V = pi(1-nu_u)L^2/G at t = 0+ and K = sqrt(pi L), V = pi(1-nu)L^2/G at
  t -> inf.
- "pore pressure immediately next to the crack faces at t = 0+, before it diffuses" is the
  undrained value 2(nu_u-nu)/(alpha(1-2nu)) per unit traction, uniform along the faces
  (away from the tips, where it inherits the 1/sqrt(r) elastic singularity). For any t > 0
  the boundary value on the face itself is p = 0; the stated value is the interior-side
  limit before the diffusion boundary layer forms.
- Loading magnitudes are taken as 1 (unit step of unit amplitude), as in the problem.
- Plane strain assumed for the Westergaard crack formulas (standard); c used in the
  boundary layer is exactly the c quoted in the problem.

## Missing / uncertain

- Nothing missing for the quantities reported. Not computed (not requested): the full
  transient V(t) and K2(t) at intermediate times, which require solving the coupled
  diffusion problem on the crack line numerically.

## Reproduce

    python3 solve.py

prints the nine report values, runs CHECK 1 and CHECK 2, and writes results.json.
