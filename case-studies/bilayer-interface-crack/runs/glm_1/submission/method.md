# Method

## Problem setting (units)

All results are dimensionless, computed with h = 1 (so H = 2), E2 = 1 (so
E1 = 9), nu = 1/3 for both layers. With this choice |K|sqrt(h)/P,
|K|h^(3/2)/M, G E2 h / P^2 etc. are pure numbers and psi = arg(K h^{i eps})
reduces to arg K (r and h are measured in the same unit).

## 1. eps (exact, closed form)

eps = (1/2pi) ln[(kappa1/mu1 + 1/mu2)/(kappa2/mu2 + 1/mu1)]
    = (1/2pi) ln(2/3) = -0.064531776207.
(kappa = 3 - 4nu = 5/3; mu1 = 9 mu2 gives (kappa mu2 + mu1)/(kappa mu1 + mu2)
= 2/3.)

## 2. G (exact, steady-state energy release)

The prescribed far-behind stress fields (given in the problem statement) are
themselves exact elasticity solutions: x-independent sigma11(y), sigma12 =
sigma22 = 0, plane strain; each layer is separately in equilibrium, its
outer face is traction free, and they exert no tractions on each other. Far
ahead of the tip the bilayer is stress free (that is how the layer-2 loading
is defined). For steady advance of the semi-infinite crack under the *dead*
loads (P, M):

  dW/da = 2 U_behind  (each load moves at its conjugate generalized strain;
                       Clapeyron per unit length),   dU/da = U_behind - 0,
  G = dW/da - dU/da = U_behind,

where U_behind is the strain energy per unit length of the two far-behind
beam states,

  U_behind = N1^2/(2 E1' h) + M^2/(2 E1' I1) + N2^2/(2 E2' H) + M2^2/(2 E2' I2),
  N1 = -P, N2 = +P, M2 = -(M + P(h+H)/2),  Ek' = Ek/(1-nu^2), I1 = h^3/12,
  I2 = H^3/12.

Numerically: G_A E2 h/P^2 = 4/81 + 2/9 + 3/2 = 1.7716049383,
G_B E2 h^3/M^2 = 16/27 + 2/3 = 1.2592592593. This reproduces the classical
moment-loaded DCB result G = 12 M^2/(E' t^3) (two arms of thickness t with
end moments M) and was confirmed by the FEM below.

## 3. |K| (exact, from G)

Interface energy relation (homogeneous limit reproduces G = K^2/E'):

  G = (Lam1 + Lam2) |K|^2 / (16 cosh^2(pi eps)),   Lam_k = (kappa_k+1)/mu_k.

|K|_A sqrt(h)/P = 1.9331267849, |K|_B h^(3/2)/M = 1.6298006013. The FEM
below reproduces these to 0.08% (A), 0.56% (B), 0.25% (combined), which is
the practical cross-check of both the energy argument and the FEM.

## 4. Phase angles psi (FEM)

The actual boundary-value problem is solved directly with a purpose-written
FEM (`fem.py`): plane strain, bilinear quads on a structured mesh graded
geometrically toward the crack tip in x and toward the interface in y;
interface nodes shared for x > 0 (bonded) and duplicated for x < 0 (crack);
the far-beam tractions sigma11(y) of the statement applied as consistent
nodal loads on x = -L (L = 40..60); all other external faces traction-free.
K is extracted by fitting the tractions ahead of the tip,

  (sigma22 + i sigma12)(r) sqrt(2 pi r) e^{-i eps ln r} = K + c1 r + c2 r^2,

using the average of the two elements adjacent to the interface (stresses at
element centres), with r in a window [8 d_min, 0.15 h] (d_min = tip element
size). Mesh refinement d_min = 1e-4 -> 1.5e-5 and fit-window/fit-order
variations change the phases by <= 0.05 deg; phases are reported from the
converged runs:

  psi_A = +52.18 deg,  psi_B = +4.34 deg,  psi_both = +30.48 deg.

## Assumptions

- Linear elastic small-strain plane strain; no crack-face contact (as stated).
- Loads are dead loads (forces/moments prescribed) far behind the tip; the
  layer-2 stresses given in the statement define the (self-equilibrated) end
  loading; far ahead stress-free (both used literally).
- Rigid-body modes pinned arbitrarily (3 dofs) - no effect on stresses.
- FEM: domain truncated at |x| <= L = 40-60 (decay is exponential; far-field
  stresses verified to match the beam fields at x = -20 to <1e-6 relative
  and to vanish at x = +20).

## Checks run

1. eps: exact formula, evaluated independently inside the FEM extractor.
2. G: reproduces classical DCB formula; FEM-consistent (item 5).
3. FEM validation, symmetric DCB (E1 = E2, h = H, pure moment): by mirror
   symmetry psi must be 0 and |K| = sqrt(12): measured psi = 0.003 deg,
   |K| error +0.26% (and decreasing with refinement).
4. FEM validation, homogeneous strip, load case A (mixed mode, no symmetry):
   |K| = sqrt(7) from energy; measured error +0.12%.
5. FEM |K| for the target cases vs exact energy: +0.08% (A), +0.56% (B),
   +0.25% (P = M), i.e. within the 1% tolerance band of the exact values.
6. Energy cross-term (exact): with K = P k_A + M k_B and G quadratic in
   (P, M),  cos(psi_A - psi_B) = (G(1,1) - G_A - G_B)/(2 sqrt(G_A G_B))
   = 0.6695130 => |psi_A - psi_B| = 47.9705 deg exactly; FEM gives
   47.844 deg (0.13 deg difference), confirming both phases jointly.
7. Linearity: FEM |K(1,1)| = |K_A + K_B| to 5 decimals (superposition).
8. psi_both consistency: direct FEM 30.478 deg vs arg(K_A + K_B) built from
   the exact |K| values and the FEM phases 30.422 deg (0.06 deg apart).
9. Mesh/domain/fit-window/fit-order sensitivity of all phase angles
   <= 0.05 deg (see fem_run2.py / fem_run3.py / fem_run4.py output in
   run_log.txt).
10. A completely independent Fourier-spectral solution route for the same
    problem was developed as well (compliance.py + sie_solve.py): its
    ingredients were all validated (per-mode layer compliance against a
    Chebyshev-collocation BVP solver to ~1e-7 relative; Betti symmetry
    A(-xi) = A(xi)^T; the xi -> 0 (beam) and xi -> infinity (bimaterial
    half-plane) limits of A; analytic basis transforms against numerical
    quadrature; exact resultant/pole-matching conditions), but the assembled
    collocation system retained a small systematic misfit (its |K| drifted
    2-10% depending on basis scale), so it was **not** used for the final
    numbers; the FEM (which agrees with the exact energy values) was used
    instead. The layer-compliance module itself is validated and retained.

## Results (results.json)

  eps                     = -0.0645317762
  caseA: G_E2_h_over_P2   =  1.7716049383
          Kabs_sqrt_h_over_P =  1.9331267849
          psi_deg           = +52.18
  caseB: G_E2_h3_over_M2  =  1.2592592593
          Kabs_h32_over_M    =  1.6298006013
          psi_deg           = +4.34
  psi_both_deg            = +30.48

Estimated uncertainties: |K| and G are exact to the digits shown (energy
argument; FEM confirms to better than 0.6%); phase angles +/- 0.1 deg
(refinement spread <= 0.05 deg; residual inconsistency with the exact
energy cross-term relation 0.13 deg).

## What remains unverified

- The phase angles rest on the FEM (converged, and consistent with the exact
  energy constraints to ~0.1 deg), but no closed-form reference exists for
  them at these parameters; an independent reference (e.g. tabulated
  Suo-Hutchinson universal functions) was not available offline.
- The relation G = (Lam1+Lam2)|K|^2/(16 cosh^2(pi eps)) is the standard
  interface result (reduces to K^2/E' for homogeneous materials); it is
  confirmed here at eps = 0 by the DCB test, and at eps != 0 only through
  the agreement of FEM-extracted |K| with the energy values.
- The spectral collocation method (sie_solve.py), despite all validated
  ingredients, was left with an unexplained ~2-10% bias and was not used.
