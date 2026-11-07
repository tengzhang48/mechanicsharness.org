# PROVISIONAL_REPORT

- UTC: 2026-10-15T13:15:57Z
- Phase: checkpoint_report (provisional). Budget continues after this.

## Problem
Periodic through-thickness cracks in an isotropic plane-stress elastic film
(thickness h, modulus E, Poisson nu, prestress sigma0 I) bonded to an elastic
foundation of stiffness k (force per unit area -k u). Energy release G per
advancing crack area; G = (KI^2+KII^2)/E. Two geometries:
- Case A: semi-infinite cracks, period H, tips at x1=0, all advance.
- Case B: existing infinite cracks period 2H + additional semi-infinite cracks
  period 2H (offset H), only additional tips advance.
Sampling at ell/H = 0.5, 1, 1.5, 2.

## What I have DERIVED and can state with confidence

1. Plane-stress dispersion (from derive3.py, symbolic determinant of the
   2x2 equilibrium symbol matrix). With C = E/(1-nu^2):
     s1^2 = q^2 + k/(h C)                 (shear-like branch)
     s2^2 = q^2 + 2 k/(h C (1-nu))        (dilatational-like branch)
   Slowest x1-decay at q=0 (branch 1), so the wake boundary-layer length is
     ell = sqrt(h C / k) = sqrt(E h / (k (1-nu^2))).
   The second length sqrt(h C (1-nu)/(2k)) ~ 0.59 ell (nu=0.3) decays faster.

2. With the foundation, there is NO x1-independent stress-free relaxed wake:
   the only bounded x1-independent equilibrium is u=0, sigma = sigma0 I
   (derive8.py). So the crack disturbance returns to the reference state far
   ahead AND far behind, and the whole effect is a boundary layer of width
   ~ell near each tip. The foundation length ell is thus the natural decay
   length requested in part 1.

3. Mode type: the geometry + equibiaxial loading are symmetric about the
   crack plane and the tip-reflection, so the near-tip field is Mode I
   (KII = 0 expected; supported by FEM field s12 -> 0 along the wake, but
   not yet proven symbolically).

4. FEM ground-truth scaffold (fem_j.py, fem_energy2.py) solves the plane-
   stress strip + foundation with prestress applied as crack-face tractions
   t = -sigma0 n on the free faces. Verified qualitatively:
   - Traction-free faces hold: far-behind wake has sigma22 -> ~0 at x2=0,H.
   - Ahead-crack-plane total sigma22 recovers to sigma0 (~1.0) within ~ell.
   - Opening has a finite plateau behind the tip (boundary layer), consistent
     with a foundation-limited wake.

## Checks ACTUALLY executed
- Fem: far behind (x1=-18ell) crack-face sigma22 = 0.0029 ~ 0 (traction free
  OK); crack-plane sigma11 ~ 0.735 (partially relaxed wake).
- Fem energy scan (nx=400): Pi decreases monotonically as crack lengthens
  (shift -4..+4): Pi = -8.4174, -8.3537, -8.2899, -8.2262, -8.1624.
  => G > 0 (physically correct sign).
- A first extractor (fem_j.py J-integral) was NOT path-independent
  (J = 0.016, 0.0015, -0.0014 for r=0.3..1.2) -> REJECTED as unreliable.
- An earlier closed-form fracture integral I wrote (compute.py Finf) was
  WRONG (gave negative/NaN); it is discarded and must NOT be used.

## Incomplete / not yet established
- Numeric G_A(ell/H), G_B(ell/H), normalized gamma and KI at the 4 sampling
  points: NOT yet produced reliably. The energy-difference FEM extractor
  (fem_energy2.py) was interrupted before a converged value; a coarse run
  suggested gamma_A(ell/H=1) ~ 0.5-0.6 (G ~ 0.5-0.6) but this is UNCONVERGED
  and unverified.
- Case B has not been computed at all.
- KII = 0 is assumed from symmetry, not proven.
- The H/ell -> 0 and H/ell -> infinity limits: not yet derived.
- No ANSWER.md / final results.json yet; results.json currently contains the
  DISCARDED wrong values from compute.py (must be overwritten).

## Uncertainty
- Whether placing ell = sqrt(Eh/(k(1-nu^2))) (branch-1) is the intended ell;
  the problem says "the decay length ell of the crack-induced wake field",
  which matches branch 1 at q=0.
- Energy-difference FEM needs mesh convergence and a larger tip-shift to beat
  discretization noise; the dPi over 2 elements at dx=0.05 was only ~0.06.
- Case B accounting (only new tips advance; old cracks fixed) needs care.

## Next (continued analysis phase)
1. Converge the FEM energy-difference extractor for G_A at multiple ell/H,
   cross-check with an independent crack-face work integral.
2. Build Case B FEM and extract G_B.
3. Derive/confirm KI and the two limits.
4. Write final ANSWER.md and overwrite results.json with verified values.