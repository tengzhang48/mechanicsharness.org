# Key notes for the grader (never given to a solver)

The reference values themselves are in `score.md` (column "Reference"). Parameters of the task: E1 = 9 E2,
nu1 = nu2 = 1/3, h/H = 0.5, plane strain. Units below: h = 1, E2 = 1.

## 1. The closed forms

- Bimaterial constant: mu1 = 27/8, mu2 = 3/8, kappa1 = kappa2 = 5/3, so eps = ln(2/3)/(2 pi) = -0.0645318. It is
  negative because the stiffer layer is above. Dundurs parameters alpha = 0.8, beta = 0.2.
- Energy release rate: the bonded part far ahead of the tip is free of stress, so G is the elastic energy per unit
  length of the two separate layers far behind the tip. With plane-strain moduli Ebar = E/(1 - nu^2), a layer of
  thickness t with axial force N and moment Q about its mid-plane holds N^2/(2 Ebar t) + 6 Q^2/(Ebar t^3). Layer 1
  carries P and M; layer 2 carries P and M + P (h + H)/2. Both layers count.
  Case A: G E2 h/P^2 = 287/162 = 1.771605. Case B: G E2 h^3/M^2 = 34/27 = 1.259259.
- Magnitude of K, for two materials: G = (c1 + c2) |K|^2 / (16 cosh^2(pi eps)), with c = (kappa + 1)/mu for each
  material. Here cosh^2(pi eps) = 25/24. Case A: |K| sqrt(h)/P = 1.933127. Case B: |K| h^(3/2)/M = 1.629801.
  The relation for one material, G = |K|^2/Ebar with an average modulus, omits the cosh factor and gives |K| 2.0 % low.

## 2. The phase angle

- No elementary formula gives it. The elasticity problem around the tip must be solved (finite elements, boundary
  elements, an integral equation), and K read from the solution in a way that respects the oscillating singularity of an
  interface crack: tractions ahead of the tip (sigma22 + i sigma12) sqrt(2 pi r) r^(-i eps) -> K, or the crack opening
  behind the tip with its two-material prefactor and r^(i eps), or an interaction integral with the two-material
  crack-tip field.
- psi = arg(K h^(i eps)) with the reference length h, the thickness of layer 1. With H instead, psi shifts by
  eps ln(H/h) = -2.56 degrees.
- Reference: psi_A = 52.249 degrees (force), psi_B = 4.279 degrees (moment). They come from two independent
  finite-element computations of ours, which agree within 0.02 degree.
- The difference of the two phase angles is exact, from energy alone: K is linear in the loads and G is quadratic in
  them, so cos(psi_A - psi_B) = b/sqrt(a c) for G E2 h = a P^2 + 2 b P M/h + c M^2/h^2. Here a = 287/162, b = 1,
  c = 34/27 and psi_A - psi_B = 47.9705 degrees. A solver may compute one phase angle and obtain the other from this
  relation; that is sound.
- Both loads with M = P h: K h^(i eps) is the sum of the two complex numbers |K_A| e^(i psi_A) and |K_B| e^(i psi_B)
  (per unit P, h = 1). Reference psi = 30.433 degrees. Averaging the two phase angles is wrong (28.26 degrees).
- In the literature the phase angle is written with an angle omega: psi_A = omega and psi_B = omega + gamma - 90
  degrees, gamma = 42.03 degrees here, from the beam constants. Published values of omega that a solver may use
  instead of computing: the printed table of the source paper gives 51.4 degrees at these parameters; a formula for
  ONE material, omega = 52.1 - 3 (h/H) degrees, gives 50.6 degrees; the one-material value for equal thicknesses is
  49.1 degrees. None of these passes: the reference is 52.249 degrees and the check allows 0.5 degree.

## 3. Usual mistakes

| Mistake | Kind | How it shows |
|---|---|---|
| Phase angle from a published or remembered value (51.4, 50.6, 49.1 degrees) instead of a computation | theory | G and the magnitudes pass; the three phase angles fail by 0.8 to 3 degrees, all by the same amount |
| Reference length H (or none) in psi | theory | the three phase angles 2.6 degrees low |
| One-material relation between G and the magnitude of K | theory | both magnitudes 2.0 % low; G passes |
| Plane stress (E for E/(1 - nu^2)) | theory | eps, G and the magnitudes fail |
| Energy of layer 1 only | theory | G of both cases too small |
| K read from the numerical solution with the formulas of one material (no r^(i eps), no two-material prefactor) | theory | phase angles off by an amount that depends on the distance used |
| Sign of K2 (sigma22 - i sigma12, or the materials exchanged) | theory or application, by where it enters | phase angles with the opposite sign |
| Right method, coarse mesh or short model or a poor fit | numerics | phase angles off by a fraction of a degree to a few degrees, and not by the same amount |
| The two load cases exchanged | application | the entries of case A and case B swapped |
| The two loads combined by averaging the phase angles | application | only the phase angle under both loads fails |
| Contact of the crack faces imposed | theory | the task says to ignore it |
