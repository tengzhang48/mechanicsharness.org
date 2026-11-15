# Key notes for the grader (never given to a solver)

The reference values themselves are in `score.md` (column "Reference"). Parameters of the task: E1 = 9 E2,
nu1 = nu2 = 1/3, h/H = 0.5, plane strain. Units below: h = 1, E2 = 1. In this version the task gives no formula for the
crack-tip field; the solver had to supply it.

## 1. The crack-tip field and the closed forms

- The crack-tip field of a crack on the interface of two materials: the traction ahead of the tip is
  T(r) = sigma22 + i sigma12 = K r^(i eps) / sqrt(2 pi r), with K complex and the bimaterial constant
  eps = (1/2 pi) ln[(kappa1/mu1 + 1/mu2)/(kappa2/mu2 + 1/mu1)], kappa = 3 - 4 nu (plane strain). Here mu1 = 27/8,
  mu2 = 3/8, kappa1 = kappa2 = 5/3, so eps = ln(2/3)/(2 pi) = -0.0645318 (negative because the stiffer layer is above).
  A solver may write eps with the opposite convention (materials exchanged, or sigma22 - i sigma12); what is graded is
  the resulting numbers.
- So the phase of T turns with r: arg T = arg K + eps ln r. When r increases tenfold it changes by eps ln 10 radians,
  here -8.514 degrees. The amplitude lim sqrt(2 pi r) |T| is |K|. The phase at r = h is arg(K h^(i eps)).
- Energy release rate: the bonded part far ahead of the tip is free of stress, so G is the elastic energy per unit
  length of the two separate layers far behind the tip, with plane-strain moduli Ebar = E/(1 - nu^2). Layer 1 carries
  P and M; layer 2 carries P and M + P (h + H)/2. Case A: G E2 h/P^2 = 287/162 = 1.771605. Case B: 34/27 = 1.259259.
- Amplitude from G, for two materials: G = (c1 + c2) |K|^2 / (16 cosh^2(pi eps)), c = (kappa + 1)/mu. Here
  cosh^2(pi eps) = 25/24. Case A: 1.933127. Case B: 1.629801. The relation of one material omits the cosh factor
  and gives the amplitude 2.0 % low.

## 2. The phase at r = h

- No elementary formula gives it. The elasticity problem must be solved (finite elements, boundary elements, an
  integral equation) and the crack-tip field read from the solution in a way that respects the turning of the phase:
  T(r) sqrt(2 pi r) r^(-i eps) -> K ahead of the tip; or the crack opening behind the tip with its two-material
  prefactor and r^(i eps); or an interaction integral with the two-material crack-tip field. A solver may also find the
  turning numerically, from the slope of arg T against ln r in its own solution; that is sound if the slope is
  -0.148588 rad per unit of ln r (eps).
- Reference: phase at r = h, 52.249 degrees (force), 4.279 degrees (moment). Two independent finite-element
  computations of ours agree within 0.02 degree.
- The difference of the two phases is exact, from energy alone: cos(psi_A - psi_B) = b/sqrt(a c) for
  G E2 h = a P^2 + 2 b P M/h + c M^2/h^2, a = 287/162, b = 1, c = 34/27, so psi_A - psi_B = 47.9705 degrees.
- Both loads with M = P h: the sum of |K_A| e^(i psi_A) and |K_B| e^(i psi_B): 30.433 degrees.
- Published values that a solver may use instead of computing: the printed table of the source paper gives
  omega = 51.4 degrees here (psi_A = omega); a formula for one material gives 50.6; one material with equal thicknesses
  49.1. None passes.

## 3. Usual mistakes

| Mistake | Kind | How it shows |
|---|---|---|
| The tip taken as that of one material: no turning of the phase | theory | change per decade 0; the phases at r = h off by eps ln(r_read/h), e.g. +25.5 degrees if read at r = h/1000 |
| The turning known, the phase reported where it was read, not at r = h | theory | the three phases off by the same amount; the change per decade right |
| The phase reported at r = H | theory | the three phases 2.56 degrees low |
| The change per factor e (eps in degrees, -3.70) or eps itself (-0.0645) in place of the change per decade | application | only the change per decade fails |
| Sign of the turning reversed (materials exchanged in eps) | theory or application, by where it enters | change per decade +8.51; the phases at r = h off |
| Relation of one material between G and the amplitude | theory | both amplitudes 2.0 % low |
| Plane stress | theory | G and the amplitudes fail |
| Energy of layer 1 only | theory | G of both cases too small |
| Right method, coarse mesh or poor fit | numerics | phases off by a fraction of a degree to a few degrees, not all by the same amount |
| The two loads combined by averaging the phases | application | only the phase under both loads fails |
