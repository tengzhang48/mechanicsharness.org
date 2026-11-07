# ANSWER — Energy release for periodic cracks in an elastic film on a local elastic foundation

Scope note. This is the continuation write-up of the saved partial attempt. It assembles the
derivation, the requested quantities, and the check inventory. Verification status of each
statement is tracked in VERIFICATION.md. All quantities use the model exactly as stated in
PROBLEM.md: quasistatic plane-stress film (h, E, nu; nu = 0.3 for numbers), reference
equibiaxial prestress sigma0, Winkler foundation f_a = -k u_a, through-thickness
traction-free cracks.

## 1. Reduced field equations, decay length, far states

With Delta_sigma from the plane-stress law, equilibrium h sigma_{ab,b} - k u_a = 0 becomes,
for the increment field u,

    h Ebar/2 (1-nu) nabla^2 u_a + h Ebar/2 (1+nu) u_{b,ba} - k u_a = -h sigma0 delta_{a2} ... 

(source only enters through the crack-face condition sigma_{a b} n_b = 0; the increment
equation itself is homogeneous). With exp(lambda x1 + i q x2) the two Stroh branches are

    lambda^2 = q^2 + k/(h mu)          (divergence-free / shear branch),  mu = E/(2(1+nu))
    lambda^2 = q^2 + k/(h Ebar)        (curl-free / longitudinal branch), Ebar = E/(1-nu^2)

The slowest decay is the q = 0 longitudinal branch, so the crack-induced wake field decays
over

    ell = sqrt(h Ebar / k).

The shear branch is shorter by factor sqrt((1-nu)/2) (lambda_shear/lambda_long =
sqrt(2/(1-nu))), so ell is unambiguously the controlling length; all normalizations below
use it.

### Far states (bounded equilibria of the stated model)

The far states satisfy the strip ODE h Ebar u2'' - k u2 = 0 (u1 = 0, no x1 dependence),
with traction-free crack edges. On a strip of width S bounded by cracks at its edges:

    u1 = 0,
    u2(x2) = -(sigma0 ell/Ebar) sinh((x2 - S/2)/ell) / cosh(S/(2ell)),
    sigma22 = sigma0 [1 - cosh((x2 - S/2)/ell)/cosh(S/(2ell))],
    sigma11 = (1-nu) sigma0 + nu sigma22,        (epsilon11 = 0 => Delta_sigma11 = nu Delta_sigma22)
    sigma12 = 0.

- Case A, far ahead (x1 -> +inf): u = 0, sigma = sigma0 I (reference state).
- Case A, far behind: width-H relaxed strips (formula with S = H on each period).
- Case B, far ahead: the equilibrium of the pre-existing infinite cracks, i.e. width-2H
  relaxed strips (S = 2H).
- Case B, far behind: the added cracks split each 2H strip into two width-H relaxed strips
  (S = H on each half).

Conditions at every strip edge:
- on every crack face (both Cases, both sides): sigma22 = 0 and sigma12 = 0 (traction free);
  displacement jumps are permitted;
- on the crack line ahead of a tip (symmetry segment): u2 = 0 and sigma12 = 0 (equivalently
  u1 even, u2 odd about the line).

## 2. Energy release G_A(H), G_B(H) — translating-tip energy argument

Strip energy. Direct integration of the DeltaPi density over one period of a relaxed
width-S strip gives (per unit x1 length)

    W(S) = -h sigma0^2 ell tanh(S/(2ell)) / Ebar.

Bookkeeping inside W(S): the linear prestress term contributes
-2 h sigma0^2 ell tanh(S/2ell)/Ebar and the quadratic (strain + foundation) terms together
+1 h sigma0^2 ell tanh(S/2ell)/Ebar — the net is minus one times the quadratic part, i.e.
the relaxed strip sits below the reference state. (SymPy-checked, see VERIFICATION.md.)

Advance accounting. In steady advance all fields translate with the tips, so the energy
released per unit advance of the translating configuration is (W_ahead - W_behind), the
difference of the strip energies of the bounded far equilibria; equivalently this is the
Eshelby integral J = int [w - h sigma_{a1} u_{a,1}] dx2 on any vertical line away from the
tip (the flux term vanishes in the far states, leaving the w-integrals above). Counting:

- Case A: period width H, exactly one advancing tip per period. Ahead W = 0 (reference
  state), behind W = W(H). Per tip: release = -W(H); divide by h (energy per crack area):

    G_A = (sigma0^2 ell/Ebar) tanh(H/(2ell)),
    gamma_A := G_A Ebar/(sigma0^2 ell) = tanh(H/(2ell)).

- Case B: period width 2H, exactly one advancing tip per period (the old cracks are fixed).
  Ahead W = W(2H) (relaxed width-2H state), behind W = 2 W(H) (two width-H strips):

    G_B = (sigma0^2 ell/Ebar) [2 tanh(H/(2ell)) - tanh(H/ell)]
        = (sigma0^2 ell/Ebar) * 2 t^3/(1 + t^2),   t = tanh(H/(2ell)),
    gamma_B = 2 tanh(H/(2ell)) - tanh(H/ell) >= 0.

Modes. Both geometries are symmetric about each crack line (loading equibiaxial, advance
straight ahead), so the local crack-tip fields are pure mode I: KII = 0; mode II is excluded
by the u2-odd symmetry, and there is no mode III in this in-plane model. With
G = KI^2/E as stated in PROBLEM.md:

    KI/(sigma0 sqrt(ell)) = sqrt(E G/(sigma0^2 ell)) = sqrt((1-nu^2) gamma).

## 3. Critical reference stress and limits

Fracture advance at G = Gc gives

    sigma0c = sqrt(Ebar Gc/(ell gamma(H/ell))).

Limits at fixed h, E, nu, k (x = H/(2ell)):

- H/ell -> 0, Case A: gamma_A ~ x = H/(2ell), so sigma0c ~ sqrt(2 Ebar Gc/H) -> inf.
  Mechanical meaning: densely spaced cracks fully relax the film; almost no stored energy
  remains to drive further advance, so the required prestress diverges like H^{-1/2}.
- H/ell -> 0, Case B: gamma_B ~ 2 x^3 = (H/ell)^3/4, so sigma0c ~ 2 ell sqrt(Ebar Gc/H^3)
  -> inf, faster (H^{-3/2}). Meaning: the pre-existing infinite cracks have already relaxed
  the film, and the new cracks only split already-relaxed strips; the extra energy
  difference is cubic in spacing.
- H/ell -> inf, both: gamma -> 1, sigma0c -> sqrt(Ebar Gc/ell). Meaning: widely spaced
  cracks do not interact; the array reduces to the isolated-crack threshold set by the
  foundation length ell.

These asymptotics are properties of this reduced model at fixed foundation stiffness; they
are not, by themselves, evidence about a real substrate (see Section 6).

## 4. Numerical values (nu = 0.3; sampling locations are normalized, not material data)

gamma = G Ebar/(sigma0^2 ell); kI = KI/(sigma0 sqrt(ell)); sigc = sigma0c sqrt(ell/(Ebar Gc))
= gamma^{-1/2}. From numeric_output.txt (values identical in run_all.py assembly):

| ell/H | H/(2ell) | gamma_A | gamma_B | kI_A   | kI_B   | sigc_A | sigc_B |
|-------|----------|---------|---------|--------|--------|--------|--------|
| 0.5   | 1.0000   | 0.761594| 0.559161| 0.832497| 0.713328| 1.145878| 1.337309|
| 1.0   | 0.5000   | 0.462117| 0.162640| 0.648480| 0.384711| 1.471038| 2.479626|
| 1.5   | 0.3333   | 0.321513| 0.060243| 0.540903| 0.234138| 1.763603| 4.074257|
| 2.0   | 0.2500   | 0.244919| 0.027720| 0.472097| 0.158825| 2.020641| 6.006232|

gamma values are nu-independent; the kI values carry the nu = 0.3 factor sqrt(1-nu^2).

## 5. Checks — what each can and cannot detect

Inventory (details, outputs, tolerances and provenance in VERIFICATION.md; single replay
command in REPRODUCE.md):

1. SymPy internal-algebra suite (analytic_checks.py, 8/8 PASS): wake ODE residual exactly 0;
   edge condition sigma22 = 0 at both strip edges; sigma12 = 0; symbolic strip integral
   equals -h sigma0^2 ell tanh(S/2ell)/Ebar (exact simplification plus 4-point numeric
   sampling); linear/quadratic decomposition (-2 : +1); the identity
   2 tanh x - tanh 2x = 2 t^3/(1+t^2). Can detect: algebra mistakes in the wake solution,
   boundary conditions, energy integration and the G_B identity. Cannot detect: a wrong
   model reduction or a wrong translating-tip bookkeeping, since it only manipulates the
   same closed forms.
2. Two-form consistency of gamma_B (numeric_eval.py, |diff| < 1e-15): internal consistency
   only.
3. Independent 2D finite-difference solve with Eshelby J extraction (fd_check.py): solves
   the full 2D membrane-foundation BVP on the symmetry half-domain with ghost-node
   enforcement of the exact edge conditions, on its own grid, and extracts
   G = (2/h)[J(ahead) - J(behind)] from line integrals. Result: |rel err| <= 0.13% on 5
   configurations (declared tolerance 3%), grid-convergence ratio 3.93 (consistent with
   second order), J-behind flux matches -G h/2 to 6e-5 absolute. This route shares with the
   analytics only the model equations and the far-field Dirichlet states (the bounded strip
   equilibria — an input to the model, not the quantity being verified). It tests: the
   translating-tip energy bookkeeping, the period-width and tip-count normalization, the
   h-division, and the claim that the closed-form G coincides with the J of the actual 2D
   field. It cannot validate: the strip-ODE far states themselves (they are prescribed as
   boundary data and the comparison target is built from the same far states), nor anything
   physical about a real film/substrate. It also shares the model PDE with the analytics,
   so a mistake in the model statement itself would pass undetected.
4. Asymptotic-limit ratios (in run_all.py): gamma_A/x -> 1 and gamma_B/(2 x^3) -> 1 as
   x -> 0 (verified to 1e-5/1e-4 at x = 1e-3), gamma_A, gamma_B -> 1 at x = 20 to machine
   precision. Detects errors in the claimed limiting expansions.
5. J path-independence statistics: reported per configuration in VERIFICATION.md, including
   the audit of the relative-spread value 8.07 (a divide-by-near-zero artifact where
   J ~ 0 ahead; measured absolute spread 7.9e-9, i.e. 2.3e-8 of the released energy).

## 6. Model limitations and supported conclusions

- The local (Winkler) foundation and the plane-stress membrane reduction are the given
  model; nothing here calibrates or justifies them for a specific specimen. The foundation
  stiffness k has not been calibrated against any real substrate.
- Conclusions supported by this reduced model: the decay length ell = sqrt(h Ebar/k); the
  closed-form G_A(H), G_B(H) and KI normalizations; pure mode-I symmetry; the critical
  prestress formula and its H/ell -> 0 / -> inf limits, as properties of the model.
- The dense-spacing asymptotics (sigma0c diverging as H^{-1/2} resp. H^{-3/2}) are
  asymptotics of this model at fixed k. They are not automatically a validation of, or a
  prediction for, a real substrate: real compliant substrates have finite thickness,
  distributed (non-local) stiffness, possible debonding and inelasticity, all outside this
  model. No independent physical or expert validation is claimed.
- Steady-state/translating-tip bookkeeping assumes quasistatic advance with the stated
  fixed old cracks (Case B) and equal advance of all tips (Case A); dynamics, crack-path
  stability and interactions with debonding are outside the scope.

## Open questions (recorded, not silently resolved)

- "Decay length" is taken as the slowest eigen-branch (q = 0 longitudinal, ell =
  sqrt(h Ebar/k)); the full wake approach behind a tip also contains strip eigenmodes with
  q = n pi/S and the faster shear branch, all of which decay at least as fast. This matches
  the requested normalization.
- G formulas rest on the steady-state energy/J argument; the FD J route confirms them
  numerically within the same model (Section 5), but no route independent of the model
  equations exists within this task.
- results.json lists normalized quantities and their sampling locations only; no
  dimensional material data were prescribed, and none are invented.
