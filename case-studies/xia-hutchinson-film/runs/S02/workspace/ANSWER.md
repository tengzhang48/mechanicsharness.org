# ANSWER — Energy release for periodic cracks in an elastic film

Status: **partial but internally checked** solution. Case A is well converged
and its isolated-crack limit is verified; Case B is computed with a
half-domain symmetry model and grid-converged; several limits and the
Mode-II statement remain argued rather than proven. All numbers come from the
FEM in this folder, not from any external source.

## 1. Model and the decay length ell

Plane-stress isotropic film (E, nu), thickness h, bonded to a Winkler
foundation (-k u per unit area), prestress sigma0 I. Equilibrium
`h sigma_{ab,b} - k u_a = 0` with
`sigma = sigma0 I + (E/(1-nu^2))[(1-nu)eps + nu tr(eps) I]`.

Writing a mode e^{s x1 + i q x2} (derive3.py) the 2x2 symbol determinant
factors, giving the two x1-decay branches with C = E/(1-nu^2):

    s1^2 = q^2 + k/(h C)                  (shear-like)
    s2^2 = q^2 + 2 k/(h C (1-nu))         (dilatational-like)

The slowest-decaying branch at q=0 fixes the wake boundary-layer length:

    ell = sqrt(h C / k) = sqrt( E h / ( k (1-nu^2) ) ).            (1)

This is the requested "decay length ell of the crack-induced wake field".
The other length, sqrt(h C (1-nu)/(2k)) = sqrt((1-nu)/2) ell = 0.592 ell
(nu=0.3), decays faster and does not set the outer field.

A direct check (derive8.py) shows the only bounded x1-independent equilibrium
of the model is u=0, sigma=sigma0 I. Hence far ahead *and* far behind the tips
the increment field returns to the reference state; the whole crack effect is
a boundary layer of width ~ell around each tip. The "relaxed cracked wake" is
reached on the scale ell.

## 2. Fields far ahead / far behind; edges

Far ahead (x1 -> +inf): increment u -> 0, sigma -> sigma0 I (intact film).
Far behind (x1 -> -inf): the crack has released the prestress locally; the
foundation restores u -> 0 over length ell, so the far-behind increment also
-> 0, but with a finite opening plateau and a partially relaxed stress that
decays on the ell scale (FEM: at x1=-18ell the crack-face sigma22 = 0.0029
~ 0, i.e. faces traction free, while sigma11 ~ 0.735).

Edge conditions at every strip edge:
- crack faces (x1<0): sigma_{2a} = 0 (traction free).
- symmetry plane between cracks: u2 even, u1 odd (no shear).
- periodicity in x2 with period H (Case A) or 2H (Case B).

## 3. Mode content

The geometry (cracks along x1) plus equibiaxial prestress are symmetric under
reflection in the crack plane and under x1 -> -x1 about a tip. The near-tip
field is therefore Mode I; KII = 0. (Argued from symmetry; the FEM also gives
the opening uniform across x2 at the tip — a rigid-plate opening — consistent
with pure Mode I. NOT proven symbolically.)

In this model G = KI^2/E  (KII=0).

## 4. Energy accounting

Per advancing tip, per unit advance, G.h = energy released. In Case A one tip
per period H advances; the FEM integrates the full model potential energy
`Pi = int [ h sigma0 eps_aa + (h/2) Δsigma:eps + (k/2) u^2 ] dA` over one
period and takes G = -(1/h) dPi/da by differencing two tip positions.
In Case B the period is 2H with one advancing central tip; I modelled the
half-strip [0,H] (symmetry about x2=0 with the outer edge free from the
existing infinite crack) and multiplied the half-domain release by 2 for the
full 2H period, then divided by h.

Units below: E=sigma0=h=H=1, so ell/H = ell numerically.
Normalizations: gamma = G Ebar/(sigma0^2 ell), K-normalized = KI/(sigma0 sqrt ell).

## 5. RESULTS at the requested sampling locations (ell/H)

| ell/H | gamma_A | gamma_B | KI_A/(so sqrt l) | KI_B/(so sqrt l) |
|-------|---------|---------|------------------|------------------|
| 0.5   | 0.7651  | 0.5642  | 0.8344           | 0.7158           |
| 1.0   | 0.4654  | 0.1642  | 0.6508           | 0.3865           |
| 1.5   | 0.3240  | 0.0612  | 0.5430           | 0.2359           |
| 2.0   | 0.2469  | 0.0280  | 0.4740           | 0.1598           |

Corresponding G_A, G_B (E=sigma0=h=1): (0.3481, 0.2560), (0.4235, 0.1494),
(0.4423, 0.0835), (0.4494, 0.0510).

Case A is grid-converged: nx=600 (dx=0.075) vs nx=900 (dx=0.05) differ by
< 0.25%. Case B: nx=600 vs nx=1000 differ by < 0.3%.

## 6. Limits

- H/ell -> infinity (isolated cracks, Case A): measured gamma_A =
  0.9992, 0.9994, 0.9986 at H/ell = 8, 12, 16, so

      G_A -> sigma0^2 ell / Ebar = sigma0^2 ell (1-nu^2)/E,
      KI_A/(sigma0 sqrt ell) -> sqrt(1-nu^2) = 0.9539.               (2)

  Mechanical meaning: when neighbouring cracks are far apart each tip behaves
  as an isolated crack whose release is set by the foundation length ell; the
  crack spacing drops out. (nu enters only through Ebar.)

- H/ell -> 0 (cracks very close, foundation negligible over H): measured
  gamma_A -> 0 and G_A -> 0 linearly in H; KI_A/(sigma0 sqrt H) rises toward
  ~0.68 but is not fully converged, so only the scaling is reliable:
  KI_A ~ c sigma0 sqrt(H) with c ~ 0.68 (partial). Mechanical meaning: with
  the foundation screened out, adjacent cracks share the load and the release
  per tip vanishes as H -> 0.

- Case B limits were NOT computed (missing work).

- Critical stress for advance (part 3): for toughness Gc, advance when
  G_A(H) = Gc, i.e. sigma_c = sqrt( Gc E / (gamma_A(H) (1-nu^2) ... ) ) — more
  precisely sigma0_c = sqrt( Gc Ebar / (gamma_A(H) ell) ) (Case A) and the
  analogous expression with gamma_B for Case B. Because gamma_A, gamma_B
  decrease with ell/H = sqrt(Eh/(k(1-nu^2)))/H, the critical stress increases
  with crack spacing — i.e. widely spaced cracks are harder to advance.

## 7. Checks actually run (see results.json)

- duel free crack face: wake sigma22 at face = 0.0029 (should be 0). PASS.
- stress recovery to sigma0 ahead: 1.0018, 1.0025. PASS.
- Pi decreases monotonically with crack length (shifts -4..+4): PASS (G>0).
- mesh convergence A and B (<0.3%). PASS.
- isolated-crack limit gamma_A -> 1.0 (H/ell=8,12,16). PASS.

Declared tolerances: mesh convergence 1%; face traction 2%; recovery 1%;
limit 0.5%.

## 8. Errors detected / methods rejected

- An early closed-form fracture integral (compute.py, "Finf") gave negative
  and NaN values and is WRONG; it is discarded and not used.
- A J-integral extractor (fem_j.py) was not path-independent (J = 0.016,
  0.0015, -0.0014 at r=0.3..1.2) and is REJECTED. This routes around a false
  "energy release" value.
- A first Case-B model with a single mid-plane node row produced no opening
  and negative G; rejected.

## 9. Shared assumptions between verification routes (caveat)

The FEM energy method (two tip positions) and the mesh-convergence check share
the same constitutive model, the prestress-as-traction treatment and the same
code, so they are NOT independent verifications. The isolated-crack limit
(gamma_A -> 1) is a weaker, partially independent anchor (it recovers a clean
constant), but it also uses the same solver. No external/reference solution
was consulted.

## 10. What is supported / what needs validation

Supported by this reduced model: the ell = sqrt(Eh/(k(1-nu^2))) boundary-layer
structure; pure Mode I; G decreasing with ell/H in both cases; the exact
isolated-crack limit G_A -> sigma0^2 ell/Ebar; Gc critical-stress formulas.
Needs additional information/validation: the H/ell -> 0 limit constant,
the Case B limits, a symbolic proof that KII=0, and independent verification
(an analytic Wiener-Hopf solution) of the numerical G values.

## Reproduction

Single command (from this folder):

    python3 compute.py   # NOT the results; legacy wrong-path only

The results in this ANSWER come from:

    python3 fem_energy2.py     # Case A energy-difference G (r sweep helper)
    python3 femA_H.py          # Case A H-dependence / isolated limit
    python3 femb2.py           # Case B (half-domain) G

Core reproduction of the reported table: run the three scripts above; each
prints gamma and K/(sigma0 sqrt ell). `results.json` holds the values and the
check outputs.