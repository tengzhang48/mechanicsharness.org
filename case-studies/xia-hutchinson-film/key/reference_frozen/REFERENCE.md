# REFERENCE.md — Energy release for periodic cracks in an elastic film (reference derivation)

Author role: independent reference author (per TASK.md). Sources: PROBLEM.md (authoritative model),
paper.txt/paper.pdf = Z.C. Xia & J.W. Hutchinson, "Crack patterns in thin films", J. Mech. Phys.
Solids 48 (2000) 1107–1131 (XH2000), especially Section 4.2, Eqs. 28 and 30; Eqs. 4–12 for model
and calibration; Eq. 27 isolated-crack limit. Scope limited to the two array geometries of PROBLEM.md
(Case A simultaneous array; Case B sequential midway cracking). No finite/curved/spiral-crack content.

**Completion status:** Sections 1–7 are complete derivations from the problem's stated functional and
boundary conditions; published formulas used as source checks. Section 8 (independent numerical strip
BVP + integrated-energy validation, controls, actual outputs) is implemented in `reproduce.py`;
its actual outputs are transcribed in Section 8 and REPRODUCE.md, and the numbers are in
`reference.json` (run completed 2026-10-15T12:23Z, exit 0). STATUS.md tracks minute-level progress.

---

## 0. Model, notation, and which variables are fixed

Film: thickness h, Young's modulus E, Poisson ratio nu (plane stress), on a Winkler foundation
f_alpha = −k u_alpha, k [Pa/m]. Reference state: uniform equibiaxial tensile prestress
sigma0·delta_alpha_beta. Ē ≡ E/(1−nu²). Incremental fields relative to the prestressed state:
u_alpha, epsilon_alpha_beta, Delta_sigma_alpha_beta = Ē[(1−nu)eps + nu·eps_gg·delta].
Equilibrium: h·sigma_alpha_beta,beta − k·u_alpha = 0 with **total** stress sigma = sigma0·delta + Delta_sigma.
Crack faces: total traction free, sigma_alpha_beta·n_beta = 0, i.e. Delta_sigma·n = −sigma0·n.

Stated incremental potential energy of a film region Omega (PROBLEM.md; this resolves XH2000 Eq. (6),
whose printed first and second terms are notationally ambiguous in the source — see Section 7):

    Delta_Pi(Omega) = ∫_Omega [ h·sigma0·eps_alpha_alpha + (h/2)·Delta_sigma_alpha_beta·eps_alpha_beta
                                 + (k/2)·u_alpha·u_alpha ] dA .

**Fixed during any increment of crack advance** (PROBLEM.md): h, E, nu, k, sigma0 (the reference
prestress and model parameters are frozen); the increment is purely geometric (tip advance).
Numbers below use nu = 0.3; the sampling values ell/H in {0.5, 1, 1.5, 2} are normalized sampling
locations, not dimensional material data. G is energy per advancing crack **area** (J/m²): release
per unit in-plane tip advance per tip divided by h. Mode relation in this model:
G = (KI² + KII²)/E.

## 1. Decay length; strip (wake) boundary-value problem and its exact solution

Reduce the field to a single transverse coordinate far from the tips, where u is x1-independent.
For mode-I symmetry u1 = 0, u2 = u(x2): eps22 = u', eps11 = 0, so
Delta_sigma22 = Ē·u', Delta_sigma11 = nu·Ē·u', Delta_sigma12 = 0, and equilibrium gives

    Ē·u''(x2) = (k/h)·u(x2)  ⟺  u'' = u/ell²,   **ell ≡ sqrt(Ē h / k)** .

ell is the interaction (decay) length of the crack-induced wake field; it is the length XH2000 call l
(Eq. 12), where the substrate spring is calibrated so that ell = (pi/2)·g(alpha,beta)·h
(≈ 1.98 h with no elastic mismatch, since (pi/2)·g(0,0) = 1.976 from Eq. (1) with Beuth's zero-mismatch
value; see CORRECTION_RECORD.md — an earlier draft mistranscribed the prefactor as sqrt(pi/2), an OCR
radical-scope error; the PDF glyph geometry places the radical only over Ebar·h/k on the left side).
All interaction effects enter only through H/ell.

Periodic wake with spacing S (cracks at x2 = mS). On the half-strip between a crack face (x2 = 0)
and the midline to the next crack (x2 = S/2), with dimensionless v = u/(sigma0·ell/Ē) and
eta = x2/ell, a = S/(2 ell):

    v'' = v on (0, a),   v'(0) = −1  (crack face: total traction free, Ē·u'(0) = −sigma0),
                          v(a) = 0   (midline: u2 = 0 by antisymmetry, and sigma12 = 0 there).

**Conditions satisfied at every strip edge** (total-traction form):
- On every crack face x2 = mS: sigma22 = sigma0 + Delta_sigma22 = 0 (i.e. Delta_sigma22 = −sigma0)
  and sigma12 = 0.
- On every midline x2 = (m+1/2)S: u2 = 0 (opening antisymmetry) and sigma12 = 0.
Far ahead of the tips: Case A reaches the bounded equilibrium u ≡ 0 with decay to the prestressed
state (the unique bounded solution); Case B reaches the spacing-2H wake below.

Exact strip solution (used as the closed-form check; the numerics never assume it):

    v(eta) = sinh(a − eta)/cosh(a),  i.e.  u(x2) = −(sigma0·ell/Ē)·sinh((x2 − S/2)/ell)/cosh(S/(2 ell)) .

Opening at the face: u(0) = (sigma0·ell/Ē)·tanh(S/(2 ell)); stress changes
Delta_sigma22(x2) = −sigma0·cosh((x2 − S/2)/ell)/cosh(S/(2 ell)), Delta_sigma11 = nu·Delta_sigma22,
Delta_sigma12 = 0. Isolated crack S → ∞ recovers XH2000 Eqs. (8)–(9) (u = ±sigma0·sqrt(h/(Ē k))·e^{−|x2|/ell}).

## 2. Energy state bookkeeping and G_A (Case A: simultaneous array, spacing H)

Advance all tips of Case A by unit d x1 in steady translation: the field translates with the tips, so
the swept strip of width H per tip (one moving tip per period H) is converted from the far-ahead state
(u ≡ 0, Delta_Pi = 0) to the far-behind spacing-H wake. Hence, per unit tip advance,

    G_A·h = −∫_0^H Delta_Pi_wake dx2 .

With the exact wake field, integrate the stated functional term by term over one period (SymPy-verified
in reproduce.py): the first-order prestress term gives −2h·sigma0·u(0) = −2h·sigma0²·ell·tanh(H/(2 ell))/Ē;
the elastic + foundation terms give +h·sigma0²·ell·tanh(H/(2 ell))/Ē exactly (their integrand difference
(hĒ/2)u'² − (k/2)u² is the first integral constant C = (h sigma0²/2Ē)·sech²(a)). Therefore

    ∫_0^H Delta_Pi_wake dx2 = −h·sigma0²·ell·tanh(H/(2 ell))/Ē,
    **G_A(H) = (sigma0²·ell/Ē)·tanh(H/(2 ell))** .

Source check: this is exactly XH2000 Eq. (28) first form. Normalized:
gamma_A ≡ G_A·Ē/(sigma0²·ell) = tanh(H/(2 ell)) (dimensionless).

Normalization bookkeeping: period H, one tip per period, thickness division by h; G [J/m²] =
[h·J/m per period per unit advance]/h.

## 3. Case B (sequential midway cracking): ahead/behind states, G_B

Case B: pre-existing infinite cracks at x2 = 2mH (spacing 2H) are fully relaxed far ahead; the new
semi-infinite cracks at x2 = (2m+1)H advance with tips at x1 = 0. Far ahead state = wake with S = 2H;
far behind state = wake with S = H. Only the new tips move: one moving tip per period 2H. Writing
Pi_bar(S) = per-area incremental energy of the spacing-S wake = −(h·sigma0²·ell/Ē)·tanh(S/(2 ell))/S,

    G_B·h = (Pi_bar(2H) − Pi_bar(H))·2H  ⟹  **G_B(H) = (sigma0²·ell/Ē)·[2 tanh(H/(2 ell)) − tanh(H/ell)]**
                                            = 2·G′(H) − G′(2H),   G′(S) ≡ (sigma0²·ell/Ē)·tanh(S/(2 ell)).

Source check: XH2000 §4.2 derive exactly this via G = 2G′(H) − G′(2H) and state Eq. (30) follows
from Eq. (28); the equivalent printed form is (2·sigma0²·ell/Ē)·[tanh(H/(2 ell)) − (1/2)tanh(H/ell)].
(The paper.txt OCR of Eqs. (28)/(30) garbles radical scopes/prefactors; the identification above is
anchored in the paper's own stated identity G = 2G′(H) − G′(2H), which my independent integration
reproduces.) Normalized: **gamma_B = 2 tanh(H/(2 ell)) − tanh(H/ell)**, and 0 < gamma_B < gamma_A for
all H > 0: sequential cracking is strictly shielded relative to simultaneous appearance (XH2000 Fig. 3).

## 4. Modes permitted; stress intensities

Both geometries are reflection-symmetric about every crack line (Case B: the new line x2 = (2m+1)H
bisects old cracks at 2mH and (2m+2)H), and the prestress is equi-biaxial with no shear on any line.
Hence KII ≡ 0 at every advancing tip: **pure mode I** is the only mode the symmetries permit. With
G = KI²/E and G·Ē/sigma0²·ell = gamma:

    **kI_A ≡ KI_A/(sigma0·sqrt(ell)) = sqrt((1 − nu²)·gamma_A) = sqrt((1−nu²)·tanh(H/(2 ell)))**,
    **kI_B ≡ KI_B/(sigma0·sqrt(ell)) = sqrt((1 − nu²)·gamma_B)** .

Both tend to sqrt(1 − nu²) = 0.9539… for nu = 0.3 as H/ell → ∞, the isolated semi-infinite crack value
(XH2000 Eq. 27; exact within the model by the calibration Eq. 12). KI units Pa·sqrt(m).

## 5. Critical prestress and the two limits (fixed h, E, nu, k ⟹ fixed ell)

G = Gc (positive mode-I film toughness) gives the critical reference stress for advance:

    sigma_c,A(H) = sqrt( Gc·Ē / (ell·tanh(H/(2 ell))) ),
    sigma_c,B(H) = sqrt( Gc·Ē / (ell·[2 tanh(H/(2 ell)) − tanh(H/ell)]) ).

- **H/ell → 0** (dense array, fixed ell): gamma_A ≈ (H/2ell)(1 − (H/ell)²/12 + …) → sigma_c,A ≈
  sqrt(2 Gc Ē / H)·(1 + O((H/ell)²)); gamma_B ≈ (H/ell)³/4 + O((H/ell)⁵) → sigma_c,B ≈
  2·sqrt(Gc Ē)·ell/H^{3/2}·(1 + O((H/ell)²)). Meaning: closely spaced cracks shield each other; the
  midway set finds the film nearly fully relaxed ahead of it, its driving force being cubic-small, so
  sequential cracking arrests at a finite saturation spacing of order ell — the model's mechanism for
  a crack-spacing floor (and sigma_c,B ≫ sigma_c,A as H → 0).
- **H/ell → ∞**: tanh → 1, gamma_A, gamma_B → 1, sigma_c → sqrt(Gc·Ē/ell) in both cases: the isolated
  semi-infinite crack result (XH2000 Eqs. 10–12); interaction is negligible beyond H ≈ 3 ell (XH2000
  §4.2), which is the model's quantitative interaction range.

## 6. Sampling values (nu = 0.3, 1 − nu² = 0.91) — analytic, numerically validated to ≤ 2.2e−10

ell/H | H/ell | gamma_A | gamma_B | kI_A | kI_B
0.5   | 2     | 0.7615941559557649 | 0.5591607318357128 | 0.8324966558009383 | 0.7133276007351032
1     | 1     | 0.4621171572600097 | 0.1626401585642546 | 0.6484802333969856 | 0.3847109880071944
1.5   | 2/3   | 0.3215127375316343 | 0.0602425297153585 | 0.5409034952316238 | 0.2341382114072289
2     | 0.5   | 0.2449186624037091 | 0.0277201675474085 | 0.4720974293378172 | 0.1588249113588349

(Reference.json carries full precision, together with the numerically validated values.)

## 7. Physical scope of the approximation; source-notation note

- The film/substrate system is reduced to a plane-stress sheet on a Winkler foundation with a single
  spring constant k. The real substrate is a half-space whose effective stiffness depends on
  wavenumber; k is calibrated (XH2000 Eq. 11) to reproduce the isolated-crack energy release
  (Eq. 1/Beuth 1992). Accuracy for arrays is correspondingly limited: XH2000 report, via Eq. (29) and
  Tada et al. plane-strain edge-crack solutions, agreement within a few percent over the whole H/ell
  range for zero mismatch. Through-thickness variations, finite film-thickness (3D) effects at scales
  ~h, interface debonding, substrate cracking/yielding are excluded; kink/curvature stability (T-stress)
  is outside the two array solutions here.
- The channeling energy argument (transfer of far-ahead strip to far-behind wake per unit advance) is
  exact within the stated model, since the steady solution merely translates with the tips; the paper's
  remark that G = 2G′(H) − G′(2H) is exact, with any error from the model approximation of the
  substrate, is consistent with this and reproduced above.
- XH2000 Eq. (6) prints the incremental energy density in a form whose first two terms are notationally
  ambiguous in the source rendering. PROBLEM.md supplies the consistent incremental functional
  (sigma0·eps first order, (1/2)Delta_sigma·eps second order, (k/2)u² foundation); this reference uses
  the problem's functional and boundary conditions as authoritative and does **not** assert any
  erratum. Similarly, paper.txt's OCR of Eqs. (28)–(30) obscures radical scopes/prefactors; my
  independent integration reproduces Eq. (28)'s tanh form and the paper's stated identity
  G_B = 2G′(H) − G′(2H), which fixes the intended content without any erratum claim.

## 8. Independent numerical validation (strip BVP + integrated stated energy) and controls

Implemented in `reproduce.py` (single command `python3 reproduce.py`); actual outputs in REPRODUCE.md
and reproduce_output.txt; numbers in reference.json. All routes share: plane-stress strip reduction with
u1 = 0 (mode-I symmetry), same fixed parameters, total-traction face condition and midline symmetry, same
period/tip/thickness normalization. The numerical routes never use the closed-form strip solution.

- Symbolic (SymPy): exact strip field satisfies v''=v, v'(0)=−1, v(a)=0; dimensional equilibrium
  h·Ē·u_{,x2x2} = k u with k = Ēh/ell²; the dimensional functional reduces to
  ell·(h sigma0²/Ē)[v' + (v'²+v²)/2]deta; period integral I = ∫_0^{2a}[...]deta = −tanh(a).
  All six checks returned True.
- Collocation strip BVP (scipy.solve_bvp): max|v_num − sinh(a−eta)/cosh(a)| = 8.5e−15 over
  S/ell ∈ {0.5, 2/3, 1, 2, 3, 4}; Simpson-integrated energy residual |I_num + tanh(a)| ≤ 5.6e−16.
- Independent 2nd-order finite-difference BVP: gamma error decays as O(N^{−2}) (N = 80…1280:
  7.6e−7 … 3.1e−9 at S/ell = 0.5; 1.1e−4 … 4.5e−7 at S/ell = 4), confirming convergence; the
  production FD route (N = 40000) plus collocation reproduce every tabulated gamma_A, gamma_B to
  within 2.2e−10 of the analytic values (asserted < 5e−9).
- Recommended comparison tolerances (from observed accuracy): gamma, kI to 1e−6 absolute; dimensional
  G to 1e−6 relative — noting that physical input uncertainty (substrate stiffness wavenumber
  dependence, k calibration) dominates in practice.

### Controls (deliberately wrong; actual outputs preserved in reproduce_output.txt; NOT production answers)
- C1 omit foundation energy ((k/2)u² dropped from the stated functional): wrong gamma_A =
  (0.932404, 0.496564, 0.332831, 0.249876) at ell/H = (0.5, 1, 1.5, 2) vs correct
  (0.761594, 0.462117, 0.321513, 0.244919); max deviation 1.7e−1 → detected (the wake would even
  mis-order states at smaller H/ell).
- C2 incremental traction Delta_sigma22 = 0 on faces instead of total traction: BVP collapses to
  v ≡ 0 (max|v| = 0.0e+00), gamma_A = −0.0 vs correct 0.462117 → detected (no cracking release at all).
- C3 initial/final strip spacings confused (states swapped: 2G′(2H) − G′(H)): wrong gamma_B =
  (1.166461, 1.061071, 0.844053, 0.679316) vs correct (0.559161, 0.162640, 0.060243, 0.027720),
  deviations up to 24× → detected.

## 9. Unresolved / out-of-scope items

- Wavenumber-dependent substrate stiffness, 3D channeling front details, crack-face friction,
  debonding/substrate fracture, and path stability (T-stress) are outside this reduced model; no
  experimental grading or physical validation of a real film is performed or implied.
- The paper's Eq. (6)/Eq. (28)/(30) printed ambiguities are noted (Section 7) without erratum claims.
