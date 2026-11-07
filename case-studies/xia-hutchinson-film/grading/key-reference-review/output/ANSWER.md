# Independent reference review — periodic film cracks (Xia–Hutchinson 2000)

Checker role. PROBLEM.md energy is authoritative. No replacement answer key.

## Verdict

`reference_status`: **valid_within_stated_model**
`target_quantities_valid`: **true**

The deciding strip equilibria, total-traction face conditions, translating-tip bookkeeping, modes, limits, and sampled (gamma, kI) values are mathematically consistent with PROBLEM.md. Source Eq. 6 typesetting and OCR of Eqs. 28/30 are cosmetic/source-transcription issues, not target-quantity errors. This is not physical validation of a real film.

## Independent derivation (PROBLEM.md)

Wake, u1=0, u2=u(x2): Delta_sigma22 = Ebar u', equilibrium Ebar u'' = (k/h) u, so

    ell = sqrt(Ebar h / k),   Ebar = E/(1-nu^2).

Total traction sigma22=0 on every crack face => u' = -sigma0/Ebar. Midline of a spacing-S cell: u2=0, sigma12=0. Full-period natural form: v'(0)=v'(L)=-1 with L=S/ell, v=Ebar u/(sigma0 ell). Exact half-cell field v=sinh(a-eta)/cosh(a), a=S/(2 ell).

Stated energy density reduces to ell (h sigma0^2/Ebar) [v' + (v'^2+v^2)/2] d eta. Period integral I(S)=-tanh(S/(2 ell)).

- Case A: period H, one moving tip, ahead u=0, behind spacing-H wake. G_A h = -int Delta_Pi_wake => **gamma_A = tanh(H/(2 ell))**.
- Case B: period 2H, one new moving tip; ahead S=2H, behind two S=H cells. **gamma_B = 2 tanh(H/(2 ell))-tanh(H/ell) = 2 G'(H)-G'(2H)**.

Reflection symmetry => KII=0. PROBLEM.md G=(KI^2+KII^2)/E => **kI = sqrt((1-nu^2) gamma)**. Critical prestress sigma_c = sqrt(Gc Ebar/(ell gamma)).

Limits at fixed h,E,nu,k (fixed ell): H/ell->0, gamma_A ~ H/(2 ell), gamma_B ~ (H/ell)^3/4, sigma_c,A ~ sqrt(2 Gc Ebar/H), sigma_c,B ~ 2 sqrt(Gc Ebar) ell/H^{3/2}. H/ell->inf, both gammas -> 1, sigma_c -> sqrt(Gc Ebar/ell) (isolated, XH2000 Eqs. 10–12, 27).

## Source anchors

- Eqs. 11–12: k is calibrated so isolated G_model = sigma0^2 ell/Ebar matches Beuth Eq. 1, giving ell=(pi/2) g(alpha,beta) h. Paper states ell=1.98 h at zero mismatch. Algebra of the calibration is confirmed; g(alpha,beta) tables were not recomputed.
- Eq. 6 printed DW mixes prestress/incremental symbols. PROBLEM.md functional is used; **not treated as a confirmed erratum**.
- Eq. 27: KI/(sigma0 sqrt(ell))=sqrt(1-nu^2)=0.953939201417 at nu=0.3. Recovered from G=KI^2/E and isolated gamma=1.
- Sec. 4.2: identity G=2G'(H)-G'(2H) is unambiguous in paper.txt. OCR of (28)/(30) garbles radicals/prefactors; intended tanh forms match the independent energy integral. Cosmetic, not a target error.
- Eq. 29 / Tada “few percent” claim is a source applicability caveat, not independently recomputed.

## Numeric samples (nu=0.3) — formula vs reference.json

ell/H | gamma_A | gamma_B | kI_A | kI_B
---|---|---|---|---
0.5 | 0.7615941559557649 | 0.5591607318357128 | 0.8324966558009383 | 0.7133276007351032
1.0 | 0.46211715726000974 | 0.16264015856425462 | 0.6484802333969856 | 0.38471098800719444
1.5 | 0.32151273753163434 | 0.060242529715358484 | 0.5409034952316238 | 0.23413821140722893
2.0 | 0.24491866240370913 | 0.027720167547408525 | 0.4720974293378172 | 0.1588249113588349

max |formula − ref| = 0 (bit-identical). Tolerance used: 1e-6 absolute on gamma and kI (author recommendation, far looser than observed).

## Independent numerical check (not the author’s functions)

Full-width P1 FEM of the PROBLEM.md strip functional; natural face tractions v'(0)=v'(L)=-1. Separate lumped-mass energy. Files: `check_independent.py`, `check_independent_output.txt`, `check_independent_results.json`.

- Convergence at L=1: n_el=40…2560, |gamma−tanh(0.5)| = 1.793e-6 … 4.688e-10; err(640)/err(1280)=4.020 (O(h^2)).
- Face slopes at n_el=2560: vp=−0.99991 (target −1).
- Production n_el=4000 vs ref: max |gamma_P1−ref|=2.709e-8; lumped n_el=8000: 1.751e-8.
- Small-H series and H/ell=20 isolated limits confirmed (gA=0.9999999959, gB=0.9999999918; kI_iso=0.9539392014).

Shared assumptions with the reference: 1D mode-I wake, Winkler k, PROBLEM.md energy, total traction, translating-tip far-field conversion, one A-tip per H and one B-tip per 2H. This check does **not** solve a 2D tip field and does **not** validate a real substrate.

## Author code replay (unchanged copy)

Copied `reference_submission/reproduce.py` to `author_replay/` (no edits of reference_submission). Live run exit 0. Outputs: `author_replay/reproduce_output.txt`, `author_replay/reference.json`.

- Symbolic six checks True (matches saved).
- Collocation max |v−exact|=8.549e-15; energy residual 5.551e-16; FD N=80…1280 clean ~4× decay; production mismatch 2.144e-10.
- Controls C1/C2/C3 failed as designed; independent reconstructions match the saved wrong values to printed digits.

Author numerical routes still use the reduced 1D ODE (half-strip essential BCs). They do not import the closed-form v into the energy, except a tanh(a) collocation *initial guess* (non-critical). Convergence and wrong-control audit: consistent with REPRODUCE.md.

## Issues

**Critical:** none for the stated model’s target quantities.

**Non-critical / check gaps:**
1. Eq. 6 and OCR of (28)/(30) are source-transcription noise; problem functional + G=2G'(H)−G'(2H) fix the content without an erratum.
2. All numerical routes (author and checker) are 1D wakes, not 2D KI extraction.
3. Beuth g(alpha,beta) and Tada Eq. 29 “few percent” not recomputed; 1.98 h is a paper statement.
4. Saved `reproduce_output.txt` records a write path `.../reference/reference.json` from the author’s original folder name; not a quantity error.
5. Collocation start uses tanh(a); does not affect the solved field at reported residuals.

**Required reference changes:** none for target quantities. Optional: note that numerical validation is 1D-wake only; keep Eq. 6 as notation, not erratum.

## What the model supports

Supported inside the reduced model: interaction length ell, simultaneous vs sequential shielding, saturation-spacing mechanism of order ell, isolated-crack limit, pure mode I for these two geometries.

Requires more than this model: wavenumber-dependent substrate stiffness, through-thickness/3D fronts, debonding, substrate cracking/yield, T-stress path stability, experimental grading. Mathematical consistency is not physical validation of a real film.
