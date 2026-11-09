# P015 — matched-test problems V1 (antiplane power-law crack) and V2 (bilinear membrane, inverse for T)

Verdict: **V1 READY** (optional: JSON output line, release metadata) · **V2 READY AFTER FIXES** (state the two
input data as given, not as outputs of unnamed calculations; say which strains the lower screen uses; drop an
unneeded stipulation). Both keys are right (recomputed independently, agreement ≤ 4 × 10⁻¹¹). Overall: READY AFTER FIXES.

Source paper: J. W. Hutchinson (1968), "Singular behaviour at the end of a tensile crack in a hardening material",
J. Mech. Phys. Solids 16, 13–31. Text layer and page images: `<runs>/library/P015-build-codex/paper.txt`,
`pages/page-NN.png` (page-NN = journal p. 12+NN). Used: antiplane crack cited, not solved (p. 19; paper.txt:340-346);
piecewise-linear law Eq. 35 (p. 24), plane-stress field Eqs. 39–41 and J relation Eq. 42 (pp. 24–25), plane-strain β
(p. 25), path independence for E_tan > 0 (p. 25), strain ahead of the crack Eq. 48 (p. 27).

Files:
- Problems: `<runs>/library/P015-build-codex/V{1,2}/{TASK_PROBLEM.md, REFERENCE.md, RUBRIC.md, solve_reference.py, reference_output.json}`, `BUILDER_NOTES.md`, `TASK.md` (builder brief), `PROCESS.md`.
- Solver input: header in `<library-pilot-folder>/pilot_entry_format/setup_matched.sh:12-29,49`; knowledge for B0/B1 from `.../library/topics/crack_tip_asymptotics` (setup_matched.sh:30-45). All 41 run `TASK.md` files rebuilt here and byte-identical; knowledge folders byte-identical to the library files (only `__pycache__` added by 8 V2 runs).
- Runs: `<matched-runs>/runs/P015-V{1,2}-*`; logs `.../library-logs/`; grades `<matched-runs>/grading/M*/GRADE.md` (map `pilot_entry_format/LABELS_coordinator_only.json`); grader brief `pilot_entry_format/GRADER_INSTRUCTIONS.md`.
- This audit: `recompute_V1.py` (+`_results.json`), `recompute_V2.py` (+`_results.json`), `scan_V1_branches.py` (+`.json/.out`), `v2_sigma_y_correction.py` (+`.out`), `mutants_V1_V2.py` (+`mutants.json/.out`), `answers_extracted.json`, `check_V1_V2.py` (+`check_V1_V2_audit.json/.out`, `key_from_recompute.json`, `check_key_from_recompute.out`), `V2_TASK_PROBLEM.md.proposed`.

## 1 Statement

**V1 (`V1/TASK_PROBLEM.md`).** Antiplane (mode III) crack in a pure power-law shear material, n = 3.5, J given.
- Complete: geometry, antisymmetry w(x,−y) = −w and the sign (upper face +z), traction-free faces (line 3); law with
  engineering shear defined, no elastic term, W = ∫τ·dγ (lines 5–14); J with contour, normal and orientation (lines
  16–24); leading field with no rigid translation (line 24); every requested quantity with sign basis and face limit
  (lines 28–30); tolerance (line 34). Only the method is open.
- Determined: finite nonzero J fixes w ∝ r^(1/(n+1)); our shooting on the exponent (stress-function form, exponent
  not assumed) finds exactly one singular symmetric eigenvalue, μ = 7/9 (stress ∝ r^(−2/9)), in (0, 1); the other roots
  (μ ≈ 1.287, 1.90) have vanishing stress and J → 0 (`scan_V1_branches.out`). |γ| ∝ r^(λ−1) is monotone on each ray, so
  each contour radius is unique; q(π) ≠ 0, so the face limit exists.
- Paper relation: a variant. The paper only cites the antiplane crack as a simpler problem (p. 19) and gives no
  solution; parameters (τ₀, γ₀, n, J) are problem parameters, not narrative. BUILDER_NOTES.md:9 says the same.
- No invented narrative. Source not named; title says "slit". Cue (minor): line 32 "whether these results **alone**
  determine the load" leans toward "no".
- No grader contested the key: no V1 grade has reference-error or grade-both-ways wording (grep of the 12 GRADE.md).

**V2 (`V2/TASK_PROBLEM.md`).** Plane-stress crack in the paper's piecewise-linear deformation law; infer the tangent
slope T from C_zz, then strain and displacement coefficients and a screened interval.
- Law: lines 7–13 are the paper's Eq. 35 (p. 24) in dimensional form (λ = (3/2)(E/E_tan − 1) ↔ (3/2)(1/T − 1/E)).
  No erratum. Plane stress stated (line 3), constants (line 21), J convention (lines 23–30), C_zz (lines 30–35).
- Determined: the leading field is the linear plane-stress K-field of the asymptotic compliance (paper Eqs. 39–41),
  with E' = T and ν_t = ½ − (½ − ν)T/E. ν_t is the paper's own β of p. 25. J = K²/T is the paper's Eq. 42 (p. 25). J is
  path independent for this law because E_tan > 0 (p. 25).
- Inverse unique (checked): |C_zz|(T) = √(J/2π)(T^(−½) − bT^(½)), b = (1−2ν)/E, decreases strictly for all T > 0,
  from ∞ to 8.09 × 10⁻⁴ mm^½ at T = E. So 0.005 has exactly one root in (0, E). Numerical scan: monotone on (10, E)
  (`recompute_V2_results.json`, "Czz_scan").
- Screens (line 44): defined, endpoints closed. On θ = 0, σ_e = σ exactly and |ε_zz| is the largest principal strain
  (2ν_t = 0.972 > 1 − ν_t = 0.514). r_min = (C_zz/0.05)² = 0.0100 mm does not depend on T; it tests only that ε_zz is used.
- **Statement — invented provenance (lines 23, 30):** "An external calculation supplies J" and "a computation with
  this constitutive law reports the leading coefficient C_zz(0) = −0.005". These values are not in the paper. They
  are the data of an inverse variant, so they are legitimate parameters. Presenting them as results of unnamed
  calculations is invented narrative under the open-book rule. The mathematics is unaffected (0.005 > 8.09 × 10⁻⁴,
  so a root exists). Fix: state them as given (`V2_TASK_PROBLEM.md.proposed`).
- **Statement — lower-screen reading (line 44, "evaluated using the leading field"):** one could instead put the
  leading stresses into the complete law (with the σ_y term). That gives r_min = 0.00927 mm, 7.3% off, which fails
  the 2% band (`mutants.out`). M22 met this reading and rejected it (M22 ANSWER.md:78-80). No graded answer used it,
  and the rubric intends the leading strains (V2/RUBRIC.md:3; REFERENCE.md:101). Fix: one clause.
- Minor — stipulation (line 37): "Matching to the surrounding body is known to make corrections to this leading
  field negligible for r ≤ 0.200 mm" is invented and needed by no answer. The complete law's own σ_y term changes
  ε_zz by −3.9% at r = 0.01 mm, −9.6% at 0.0613 mm and −17% at 0.2 mm (`v2_sigma_y_correction.out`). Line 44 then
  disclaims accuracy, which is consistent with this. Replace it with "All requested quantities refer to this leading
  field."
- Paper relation: close to the paper's own Sec. 6. The forward formula is Eq. 48 (p. 27, plane stress, ε_z ahead of
  the crack) without its incompressible approximation (ν_t → ½). Inverting the ν_t = ½ form of Eq. 48 (ε_z = −σ/E_tan) gives T = 5093 MPa,
  5.7% high, and fails the 1% band (`mutants.out`, "incompressible"). The key's exact compressible form is right
  (recomputed below). REFERENCE.md should say that it departs from Eq. 48 on purpose. The transfer being tested is
  the inversion with the exact ν_t, the polar/displacement coefficients and the screens.
- Cues (minor): line 21 "There is no unspecified hardening exponent" rules out an HRR power-law field; line 46 "including
  the treatment of ε_zz" and line 44 "the three principal strains" point at the two decisive steps (ν_t ≠ ν; ε_zz
  largest). They help; they do not give numbers. M08 still missed the second.
- No grader contested the key (same grep over the 12 V2 GRADE.md).

## 2 Key

Recomputed with our own code, written from the statements. Each route differs from the builder's in the step at
risk. The case page's "our script" (`assets/solution-evidence-2026-10-25/hrr_check.py:192-287`) re-implements the
builder's route A (same displacement ODE and polar J integrand for V1; same closed formulas ν_t, J = K²/T, Kolosov
displacement for V2). It confirms the arithmetic, not the formulation.

- **V1** (`recompute_V1.py`): Airy-type stress function (equilibrium exact) with compatibility as the governing
  equation. Exponent found by shooting: μ = 0.77777777777778, so λ = n(μ−1)+1 = 0.2222222222222 = 1/(n+1). W comes from
  quadrature of ∫τ(γ)dγ, which agrees with nτγ/(n+1) to 5 × 10⁻¹⁶. J is computed on a square contour (h = 1, 0.37) and on
  a circle, in Cartesian components, with w_,x from finite differences of w(x,y). The three contours agree to
  1.5 × 10⁻⁹. The final field gives J = 0.11999999999838 on a third square (h = 0.05).
- **V2** (`recompute_V2.py`): leading compliance computed numerically as lim ε(Lσ)/L of the complete law (mpmath,
  L = 10²⁵). Uniaxial check at T = 5000: slope 5000.0, lateral ratio 0.485714 = ν_t. J is computed on square contours,
  with the rotation integrated from compatibility, so neither J = K²/T nor the Kolosov formula is used: J·T/K² =
  1.00000000006. Then brentq on C_zz(T); C_y by integrating ε_yy along x = 0; screens by root finding on σ_e and on
  the eigenvalues of the 3D strain tensor.

| Item | Key (`reference_output.json`, route A) | Ours | Rel. diff | Origin |
|---|---|---|---|---|
| V1 Δw(0.02) mm | 0.001318126778 | 0.001318126778 | −1.7e−12 | computation (variant) |
| V1 Δw(0.08) mm | 0.001793693738 | 0.001793693738 | −1.7e−12 | computation |
| V1 τ_rz MPa | 38.58659686 | 38.58659686 | −4.2e−13 | computation |
| V1 τ_θz MPa | 138.0106166 | 138.0106166 | −4.5e−13 | computation |
| V1 r(0°), r(90°), r(180°) mm | 0.04689442305, 0.02506612349, 0.01339840659 | same | ≤ 2.4e−12 | computation |
| V2 T MPa | 4816.472187 | 4816.472187 | −1.4e−11 | computation; paper Eqs. 35, 42, β (p. 25) |
| V2 C_rr, C_θθ, C_γrθ mm^½ | 0.003942048238, 0.0006331781351, 0.005731131133 | same | ≤ 5.7e−12 | computation |
| V2 C_y mm^½ | 0.009139018941 | 0.009139018941 | −4.0e−12 | computation |
| V2 r_min, r_max mm | 0.0100000, 0.06132522855 | same | ≤ 3.6e−11 | computation (definitions) |

- The builder's two routes share more than REFERENCE.md suggests. V1 routes A and B use the same amplitude function
  (`V1/solve_reference.py:66-109`); route B's exact I = π(n+1)/(2n) comes from the same integrand. Both V2 routes
  assume K = √(JT) (`V2/solve_reference.py:60,82,85`); J is verified only afterwards, on circles. Our routes remove
  that shared step. Numerically all routes agree.
- Case-page claims (`cases/HRR_Field_Case_Study_2026-10-25.html:1808`). "Builder's two routes": true, partly
  independent (above). "Our script to 3 × 10⁻¹⁰": true (max 2.3 × 10⁻¹⁰ against the 10-digit REFERENCE table,
  `check-outputs/hrr_check.txt:48-66`), but it repeats the builder's formulation. "V1 also by Opus on 7 October":
  recorded only as a log line (`<study-folder>/crack_front/PROCESS_LOG.md:756`, "within 2e-5"); no script or
  output found, so not verified.
- No errata found in REFERENCE.md or RUBRIC.md (RUBRIC tables V1:25-31 and V2:25-31 equal the keys).

## 3 Checker

No automatic checker existed; Codex graders applied `RUBRIC.md`. `check_V1_V2.py` applies the stated tolerances
(V1 all 2%; V2 T 1%, others 2%; relative, signed) to a `{"V1": {...}, "V2": {...}}` file.

(a) Key: 14/14 pass, both from the embedded key and from our recomputed values (`check_key_from_recompute.out`, exit 0).

(b) Graded answers (numbers transcribed from each graded `submission/ANSWER.md`, byte-identical to the run copy;
`answers_extracted.json` gives the lines). The checker agrees with every grader's tolerance verdict, 24/24
(`check_V1_V2_audit.out`). 22 have 7/7. **M20** (V1, C r2, DeepSeek) has 0/7, errors 9.9–52.9%. **M08** (V2, B0 r1, GLM)
has 6/7: r_min 0.02336 against 0.0100. The three complete Kimi K3 V1 answers (ungraded) are 7/7. The fourth K3 V1
folder holds a placeholder draft (`P015-V1-C-r1-k3az/ANSWER.md:1-17`). Case-page leads confirmed:
- 22 of 24 V1/V2 answers have every number within tolerance. With the 4 of 6 probe-round answers from the T1–T3 audit,
  this gives the page's 26 of 30 (`HRR_Field_Case_Study_2026-10-25.html:379`).
- M20's parity error: extending f′ with odd parity cancels the cross term between the half-planes. Our computed mutant
  gives J ratio 0.654321, the grader's ratio 624.136/953.869, and reproduces all seven of M20's printed numbers to their
  rounding (`mutants.out`). The bug is at `P015-V1-C-r2-deepseek/solve.py:44-46`.
- The "GLM V2 screening endpoint" is M08 (labels map). It left the trace term out of the in-plane strains
  (M08 ANSWER.md:34) and so took them as the largest.

(c) Mutants (`mutants.json`; every one fails the "all items pass" rule):

| Mutant | Changed values (key otherwise) | Result |
|---|---|---|
| V1 parity error (M20 class), computed / as reported | Δw 1.833e−3, 2.495e−3; τ 42.40, 151.65; r 0.0717, 0.0383, 0.0205 | 0/7 |
| V1 J: traction term dropped / W term dropped | J₁ × 0.432 / × 0.568 → all seven off | 0/7 |
| V1 complementary work in J / W = τγ/2 | J₁ × 0.691 / × 0.846 | 0/7 |
| V1 law applied to tensor shear (compliance × 2) | Δw 1.538e−3, τ 33.08, 118.3, r 0.0572… | 0/7 |
| V1 contour on tensor shear (abs. γ = 0.02) | r 0.0192, 0.0103, 0.0055 | 4/7 |
| V1 τ_θz sign / τ_rz sign | −138.01 / −38.59 | 6/7 |
| V1 Cartesian components for polar | τ_xz −100.2, τ_yz 102.4 | 5/7 |
| V1 one-face slip | half Δw | 5/7 |
| V1 wrong branch: linear K_III with G = τ₀/γ₀ | Δw 5.53e−4, r 0.00955 (all rays) | 0/7 |
| V2 ν held at 0.30, E → T | T 1833.5 | 0/7 |
| V2 ν_t = ½ (paper Eq. 48 inverted) | T 5093.0; C_θθ, C_rr, C_y, r_max off | 2/7 (**C_γrθ passes**, 1.85%; r_min unchanged) |
| V2 plane-strain J–K relation | T 6176.9 | 4/7 (C_rr, C_γrθ, C_y pass; r_min unchanged) |
| V2 plane-strain field, T right | C_θθ < 0, κ = 3 − 4ν_t, empty interval | 2/7 |
| V2 elastic E in J = K²/E | T 17043 | 1/7 (r_min unchanged) |
| V2 tensor shear / C_γrθ sign / C_rr–C_θθ swapped | 0.002866 / −0.005731 / swapped | 6/7, 6/7, 5/7 |
| V2 κ = 3 − 4ν_t in C_y only | 0.005701 | 6/7 |
| V2 r_min: trace term dropped (M08) / in-plane only / complete law at leading stress | 0.02336 / 0.00279 / 0.00927 | 6/7 each |
| V2 σ_e ≥ σ_y (factor 10 dropped) | r_max 0.200 | 6/7 |

(d) Tolerances. Key precision is ~10⁻¹¹, against bands of 1–2%. V1: τ ∝ J^(1/(n+1)), so the 2% τ bands admit a J error
up to ≈ 9%. Δw ∝ J^(n/(n+1)) (2.6%) and the radii ∝ J (2%) catch it; every J mutant fails Δw and r. V2: T has 1%
(d ln T/d ln |C_zz| = −1.89); r_max ∝ T so its 2% band is weaker than T's; r_min is independent of T and tests only
the ε_zz identification. Single items pass for some wrong models (ν_t = ½: C_γrθ; plane-strain J–K: C_rr, C_γrθ, C_y).
The set rejects all of them, so "all seven within tolerance" is the right pass rule.

(e) Rubric. Numeric rows are decidable and equal the key. Rows not covered by any number: model justification,
verification, initiation/scope (V1/RUBRIC.md:16-19,42,46; V2/RUBRIC.md:16-19,42,46-47) and "interval nonempty"
(V2/RUBRIC.md:33). These are judgement rows.
- **Grade:** the consequential-credit rule (V1/RUBRIC.md:33, V1/RUBRIC.md:19) gives M20 76/100 with 0/7 numbers in
  tolerance. A correct answer with weak justification and no verification could score below that (e.g. 40 + 10 + 5 + 15 = 70). M08's 134% endpoint
  error costs 2 points (98). This is the same weakness as in the T1–T3 audit (fix 5 there). Report numbers pass/fail
  separately from the total, as the case page already does.

## 4 For outside solvers

- Received: `TASK.md` = header (`setup_matched.sh:12-29`) + `V?/TASK_PROBLEM.md`. Conditions B0/B1 add `knowledge/`:
  README, TOPIC.md, checks (crack_tip_checks.py, mutant_detection.txt), hutchinson_1968/ENTRY.md (section 7
  replaced in B0) and reference_data.json; B1 adds reference_code/hrr.py. The knowledge covers the in-plane crack only;
  ENTRY.md:33 names Eq. 42's √(E_tan/E) factor (relevant to V2). No absolute paths in any solver file (grep).
- Output: free-form `ANSWER.md` + `PROCESS.md`. Add the JSON line (V2 proposal; V1 likewise:
  `{"slip_0.02", "slip_0.08", "tau_rz", "tau_tz", "r_0", "r_90", "r_180"}`) so `check_V1_V2.py` can run.
- Budget and tools: 45 min stated; hard timeout 3000 s (`queue_matched.sh:7`); python3 with numpy/scipy; no web.
  All 24 GLM/DeepSeek runs exited 0 with an answer. K3: 3 complete V1 answers; 13 launches stopped by the provider's
  rate limit (`.exit` = 1 and "token rate" in each `.err`; five kept as `.attempt*` folders). These are infrastructure failures, not results.
- Grading: blind Codex with `GRADER_INSTRUCTIONS.md` (= each grade's `TASK.md`), `TASK_PROBLEM.md`, `RUBRIC.md`,
  `REFERENCE.md`, `solve_reference.py`, `reference_output.json`. All copies are byte-identical to the build files. Graders
  may rerun submission code. Record integrity: no run file is newer than its `.exit` time; each graded submission's
  ANSWER.md is byte-identical to the run copy.
- Public: `V1/TASK_PROBLEM.md`; `V2/TASK_PROBLEM.md` after fixes; the header; the knowledge variants (library topic,
  frozen v0.2; not audited here); `GRADER_INSTRUCTIONS.md`; the output schema; scores.
- Hidden: `V{1,2}/REFERENCE.md`, `RUBRIC.md` (key numbers), `solve_reference.py`, `reference_output.json`;
  `BUILDER_NOTES.md` (states ν_t and the decisive steps), build `PROCESS.md`, `scratch/`; `A1/` (out of scope); all
  `grading/M*`; all `runs/P015-V*` (answers hold the keys); `LABELS_coordinator_only.json`; this folder's
  `check_V1_V2.py` (embeds the keys) and all outputs. Never shipped: `paper.txt`, `pages/`.

## 5 Fixes

1. V2/TASK_PROBLEM.md:23 "An external calculation supplies J = 0.800 N/mm. Its convention is" → "The loading gives
   J = 0.800 N/mm, with the convention". Line 30 "a computation with this constitutive law reports the leading
   coefficient" → "the leading coefficient of the transverse strain is given:". Reason: invented provenance (open-book rule).
2. V2/TASK_PROBLEM.md:44 "evaluated using the leading field" → "evaluated with the stresses and strains of the leading
   field (not the complete law at finite r)". Reason: the alternative reading gives r_min = 0.00927 (fails).
3. V2/TASK_PROBLEM.md:37 delete "Matching to the surrounding body is known to make corrections to this leading field
   negligible for r ≤ 0.200 mm." and add "All requested quantities refer to this leading (r^−1/2 strain) field."
   Reason: invented and inconsistent with the law's own 17% correction at 0.2 mm; no answer needs it.
   (Fixes 1–3 and 4 for V2 are in `V2_TASK_PROBLEM.md.proposed`.)
4. V1 and V2 statements: append a JSON output line (fields in section 4). Reason: automatic checking with `check_V1_V2.py`.
5. Release metadata. V1: "Variant: antiplane crack cited on p. 19 of Hutchinson (1968), not solved there; key
   computed (builder's two routes + audit route, ≤ 2.4e−12)". V2: "Adapted from Sec. 6 (Eqs. 35, 40–42, 48; β on
   p. 25): inverse for E_tan with exact compressible ν_t; key departs on purpose from Eq. 48's ν_t = ½
   approximation (5.7% in T)". Add the same Eq. 48 note to V2/REFERENCE.md.
6. Scoring: report "all requested numbers within tolerance" (checker) separately from the rubric total. Reason: M20
   scored 76 with 0/7 (V1/RUBRIC.md:33).
7. Optional, V1/TASK_PROBLEM.md:32 "these results alone determine" → "these results determine". Reason: minor cue.

## 6 Not verified

- Opus's V1 recomputation of 7 October: only a log line (PROCESS_LOG.md:756); no script found.
- Session logs not reread for rule compliance. PROCESS_LOG P256 records FAILURES 58 (P015-V1-B0-r1-glm wrote into
  the B1-r1-glm folder and moved the file back) and judges it uncontaminated; I checked only file identity and times.
- The probe-round half of "26 of 30" rests on the T1–T3 audit, not on a new check here.
- Judgement rows of the rubrics (model, verification, scope) were not re-graded.
- The knowledge files (TOPIC.md, checks, ENTRY.md) were not audited for correctness here.
- V1 branch uniqueness is shown only on μ ∈ (0.02, 2.5) by scan. For V2, the non-uniformity of the leading field at
  the faces (σ_e → 0 at θ = ±π, paper p. 25) is argued, not computed. It does not touch the requested angles or the
  limit J.
- A1 (out of scope).
