# GRADE

## Judgments

1. Wake formulation and fields: **right**
   Reason: Decay length ell=sqrt(h Ebar/k); translation-invariant strips with u1=0, u2 hyperbolic; total traction-free faces (u2'=-sigma0/Ebar); A ahead intact prestress, A behind width H; B ahead width 2H, B behind width H; residual sigma11=sigma0(1-nu).
   Evidence: submission/ANSWER.md Strip Fields.
   Check: independent BVP v''-v=0, v'(±W/2)=-1; faces match; U≈-tanh(W/2) within 1.1e-11.

2. Energy release: **right**
   Reason: Translating-tip argument; G=(U_ahead-U_behind)/(N h); A: period H, N=1, U_ahead=0; B: period 2H, N=1, U_ahead=U(2H), U_behind=2U(H); gamma_A=tanh(H/(2ell)), gamma_B=2tanh(H/(2ell))-tanh(H/ell).
   Evidence: ANSWER.md Energy Accounting.
   Check: independent quadrature of stated functional on BVP fields; 16-sample match.

3. Fracture modes and intensities: **right**
   Reason: Reflection symmetry => Mode I, KII=0; stipulated G=(KI^2+KII^2)/E => KI/sqrt=sqrt((1-nu^2) gamma). Intensity inferred from energy, not tip-field extraction (solver states this).
   Evidence: ANSWER.md p.170-176; results.json kI.
   Check: kI reconstructed from independent gamma; all 8 kI samples agree to 1e-6.

4. Propagation thresholds and limits: **right**
   Reason: sigma_c=sqrt(Ebar Gc/(ell gamma)); H/ell->0: GA~sigma0^2 H/(2 Ebar), GB~sigma0^2 H^3/(4 Ebar ell^2); H/ell->inf both G->sigma0^2 ell/Ebar; shielding of B vs A; fixed h,E,nu,k.
   Evidence: ANSWER.md Critical Stress, Spacing Limits.
   Check: independent U for small/large W consistent with those powers (closed-form vs BVP).

5. Requested numerical values: **right**
   Numerical agreement: **16/16** at abs tol 1e-6 (observed diffs ~1e-11 to 1e-12).
   Mapping: ell/H = 0.5,1,1.5,2 as requested; nu=0.3.

6. Verification and reproducibility: **right**
   Solver executed check.py (8 checks, OVERALL PASS); coordinator replay exit_code 0, script_unchanged, numeric_rows_equal_original.
   Independent reviewer check is BVP+Simpson energy, not a rerun of submitted tanh.
   Shared assumption: 1D wake/foundation law; no 2D tip extraction (solver and reviewer).
   Command success is not itself a physical verdict.

7. Physical scope and uncertainty: **right**
   Local foundation, plane stress, imposed equal-advance patterns; no nucleation, substrate 3D nonlocality, plasticity, debonding. Small-spacing divergence not claimed as real-film validation.

## Independent check
File: review_independent_check.py / review_independent_check_results.json
Route: scipy.solve_bvp on v''=v with natural faces; Simpson of h sigma0 eps + (h/2)Ebar eps^2 + (k/2)u^2.
agreement_count: 16/16, tol 1e-6.
closed_form_vs_bvp_max: 1.098e-11.
FD residual on nonuniform mesh is not used as a deciding check.

## Coordinator replay
REPLAY_RECORD.json: timeout 120s python check.py, exit 0; factual replay only.

## Source recognition
Submission states no sources consulted. Formulas match Xia-Hutchinson 2000 periodic-crack reduction (reference used only to check source/model claims). No explicit paper citation in submission.

## Uncertainty / incomplete
No missing requested delivery. Unresolved by design: 2D tip K extraction, experimental k calibration.

## Solver vs reviewer vs replay
Solver: analytic + FE/BVP self-checks of own formulas.
Reviewer: separate BVP energy reduction.
Replay: unchanged command success, not a scientific grade.
