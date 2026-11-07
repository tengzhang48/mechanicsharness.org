# GRADE.md — independent review (R04)

## Item 1. Wake formulation and fields — **right**
Reason: Decay length ell=sqrt(h Ebar/k); 1D wake ODE, total-traction crack faces, Case A ahead u=0 / behind strip-H, Case B ahead spacing-2H / behind two H-strips; sigma12=0, sigma22=0 on faces.
Evidence: submission/PROVISIONAL_REPORT.md; analytic_checks.py (ODE residual, C, sigma22(S)).
Deciding check: reviewer quadrature of wake energy matches tanh W(S) to ~1e-8 (review_independent_check.py).

## Item 2. Energy release — **right**
Reason: Stated DeltaPi used; translating-tip: G=(W_ahead-W_behind)/(h n_tips); A: period H, 1 tip; B: period 2H, 1 moving tip; G_A, G_B match reference.
Evidence: PROVISIONAL_REPORT.md; run_all.py definitions; replay_output/results.json.
Deciding check: independent (Wa-Wb)/h from PROBLEM.md density; gamma_A, gamma_B agree with tanh.

## Item 3. Fracture modes and intensities — **right**
Reason: Reflection symmetry => mode I, KII=0; G=KI^2/E (PROBLEM.md plane-stress); kI=sqrt((1-nu^2) gamma), nu=0.3.
Evidence: numeric_eval.py; numeric_rows.json.
Deciding check: 8 kI values vs reference.json within 3e-16.

## Item 4. Propagation thresholds and limits — **right**
Reason: sigma0c=sqrt(Ebar Gc/(ell gamma)); H/ell->0: gamma_A~H/(2ell), gamma_B~(H/ell)^3/4; sigma_c,A~sqrt(2 Ebar Gc/H), sigma_c,B~2 ell sqrt(Ebar Gc/H^3); H/ell->inf both gamma->1, isolated threshold sqrt(Ebar Gc/ell). Mechanical shielding interpretation present.
Evidence: PROVISIONAL_REPORT.md; run_all.py critical_stress and asymptotic_limits check (replay).
Deciding check: independent tanh series (review_independent_check.py).

## Item 5. Requested numerical values — **right**
Sixteen samples at ell/H=0.5,1,1.5,2. Agreement **16/16** vs reference at 1e-12 (two 1e-16 kI_B floats).
Evidence: submission/numeric_rows.json; coordinator replay_output/results.json (not original results.json).
Missing: solver never delivered ANSWER.md or results.json; numbers exist in numeric_rows.json. Not repaired.

## Item 6. Verification and reproducibility — **partly**
Solver executed: sympy 8/8 (analytic_output.txt); gamma_B identity; FD J-integral |rel err|<0.13% (fd_output.txt, 3% gate); asymptotic limits. Coordinator: python3 run_all.py exit 0, max |diff| 2.22e-16 vs original table. Reviewer: independent 1D energy quadrature, not formula copy; did not re-execute 2D FD.
Gaps: no REPRODUCE.md/ANSWER.md in submission; STATUS.md still lists FD incomplete while fd_check.py exists; Case A ahead J relative path_spread up to 8.07 (J~0); FD compares to same tanh G_exact (far-field Dirichlet is shared assumption). Successful command is not physical validation.

## Item 7. Physical scope and uncertainty — **partly**
Provisional report states model assumptions and that G rests on J/steady-state pending FD; notes ell vs faster eigenmodes. Little on k-calibration, substrate wavelength, or that small-H limits are not automatically real-substrate validated. No source recognition (Xia/Hutchinson) in submission files.

## Source recognition
None found in submission/.

## Numerical agreement
16/16 (tol 1e-12) for gamma_A, gamma_B, kI_A, kI_B at four ell/H.

## Incomplete obligations
None remaining after this grade file. Solver missing ANSWER.md/results.json recorded as missing delivery, not filled in.
