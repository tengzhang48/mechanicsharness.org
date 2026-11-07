# GRADE — independent review (R01)

Reviewer checks: `independent_check.py` (1D ODE energy, not solver FEM, not `reference/reproduce.py`).
Coordinator replay: `REPLAY_RECORD.json` (factual only). Solver verification: `submission/results.json` checks.

## Item 1. Wake formulation and fields — **wrong**

Reason: ell = sqrt(Ebar h/k) is correct, but the solver concluded (derive8.py, ANSWER §1–2) that the only bounded x1-independent equilibrium is u=0, so far-ahead and far-behind both return to the reference. That is false: cracked-strip wakes with u1=0, u2=u(x2), Ebar u''=(k/h)u, and total traction sigma22=0 on faces admit the bounded field v=sinh(a-eta)/cosh(a). Far ahead of Case A is intact (u=0); far behind is a nonzero 1D wake. Case B ahead is the spacing-2H wake, behind two spacing-H cells. “Total crack-face traction” is applied as prestress cancellation in FEM (good), but the analytic far-field states are incorrect.

Evidence: submission/ANSWER.md L30–42; submission/derive8.py L20–26.
Deciding check: independent tanh/FD wake (independent_check_output.json); fd_err_A ≤ 7e-8.

## Item 2. Energy release — **partly**

Reason: Period/tip/h bookkeeping is stated (A: one tip per H; B: one new tip per 2H; G = −(1/h) dPi/da). The stated Pi functional is used. Translating-tip reduction to wake energy difference is not derived; G is 2D FEM differencing. Isolated gamma_A→1 is recovered. Closed forms gamma_A=tanh(H/(2 ell)), gamma_B=2tanh(H/(2 ell))−tanh(H/ell) are absent. Case B half-domain ×2 is hand-applied (solver unresolved).

Evidence: ANSWER.md §4–5; femb2.py L96–98.
Deciding check: reviewer gamma vs submitted; diffs ~0.002–0.005, not the exact identities.

## Item 3. Fracture modes and intensities — **partly**

Reason: Reflection symmetry ⇒ KII=0 and G=KI^2/E (PROBLEM plane-stress) is justified in outline. kI=sqrt((1-nu^2) gamma) is used numerically as sqrt(G E)/… with E=1. Isolated kI→sqrt(1-nu^2)=0.9539 is correct. Solver marks KII as unproven. FEM KI samples miss 1e-3 vs exact except one B point.

Evidence: ANSWER.md §3; results.json kI_A/kI_B.
Deciding check: kI_exact = sqrt((1-nu^2)*gamma_exact); 0/16 at 1e-6.

## Item 4. Propagation thresholds and limits — **partly**

Reason: sigma_c = sqrt(Gc Ebar/(gamma ell)) is given. H/ell→∞: gamma_A→1 and isolated mechanical meaning correct; gamma_B→1 not computed (missing). H/ell→0: gamma_A→0 linearly is the right scaling (gamma_A~H/(2ell)); they did not state gamma_B~(H/ell)^3/4 nor sequential shielding. Critical-stress vs spacing: they say widely spaced cracks are harder; isolated sigma_c is actually the saturation floor while dense A-cracks have lower sigma_c,A~sqrt(2 Gc Ebar/H). Mixed.

Evidence: ANSWER.md §6; results.json limits.B_limits="not computed".
Deciding check: series at ell/H=2: gamma_A~0.25 vs exact 0.24492; gamma_B~0.03125 vs 0.02772.

## Item 5. Requested numerical values — **partly**

Sixteen submitted samples exist with requested names/locations. None agree with the independent exact values at 1e-6. At 1e-3, 3/16 agree (gamma_B 1.5, gamma_B 2.0, kI_B 2.0). Trend and isolated limit are close; FEM bias ~0.3–0.9% relative. Coordinator replay: 11/16 match submitted JSON to display rounding; 5 Case-B table entries do not.

Agreement count: **0 / 16** at 1e-6 (scientific vs exact 1D reduction). Display-rounding solver-vs-replay: 11/16 (not a scientific verdict).

## Item 6. Verification and reproducibility — **partly**

Solver checks (mesh, face traction, Pi monotonic, isolated gamma_A) ran internally; they share the same FEM and cannot validate the wake reduction. Solver correctly rejected compute.py Finf and non-path-independent J. Declared command in ANSWER.md is split/wrong (`python3 compute.py` marked NOT the results). Coordinator ran `python3 reproduce_all.py`, exit 0, with unguarded fem_energy2 negative G / nan KI on import. Reviewer independent FD confirms tanh identities, not the FEM table.

Evidence: REPLAY_RECORD.json; REPLAY_STDERR.txt; submission/results.json checks; independent_check_output.json.

## Item 7. Physical scope and uncertainty — **partly**

Solver notes the model is Winkler/plane-stress, needs Wiener–Hopf/analytic check, and that H→0 constant and Case B limits are open. They overclaim that far-field u→0 is the physical wake. They do not distinguish k-calibration (Beuth/isolated G) from real-substrate wavelength dependence, 3D fronts, or debonding. Small-spacing limits are not treated as automatically real-film validated.

Source recognition: submission does not cite Xia–Hutchinson 2000; numbers are FEM, not the paper closed form. Reviewer used PROBLEM.md + independent ODE; paper.txt only to confirm source-claim check is unused by the solver.

## Independent-check details

Command: `python3 independent_check.py`
FD vs tanh: max |Δgamma| ≈ 7.0e-8.
Submitted vs exact (1e-6): 0/16; (1e-3): 3/16.
Tolerance used for “agreement count out of sixteen”: 1e-6 absolute on gamma and kI (reference recommendation).

Unresolved: no 2D tip-field KI extraction by reviewer; Case B FEM energy dump not re-integrated; solver 2D vs 1D residual not mesh-extrapolated to the tanh limit.

Uncertainty: FEM table is a coarse approximation of the correct 1D-wake G, not a validation of a real film.
