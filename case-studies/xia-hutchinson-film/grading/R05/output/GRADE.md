# GRADE.md

Independent scientific review. UTC 2026-10-16T00:07:56Z. No desired verdict supplied.
submission/ and reference/ were not edited.

Numerical agreement: **16 / 16** at 1e-12 (reviewer tanh/sqrt formulas vs submission/results.json).
Independent P1 FEM vs the same 16 values: **16 / 16** at 1e-6 (max |d gamma| = 2.7e-8 at n_el=4000).

| # | Item | Judgment | Reason | Evidence | Deciding check |
|---|------|----------|--------|----------|----------------|
| 1 | Wake formulation and fields | right | ell=sqrt(h Ebar/k); A ahead u=0; A behind sinh strip with total-traction faces; B ahead S=2H, behind two S=H cells. | submission/ANSWER.md §1–2; compute.py | replay u2'(0)=u2'(H)=-sigma0/Ebar; FEM natural slopes -> -1 |
| 2 | Energy release | right | GA, GB, one moving tip per period (H vs 2H), divide by h. ANSWER.md L88–90 wrongly says int h sigma0 eps_gg vanishes; compute.py includes it and the tanh form is correct. | ANSWER.md §4; compute.py W_dens | FEM I[v]=v(L)-v(0)+int (v'^2+v^2)/2; gamma_A=-I(H); gamma_B=I(2H)-2I(H); int v' != 0 |
| 3 | Fracture modes and intensities | right | Reflection: KII=0; G=KI^2/E; kI=sqrt((1-nu^2) gamma). | ANSWER.md §3,5 | kI vs sqrt(0.91*gamma), 8/8 |
| 4 | Thresholds and limits | right | sigma_c=sqrt(Gc Ebar/(ell gamma)); H/ell->0: gamma_A~H/(2ell), gamma_B~H^3/(4 ell^3); both ->1 as H/ell->inf at fixed ell. Meaning consistent. | ANSWER.md §6; results.json limits | FEM at H/ell=0.02,0.05 and 8,16 |
| 5 | Requested numerical values | right | Four ell/H × {gamma_A, gamma_B, kI_A, kI_B}. | submission/results.json rows | 16/16 vs tanh; FEM 16/16 at 1e-6 |
| 6 | Verification and reproducibility | right | Solver 1–8 are 1D consistency (shared reduction). Solver 9–10 are a 2D J route; coordinator replay exit 0 with IntegrationWarning. Reviewer FEM is a separate variational check of the stated energy. Command success is not a scientific verdict. | checks.py; factual_replay/*; review_check_independent_output.txt | C1 omit foundation shifts gamma_A by up to 0.17; C2 incremental traction => v=0 |
| 7 | Physical scope and uncertainty | partly | Lists plane-stress / bonding / linear limits and missing Gc. Omits that k is a Winkler stand-in for a wavenumber-dependent half-space, what isolated-crack calibration would establish, and that dense-array limits are not automatically real-substrate physics. | ANSWER.md §8 | PROBLEM.md item 5 vs submitted scope |

## Solver verification vs reviewer vs coordinator

- Solver: sympy of the 1D wake; quadrature/BVP of the same field; later 2D FD + J-arc (Case A ~1e-4; Case B ~1% at dx=0.02, first-order).
- Coordinator: FACTUAL_REPLAY.json, all listed commands exit 0; IntegrationWarning on j2d.py quad. Not a verdict.
- Reviewer: PROBLEM.md energy assembled as full-width P1 FEM with natural total-traction BCs; no closed-form field inside the solver. Controls C1/C2 detect omitted foundation energy and wrong face BC.

## Independent check details

Code: `review_check_independent.py`. Output: `review_check_independent_output.txt`, `review_check_independent_results.json`.
Command: `timeout 60 python3 review_check_independent.py` (PID bounded by timeout).

- P1 energy I = v(L)-v(0)+int (v'^2+v^2)/2; stationarity => v''=v, v'(0)=v'(L)=-1.
- Convergence at L=1: err 1.79e-6 … 1.74e-9 for n_el=40…1280; ratio 4.02 (O(h^2)).
- Production n_el=4000 vs submission: |d gamma_A| <= 3.6e-9; |d gamma_B| <= 2.7e-8.
- Linear prestress term int v' = -1.52, -0.92, -0.64, -0.49 at the four samples (not zero).
- C1 (drop (k/2)u^2): gamma_A = 0.932, 0.497, 0.333, 0.250 vs 0.762, 0.462, 0.322, 0.245.
- C2 (incremental traction): v≡0, gamma=0.

## Sixteen samples (submission = formula at 1e-12)

ell/H | gamma_A | gamma_B | kI_A | kI_B
0.5 | 0.7615941559557649 | 0.5591607318357129 | 0.8324966558009383 | 0.7133276007351032
1.0 | 0.46211715726000974 | 0.16264015856425462 | 0.6484802333969856 | 0.38471098800719444
1.5 | 0.32151273753163434 | 0.06024252971535857 | 0.5409034952316237 | 0.23413821140722907
2.0 | 0.24491866240370913 | 0.027720167547408525 | 0.4720974293378172 | 0.1588249113588349

## Source recognition

Solver: “No external source was consulted.” Closed forms match Xia & Hutchinson, JMPS 48 (2000) §4.2 (Eqs. 28, 30) as frozen in reference/. Used only to check source/model claims; PROBLEM.md is authoritative.

## Errors / uncertainty

- ANSWER.md incorrectly states that the first-order prestress integral vanishes; the quadratic remainder alone would not give tanh. Code still uses the full functional.
- 2D J replay: IntegrationWarning; Case B not path-grid converged at the tabulated dx.
- Model: Winkler k; no real-film validation. Small-H limits are of this reduction only.
- Model-internal consistency is not real-film validation.

Unresolved: none of the sixteen samples; physical k-calibration / 3D channeling / debonding remain outside the reduced model.

Incomplete obligations: none. Review finished.
