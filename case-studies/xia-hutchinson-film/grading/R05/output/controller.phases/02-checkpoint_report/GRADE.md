# GRADE.md (partial, provisional)

Reviewer independent scientific review. No desired verdict supplied.

Numerical agreement: **16 / 16** (tol 1e-12 abs vs tanh formulas and reference.json rows).

| # | Item | Judgment | Reason | Evidence | Deciding check |
|---|------|----------|--------|----------|----------------|
| 1 | Wake formulation and fields | right | ell=sqrt(h Ebar/k); fields and total-traction edges match PROBLEM.md | submission/ANSWER.md §1–2; compute.py wake BCs | replay: u2'(0)=u2'(H)=-sigma0/Ebar |
| 2 | Energy release | right | GA, GB, one tip/period, /h correct. ANSWER.md wrongly claims linear term integrates to 0; compute.py includes it | ANSWER.md §4; compute.py W_dens | replay GA,GB = tanh forms |
| 3 | Fracture modes and intensities | right | Mode I, G=KI^2/E, kI=sqrt((1-nu^2)gamma) | ANSWER.md §3,5 | kI vs sqrt(0.91*gamma) |
| 4 | Thresholds and limits | right | sigma_c and H->0, inf at fixed h,E,nu,k | ANSWER.md §6; results.json limits | series 2 beta^3 = H^3/(4 ell^3) |
| 5 | Requested numerical values | right | four locations × four quantities | results.json rows | 16/16 vs reference.json |
| 6 | Verification and reproducibility | partly | 1D checks share reduction; 2D J is a separate route (replay OK, IntegrationWarning); reviewer FEM not yet run | checks.py; factual_replay | pending independent FEM |
| 7 | Physical scope | partly | states linear/Winkler limits; omits k-calibration / wavenumber-dependent substrate | ANSWER.md §8 | vs PROBLEM.md item 5 / REFERENCE §7 |

Solver verification vs reviewer vs coordinator:
- Solver: sympy + quadrature/BVP of same 1D wake; later 2D J.
- Coordinator: commands exit 0; not a scientific verdict.
- Reviewer: formula re-derivation and 16-value match; FEM incomplete.

Source recognition: solver claims PROBLEM.md only. Closed forms match Xia & Hutchinson JMPS 2000 §4.2 as in reference/ (used only to check source claims).

Unresolved: independent FEM; foundation vs real substrate.

Incomplete obligations: REPRODUCE.md for reviewer code; final reporting.
