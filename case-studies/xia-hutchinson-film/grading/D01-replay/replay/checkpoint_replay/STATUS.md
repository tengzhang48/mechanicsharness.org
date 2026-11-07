# STATUS.md

UTC time of this update: 2026-10-15T21:25:00Z (session start 21:18:42Z)
Phase: checkpoint_report (minute-5-to-6 reserve).

## Task
Solve PROBLEM.md (periodic cracks in film on elastic foundation). Deliver
ANSWER.md, results.json, checks, REPRODUCE.md; schedule per RUN_PROTOCOL.md.

## Authorized inputs read
TASK.md, PROBLEM.md, RUN_PROTOCOL.md only.

## Derived quantities (done)
- Decay length ell = sqrt(h Ebar/k) = sqrt(h E/(k(1-nu^2))), Ebar=E/(1-nu^2).
- Governing reduced PDE; wake: u1=0, u2=B sinh((x2-W/2)/ell),
  B=-sigma0 ell/(Ebar cosh(W/(2ell))); edges traction free
  (u2'=-sigma0/Ebar, u1'=0). Far-ahead u=0.
- Pure Mode I (KII=0).
- GA=(sigma0^2 ell/Ebar) tanh(H/(2ell)).
- GB=(sigma0^2 ell/Ebar)[2 tanh(H/(2ell))-tanh(H/ell)].
- Normalized: gamma_A=tanh(beta), gamma_B=2 tanh(beta)-tanh(2beta),
  beta=H/(2ell); kI=sqrt(gamma(1-nu^2)).
- Limits: GA~sigma0^2 H/(2Ebar), GB~sigma0^2 H^3/(4 ell^2 Ebar) as H->0;
  both ->sigma0^2 ell/Ebar as H->inf. Critical stresses derived.

## Sampling table (nu=0.3) — saved to results.json (fully populated)
| ell/H | gamma_A | gamma_B | kI_A | kI_B |
|------|---------|---------|------|------|
| 0.5 | 0.7615941560 | 0.5591607318 | 0.8324966558 | 0.7133276007 |
| 1.0 | 0.4621171573 | 0.1626401586 | 0.6484802334 | 0.3847109880 |
| 1.5 | 0.3215127375 | 0.0602425297 | 0.5409034952 | 0.2341382114 |
| 2.0 | 0.2449186624 | 0.0277201675 | 0.4720974293 | 0.1588249114 |

## Actual checks executed and saved outputs (all PASS)
Saved in check_output.txt (checks.py, timeout 180):
1. energy quadrature vs closed form: max rel 1.5e-16 (tol 1e-9)
2. wake ODE residual: <7e-7 (tol 1e-4)
3. traction-free edges: 0.0 (tol 1e-10)
4. independent numerical BVP (solve_bvp) vs closed form: <1e-10 (tol 1e-6)
5. integration-by-parts identity: rel 1e-16 (tol 1e-8)
6. small-H limits gamma_A~beta, gamma_B~2beta^3: match (tol 1e-6)
7. large-H limits ->1: match (tol 1e-9)
8. G=K^2/E consistency: exact (tol 1e-12)
All 8 checks share the same 1D reduction: consistency checks, not an
independent solution route.

## Files saved
- STATUS.md, PROVISIONAL_REPORT.md
- ANSWER.md (full derivation, results, limits, checks, limitations)
- results.json (rows + formulas + limits + checks, fully populated)
- compute.py, checks.py, check_output.txt, results_analytic.json
- REPRODUCE.md (single command)

## Incomplete / uncertainty
- Independent 2D verification (J-integral / boundary element) of the energy
  coefficient NOT executed. The energy prefactor rests on the 1D steady-state
  energy-difference argument (sound, but not independently cross-solved in
  2D). Explicitly flagged in ANSWER.md §8 and results.json "unresolved".
- No known numerical discrepancies; all executed checks PASS.
