# PROVISIONAL REPORT

- UTC time: 2026-10-15 14:55:27 UTC (session started 14:37:25 UTC; wall-clock gap due
  to a session interruption between the numeric evaluation and this report).
- Session phase: checkpoint_report (reserved reporting window, before minute 15).

## What has been derived (analytic, complete derivations in progress of write-up)

Model: plane-stress film (h, E, nu=0.3), prestress sigma0 I, Winkler foundation
f = -k u. Equilibrium h sigma_{ab,b} - k u_a = 0 reduces to
h Ebar/2 (1-nu) nabla^2 u + h Ebar/2 (1+nu) grad div u - k u = 0, Ebar = E/(1-nu^2).

1. Decay length: fundamental solutions exp(lambda x1 + i q x2) give two branches,
   lambda^2 = q^2 + 2k/(h Ebar (1-nu)) = q^2 + k/(h mu) (divergence-free) and
   lambda^2 = q^2 + k/(h Ebar) (curl-free). The slowest decay (q=0, curl-free,
   compatible with the symmetry conditions ahead of the tips) sets
   ell = sqrt(h Ebar/k). This is the unique length in the q=0 longitudinal equation
   and is used for all normalizations.
2. Far-ahead state Case A: u = 0 (prestressed). Far-behind: 1D relaxed wake on each
   strip of width H: u1 = 0, u2(x2) = -(sigma0 ell/Ebar) sinh((x2-H/2)/ell)/cosh(H/(2ell)),
   sigma22 = sigma0[1 - cosh((x2-H/2)/ell)/cosh(H/(2ell))], sigma11 = sigma0 - nu(sigma0 - sigma22),
   sigma12 = 0. Edge conditions: both crack faces traction-free (sigma22 = sigma12 = 0);
   ahead of the tip on the crack line: symmetry (u2 = 0, sigma12 = 0).
3. Energy release via steady-state translation / Eshelby J-integral argument:
   G per advancing tip = (W_ahead - W_behind)/h with W = strip integral of the
   DeltaPi density over one period. Strip energy: W(S) = -h sigma0^2 ell
   tanh(S/(2ell))/Ebar (linear prestress term contributes -2x the quadratic term).
   - Case A (period H, one tip per period):
     G_A = (sigma0^2 ell/Ebar) tanh(H/(2ell)).
   - Case B (period 2H, one advancing tip per period; ahead = relaxed width-2H state,
     behind = two relaxed width-H strips):
     G_B = (sigma0^2 ell/Ebar)[2 tanh(H/(2ell)) - tanh(H/ell)]
         = (sigma0^2 ell/Ebar) 2 t^3/(1+t^2), t = tanh(H/(2ell)).
4. Symmetry => pure mode I, KII = 0; KI = sqrt(E G), so
   KI/(sigma0 sqrt(ell)) = sqrt((1-nu^2) gamma), gamma = G Ebar/(sigma0^2 ell).
5. Critical stress sigma0c = sqrt(Ebar Gc/(ell gamma)). Limits:
   H/ell -> 0: Case A gamma ~ H/(2ell) (full relaxation, sigma0c ~ H^{-1/2} diverges);
   Case B gamma ~ (H/ell)^3/4 (already-relaxed film, sigma0c ~ H^{-3/2} diverges faster).
   H/ell -> infinity: both gamma -> 1 (isolated crack, sigma0c -> sqrt(Ebar Gc/ell)).

## Checks actually executed and outcomes

- analytic_checks.py (SymPy, ran at ~14:40 UTC, output in analytic_output.txt):
  ALL 8 CHECKS PASS: wake ODE residual = 0; amplitude C = -ell sigma0/(Ebar
  cosh(S/2ell)); sigma22(S) = 0 exactly; sigma12 = 0; direct symbolic strip integral
  W(S) equals -h sigma0^2 ell tanh(S/(2ell))/Ebar (exact simplify + 4-point numeric
  sampling); linear term = -2 h sigma0^2 ell tanh/Ebar; quadratic term = +1 h sigma0^2
  ell tanh/Ebar; identity 2tanh(x)-tanh(2x) = 2t^3/(1+t^2).
- numeric_eval.py (NumPy, output in numeric_output.txt, rows in numeric_rows.json):
  ell/H = 0.5: gamma_A = 0.761594, gamma_B = 0.559161, kI_A = 0.832497, kI_B = 0.713328
  ell/H = 1.0: gamma_A = 0.462117, gamma_B = 0.162640, kI_A = 0.648480, kI_B = 0.384711
  ell/H = 1.5: gamma_A = 0.321513, gamma_B = 0.060243, kI_A = 0.540903, kI_B = 0.234138
  ell/H = 2.0: gamma_A = 0.244919, gamma_B = 0.027720, kI_A = 0.472097, kI_B = 0.158825
  (internal assert: two algebraic forms of gamma_B agree to < 1e-15)
- NOT yet executed: independent 2D finite-difference solution of the crack problem
  with J-integral evaluation (the main independent verification route). The symbolic
  checks verify internal algebra only; they share the same analytic assumptions, so G
  itself is not yet independently confirmed.

## Incomplete work

- Finite-difference 2D solve (Cases A and B) + path-independence of J (script not yet
  written/run).
- ANSWER.md (full derivation write-up), results.json (final assembly with checks list),
  REPRODUCE.md, QUESTIONS.md (interpretation notes).
- Limit-behavior numeric cross-checks against the tanh formulas.
- Final reporting pass.

## Remaining uncertainty / open items

- "Decay length ell" interpreted as sqrt(h Ebar/k); the wake-approach behind the tip
  involves strip eigenmodes with faster-or-equal decay — will be stated explicitly.
- G formulas rest on the steady-state/J argument; pending independent FD confirmation.
- KI normalizations use G = KI^2/E exactly as stated in PROBLEM.md (mode I only).
- All numbers use nu = 0.3 per instructions; gamma values are nu-independent.
