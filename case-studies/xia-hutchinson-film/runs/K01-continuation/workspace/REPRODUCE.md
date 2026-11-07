# REPRODUCE

Single reproduction command (from this working directory):

    python3 run_all.py

This runs, in order:
1. `python3 analytic_checks.py` — SymPy internal-algebra suite (8 checks); output saved to
   `analytic_output.txt`. Detects algebra errors in the wake solution, strip-edge
   conditions, energy integral and the gamma_B identity. Shared assumption with the
   analytics: none beyond the closed forms themselves.
2. `python3 numeric_eval.py` — normalized table; output to `numeric_output.txt` and
   `numeric_rows.json`; asserts the two algebraic forms of gamma_B agree to < 1e-15.
3. `python3 fd_check.py` — independent 2D finite-difference solve (NumPy/SciPy sparse
   direct solve) of the membrane-foundation BVP on the symmetry half-domain; G extracted
   from Eshelby line integrals; 5 configurations plus a d vs d/2 grid-convergence pair;
   output to `fd_output.txt`. Declared tolerance: |rel err| < 3% vs the closed form.
   This is the only expensive step (~tens of seconds to ~2 min depending on the machine).
4. Asymptotic-limit ratios (computed inline in run_all.py): gamma_A/x -> 1,
   gamma_B/(2x^3) -> 1 as x -> 0; gamma_A, gamma_B -> 1 at x = 20.

`run_all.py` then assembles `results.json` (overwriting it with fresh check outputs and
provenance `generated_by: run_all.py`).

Requirements: python3 with numpy, scipy, sympy. No package installation was performed in
either session; the environment already provided these.

Provenance notes:
- The analytic results (ell, wake states, W(S), G_A, G_B, KI normalizations) are closed
  forms; `run_all.py` recomputes every sampled number from them.
- The FD route shares with the analytics only (a) the model PDE and (b) the far-field
  Dirichlet states (bounded strip equilibria — an input to the model). The comparison
  target G_exact in fd_check.py is the closed-form tanh result, so the FD check validates
  the translating-tip energy bookkeeping and the 2D field consistency, not the far states
  themselves and not the physical model.
- Individual steps can be rerun separately with the three commands above; each writes its
  own output file.
