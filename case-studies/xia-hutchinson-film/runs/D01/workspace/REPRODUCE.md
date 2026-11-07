# REPRODUCE.md

Requirements: Python 3 with `numpy`, `scipy`, `sympy` (all standard).

Single reproduction command (run in this folder):

    timeout 180 python3 compute.py && timeout 180 python3 checks.py

- `compute.py` derives `ell`, the wake fields, `G_A`, `G_B`, the normalized
  `gamma_A/gamma_B/kI_A/kI_B` and the limits symbolically; writes
  `results_analytic.json`; prints the sampling table.
- `checks.py` independently re-checks the wake ODE, boundary conditions, the
  energy integral (direct quadrature), an independent numerical BVP solve,
  the integration-by-parts identity and the limits; prints PASS/FAIL and
  writes to `check_output.txt` when piped through `tee`.

Both scripts are bounded by `timeout 180` (no free-running processes). The
final machine-readable answer is `results.json`; the derivation is `ANSWER.md`.

Independent 2D J-integral verification (checks 9-10 in ANSWER.md):

    timeout 300 python3 check_j2d_cases.py     # Case A + Case B
    timeout 600 python3 check_B_conv.py        # Case B grid convergence
    timeout 500 python3 check_A_H.py           # Case A H-dependence (tanh shape)

`j2d.py` is shared by those. Outputs are saved to
`check_j2d_cases_output.txt`, `check_B_conv_output.txt` and
`check_A_H_output.txt`.
