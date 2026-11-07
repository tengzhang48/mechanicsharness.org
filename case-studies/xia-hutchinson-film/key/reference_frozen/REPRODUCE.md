# REPRODUCE.md — how to reproduce the reference numbers

## Single deterministic command

```
python3 reproduce.py
```

Run from this folder (the script writes `reference.json` and `reproduce_output.txt` here and prints
the full report to stdout). Runtime ≈ 1.5 s. Requirements: python3 with numpy, scipy (solve_bvp,
solve_banded); sympy (optional — symbolic pre-check; skipped gracefully if absent). No network, no
randomness, no arguments. Exit code 0 iff all production checks pass AND every deliberately-wrong
control is confirmed to deviate by a clear margin.

## What it does (verification routes; shared assumptions)

1. **Symbolic** (SymPy): verifies the exact strip field `v = sinh(a-eta)/cosh(a)` satisfies `v''=v`,
   `v'(0) = -1` (total traction-free crack face: `Delta_sigma22 = -sigma0`), `v(a) = 0` (midline
   antisymmetry); that `u = (sigma0*ell/Ebar)v`, `x2 = ell*eta` solves the dimensional equilibrium
   `h*Ebar*u_{,x2x2} = k u` with `k = Ebar*h/ell^2`; that the problem's stated incremental functional
   reduces exactly to `ell*(h sigma0^2/Ebar)*[v' + (v'^2+v^2)/2] deta`; and that the period integral
   `I = Int_0^{2a}[v' + (v'^2+v^2)/2] deta = -tanh(a)`.
2. **Collocation strip BVP** (`scipy.integrate.solve_bvp`): solves `v''=v, v'(0)=-1, v(a)=0` for
   `a = S/(2 ell)`; energy from Simpson quadrature of the stated functional using only the numerical
   `v`. Then `gamma_A(S) = -I(S)` and `gamma_B = 2*gamma_A_num(H/ell) - gamma_A_num(2*H/ell)`.
3. **Independent finite-difference BVP** (2nd-order, ghost node for `v'(0)=-1`, tridiagonal solve):
   same energy quadrature; used at N = 80…1280 for a convergence study (clean O(N^-2) decay) and at
   N = 40000 in the production table.

Shared assumptions between routes: plane-stress film on Winkler foundation; `u1 = 0` far from tips
(mode-I symmetry); **total** traction `sigma.n = 0` on faces (not incremental); ahead/behind states =
bounded equilibria of the model; period/tip/thickness normalization per REFERENCE.md Secs. 2–3;
fixed h, E, nu, k, sigma0. The numerical routes never use the closed-form strip solution.

## Actual outputs (2026-10-15T12:23Z run, exit 0; also saved in reproduce_output.txt)

```
[sym] symbolic checks: ODE=True face=True midline=True equilibrium(k=Ebar*h/ell^2)=True functional_reduction=True period_integral=-tanh(a):True

  S/ell=  0.5000  colloc_max_err=9.714e-16  gamma_num_err=0.000e+00  fd_gamma_err(N=80..1280)=7.61e-07 1.93e-07 4.85e-08 1.22e-08 3.05e-09
  S/ell=  0.6667  colloc_max_err=3.331e-16  gamma_num_err=0.000e+00  fd_gamma_err(N=80..1280)=1.75e-06 4.43e-07 1.11e-07 2.80e-08 7.00e-09
  S/ell=  1.0000  colloc_max_err=4.607e-15  gamma_num_err=5.551e-17  fd_gamma_err(N=80..1280)=5.42e-06 1.37e-06 3.46e-07 8.67e-08 2.17e-08
  S/ell=  2.0000  colloc_max_err=7.994e-15  gamma_num_err=2.220e-16  fd_gamma_err(N=80..1280)=2.98e-05 7.58e-06 1.91e-06 4.79e-07 1.20e-07
  S/ell=  3.0000  colloc_max_err=8.549e-15  gamma_num_err=5.551e-16  fd_gamma_err(N=80..1280)=6.63e-05 1.69e-05 4.27e-06 1.07e-06 2.69e-07
  S/ell=  4.0000  colloc_max_err=6.217e-15  gamma_num_err=2.220e-16  fd_gamma_err(N=80..1280)=1.10e-04 2.81e-05 7.11e-06 1.79e-06 4.48e-07
[BVP] checks passed: collocation max error 8.549e-15 < 1e-08; energy residual 5.551e-16; FD(1280) 4.484e-07

    ell/H     H/ell              gamma_A              gamma_B                 kI_A                 kI_B
   0.5000  2.000000   0.7615941559557649   0.5591607318357128   0.8324966558009383   0.7133276007351032
   1.0000  1.000000   0.4621171572600097   0.1626401585642546   0.6484802333969856   0.3847109880071944
   1.5000  0.666667   0.3215127375316343   0.0602425297153585   0.5409034952316238   0.2341382114072289
   2.0000  0.500000   0.2449186624037091   0.0277201675474085   0.4720974293378172   0.1588249113588349
[MAIN] max |numerical - analytic| over all samples and routes: 2.144e-10 < 5e-09  -> PASS
```

## Controls (deliberately wrong variants — actual failure outputs preserved, NOT production answers)

```
  C1 (omit foundation energy): wrong gamma_A at ell/H=[0.5, 1.0, 1.5, 2.0] -> ['0.932404', '0.496564', '0.332831', '0.249876']
     correct gamma_A = ['0.761594', '0.462117', '0.321513', '0.244919'] ; max deviation 1.708e-01  -> CONTROL FAILS AS EXPECTED
  C2 (incremental instead of total traction on faces): collapses to v = 0,
     wrong gamma_A = -0.000e+00 (correct 0.462117); max|v| = 0.000e+00  -> CONTROL FAILS AS EXPECTED
  C3 (initial/final spacings swapped: 2*G'(2H) - G'(H)): wrong gamma_B at ell/H=[0.5, 1.0, 1.5, 2.0] -> ['1.166461', '1.061071', '0.844053', '0.679316']
     correct gamma_B = ['0.559161', '0.162640', '0.060243', '0.027720']  -> CONTROL FAILS AS EXPECTED
```

What each control can detect (error classes):

- **C1** — dropping the foundation energy `(k/2)u^2` from the stated functional (e.g., treating the
  film as free-standing, or forgetting the substrate's stored energy): biases gamma_A by up to +22%
  here and would reverse the sign of the wake energy density at smaller H/ell.
- **C2** — imposing the incremental traction `Delta_sigma.n = 0` instead of the total
  `sigma.n = 0` on crack faces: the strip problem collapses to `v ≡ 0` and zero energy release;
  catches any confusion between prestressed reference and total stress.
- **C3** — confusing the initial (ahead, spacing 2H) and final (behind, spacing H) strip states in
  the Case-B bookkeeping (computing `2*G'(2H) - G'(H)`): errors up to a factor ~24 in gamma_B; also
  catches misidentification of which tips move and of the period (2H per advancing new tip).

## Recommended comparison tolerances (from observed numerical accuracy)

- Normalized `gamma_A`, `gamma_B`, `kI_A`, `kI_B`: **1e-6 absolute** (worst observed cross-route
  mismatch 2.2e-10; sqrt roughly halves relative error for kI).
- Dimensional G: **1e-6 relative** on the model formula; in practice physical-input uncertainty
  (substrate stiffness wavenumber dependence, k calibration, ell/h determination) dominates.
- No experimental grading or physical validation of a real film is performed; this is reference
  preparation only (see REFERENCE.md Secs. 7 and 9).
