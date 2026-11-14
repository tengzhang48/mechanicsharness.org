```json
{
  "theory_source": {
    "answer": "stated_without_source",
    "evidence": [
      {
        "file": "submission/method.md",
        "lines": "5-6",
        "quote": "No web sources or external data were used."
      },
      {
        "file": "submission/method.md",
        "lines": "21-39",
        "quote": "The first-order in-plane crack-front weight-function relations in an infinite isotropic solid are δK1 = K1^0 J[h], δK2 = p2 K2^0 J[h] - b2 K3^0 h', δK3 = p3 K3^0 J[h] + b3 K2^0 h', where p2=(2-3ν)/(2-ν), p3=(2+ν)/(2-ν), b2=2/(2-ν), b3=2(1-ν)/(2-ν)."
      },
      {
        "file": "submission/method.md",
        "lines": "41-46",
        "quote": "These are the planar crack's linear elastic response, including the change to the local front frame. ... As a sign check, when nu = 0 their coefficients reduce to the local rotation of the two shear components."
      },
      {
        "file": "submission/solve.py",
        "lines": "23-27",
        "quote": "P2 = (2.0 - 3.0 * NU) / (2.0 - NU); P3 = (2.0 + NU) / (2.0 - NU); B2 = 2.0 / (2.0 - NU); B3 = 2.0 * (1.0 - NU) / (2.0 - NU)"
      }
    ],
    "coefficients": {
      "p1": {"used": "1 (δK1 = K1^0 J[h], with J the 1/(2π) PV operator)", "equals_key": "yes"},
      "p2": {"used": "(2-3ν)/(2-ν) = 11/17 at ν=0.3", "equals_key": "yes"},
      "p3": {"used": "(2+ν)/(2-ν) = 23/17 at ν=0.3", "equals_key": "yes"},
      "b2": {"used": "2/(2-ν), sign minus K3^0 h'", "equals_key": "yes"},
      "b3": {"used": "2(1-ν)/(2-ν), sign plus K2^0 h'", "equals_key": "yes"},
      "m3": {"used": "1/(1-ν), via r3^2/(1-ν) in D and in the linearised Irwin relation", "equals_key": "yes"}
    }
  },
  "claimed_checks": {
    "answer": "all_supported",
    "checks": [
      {
        "claim": "Reconstructs CG and SG directly from the SIF coefficients and Irwin's relation, including cancellation of the energy sine coefficient.",
        "claimed_at": "submission/method.md:169-170",
        "status": "supported",
        "evidence": "submission/solve.py:140-151 computes energy_cos and energy_sin from C1,C2,C3 and S1,S2,S3 and asserts equality with CG and zero. No saved run log; the check is in the code, which prints only on success."
      },
      {
        "claim": "Checks the straight-front loading normalization, zero means, periodicity, and reflection symmetry of the front.",
        "claimed_at": "submission/method.md:171-172",
        "status": "supported",
        "evidence": "submission/solve.py:152 and 154-159: D*K10_scaled**2 == 1, mean of f and of the three kappa profiles, profiles(s+1) vs profiles(s), and f(-s) vs f(s), on 4096 points."
      },
      {
        "claim": "Independently sums and solves 96 Fourier harmonics and compares every profile against the closed forms at 4096 positions.",
        "claimed_at": "submission/method.md:173-175",
        "status": "supported",
        "evidence": "submission/solve.py:161-177 builds f, J[f], f' and the three kappa profiles from 96 cosine harmonics and asserts agreement with profiles() within 2e-14. The harmonic coefficients are taken from the closed-form CG, so this re-checks the Fourier sum against the log closed form rather than a fully independent solve."
      },
      {
        "claim": "Evaluates the periodic real-space line integral by adaptive quadrature at all sixteen requested positions and checks the normalization and sign of the line operator against its Fourier implementation.",
        "claimed_at": "submission/method.md:176-180",
        "status": "supported",
        "evidence": "submission/solve.py:102-126 and 179-185: scipy.integrate.quad of the symmetrized kernel (pi/2)/sin(pi u)^2, compared with kappa1 at s=k/16, tolerance 2e-11."
      },
      {
        "claim": "Reconstructs the pointwise first-order energy release rate from all three SIF profiles at 4096 positions and checks it against g.",
        "claimed_at": "submission/method.md:181-182",
        "status": "supported",
        "evidence": "submission/solve.py:187-190: energy_profile from kappa1, kappa2, kappa3 compared with toughness(s), tolerance 2e-14."
      },
      {
        "claim": "Checks the analytic peak-to-peak displacement against a grid containing both front extrema.",
        "claimed_at": "submission/method.md:183-184",
        "status": "supported",
        "evidence": "submission/solve.py:191-195: np.ptp of f on the 4096-point grid (which includes s=0 and s=1/2) versus the closed-form log((1+q)/(1-q))/(pi L)."
      },
      {
        "claim": "Locates the largest mode III value on a 65536-point grid over a full period, refines it with scalar numerical optimization, and compares it with the analytic stationary-point calculation.",
        "claimed_at": "submission/method.md:185-187",
        "status": "supported",
        "evidence": "submission/solve.py:55-72 and 197-217: mode3_extrema solves the stationary equation analytically; the grid argmax is refined with scipy.optimize.minimize_scalar and compared in value (1e-12) and location (1e-8)."
      },
      {
        "claim": "Checks the requested JSON schema, finite numerical values, and sixteen entries per sampled profile.",
        "claimed_at": "submission/method.md:188-189",
        "status": "supported",
        "evidence": "submission/solve.py:219-234: json round-trip with allow_nan=False and exact key/length checks. results.json matches that schema."
      },
      {
        "claim": "All checks passed; Fourier/profile and equilibrium residuals below 2e-14; real-space integral discrepancies below 2e-11; mode III maximum within 1e-12 and location within 1e-8.",
        "claimed_at": "submission/method.md:191-194",
        "status": "supported",
        "evidence": "Those tolerances are the ones asserted in solve.py:177, 184, 190, 212-213. No printed log was saved (MANIFEST lists only method.md, solve.py, results.json, PUBLIC_TASK.md), so the numerical residuals themselves are not in the folder; the code that enforces them is."
      }
    ]
  },
  "error_location": {
    "answer": "none",
    "first_wrong_step": "",
    "where": "",
    "explains_groups": []
  },
  "notes": "score.md is 76/76. The six coefficients and Irwin weight match KEY_NOTES.md section 1, including signs, and the code applies them as in section 2 (C1=-π, CG=-2π L with L=1/ρ, f=log(d)/(2π L), slope terms cancel). They are written down with a local-frame sign argument and a ν=0 rotation check, but there is no derivation from the elasticity equations and no named source (method.md explicitly says none were used), so this is stated_without_source rather than derivation_shown. The internal checks in solve.py are consistency checks of the closed forms, the line operator, equilibrium, and the mode-III maximum; method.md line 197 states that no independent 3D elastic computation was run. No output log was copied, but the checks live in solve.py and the grading rules treat code alone as support."
}
```
