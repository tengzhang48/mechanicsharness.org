The first-order response was meant to come from the solver’s own three-dimensional crack code, not from a written formula. No coefficient of the first-order theory is ever formed, Part B is never attempted, and `results.json` is absent, which is why every scored number is missing.

```json
{
  "theory_source": {
    "answer": "own_computation",
    "evidence": [
      {
        "file": "submission/crack3d.py",
        "lines": "1-10",
        "quote": "Spectral dislocation solver v6 for a half-plane crack x < a(z), y = 0 in an infinite isotropic solid ... Exact quasistatic traction-dislocation operator"
      },
      {
        "file": "submission/crack3d.py",
        "lines": "289-309",
        "quote": "K2 = (self.Ep * s2p / 8.0) * (C[0] - aph * C[2]) * rot ... out['C' + name] = coef[1] * self.lam / (self.A * self.K10)"
      }
    ],
    "coefficients": {
      "p1": {"used": "not formed; K1 taken from the opening density of the 3D solve", "equals_key": "not_stated"},
      "p2": {"used": "not formed; K2 taken from the local-frame shear density", "equals_key": "not_stated"},
      "p3": {"used": "not formed; K3 taken from the local-frame shear density", "equals_key": "not_stated"},
      "b2": {"used": "not formed; slope coupling enters only through the tilted frame n, t", "equals_key": "not_stated"},
      "b3": {"used": "not formed; slope coupling enters only through the tilted frame n, t", "equals_key": "not_stated"},
      "m3": {"used": "not formed; energy release rate is never computed", "equals_key": "not_stated"}
    }
  },
  "claimed_checks": {
    "answer": "some_unsupported",
    "checks": [
      {
        "claim": "Straight crack a=0 recovers K1=1, K2=0.6, K3=0.8, uniform along the front",
        "claimed_at": "submission/validate1.py:7-14",
        "status": "supported",
        "evidence": "validate1.py:7-14 calls solve and prints the three means against the exact values; no output of this run is saved, and sweep3.json line 1 shows the later solver still misses them (K0_1=0.981, K0_2=0.589, K0_3=0.762)"
      },
      {
        "claim": "A constant shift of the front leaves the three factors unchanged (shift invariance)",
        "claimed_at": "submission/validate1.py:18-25",
        "status": "supported",
        "evidence": "validate1.py:18-25, validate2.py:71-72, validate3.py:21-23 and validate4.py:19-22 all run a shifted straight front and compare the extracted K with the exact values; none of those runs is saved"
      },
      {
        "claim": "Straight-crack factors converge as M, Nxf and the collocation margin are varied",
        "claimed_at": "submission/validate2.py:15",
        "status": "supported",
        "evidence": "validate2.py:16-20 and validate3.py:15-20 loop over M, Nxf and the taper and print the K errors; output not saved"
      },
      {
        "claim": "Windowed straight-crack density produces the known 2D traction, ratio ~1 on the crack and ~0 ahead",
        "claimed_at": "submission/dbg1.py:28-32",
        "status": "unsupported",
        "evidence": "dbg1.py imports smooth_chi, which crack3d.py does not define (it defines smootherstep), so the script cannot run as written and no output is saved"
      },
      {
        "claim": "Chebyshev projection of the exact straight-crack density recovers C0",
        "claimed_at": "submission/dbg2.py:25-28",
        "status": "unsupported",
        "evidence": "dbg2.py:3 also imports the missing smooth_chi; no output file"
      },
      {
        "claim": "Direct Chebyshev projection of the windowed density recovers C0",
        "claimed_at": "submission/dbg3.py:24-27",
        "status": "unsupported",
        "evidence": "dbg3.py:3 imports the missing smooth_chi; no output file"
      },
      {
        "claim": "Solved density on a shifted straight front matches the exact sqrt(a-x) density, and the collocation residual is small",
        "claimed_at": "submission/dbg4.py:16-21",
        "status": "unsupported",
        "evidence": "dbg4.py:34-43 uses p.c as a 3D array and calls basis shapes the current solver does not have, so it does not run against crack3d.py; no output saved"
      },
      {
        "claim": "Basis, collocation range and the straight-front K1 are sensible at low resolution",
        "claimed_at": "submission/dbg5.py:8-19",
        "status": "unsupported",
        "evidence": "dbg5.py:10 calls p.basis_grid, which crack3d.py does not define; no output saved"
      },
      {
        "claim": "Pure mode I gives C2=C3=S2=S3=0 and S1=0; pure mode II gives no K1/K3 coupling",
        "claimed_at": "submission/puremodes.py:19-21",
        "status": "supported",
        "evidence": "puremodes.py:18-26 runs the three pure modes and the mixed loading through fit_partA; it writes /work/puremodes.json, which is not in the folder, so the numbers are not shown"
      },
      {
        "claim": "Wavy-front coefficients should approach C1~-pi, S1~0, C2~-0.6 pi, S2~+1.6 pi, C3~-0.8 pi, S3~-1.2 pi",
        "claimed_at": "submission/validate4.py:23",
        "status": "contradicted",
        "evidence": "sweep1.log and sweep3.json are the outputs of that computation and do not match the prediction: sweep3.json at A=0.02 gives C1=-0.583, C2=-0.255, S2=6.759, C3=-0.424, S3=-2.587, and sweep1.log scatters from C1=-14 to +189"
      }
    ]
  },
  "error_location": {
    "answer": "incomplete",
    "first_wrong_step": "No results were ever written. The 3D solver only fits Part A coefficients inside parameter sweeps, those fits are unstable and were not turned into the eight reported numbers, and Part B (the equilibrium front, the kappas and K10_scaled) is not attempted anywhere.",
    "where": "MANIFEST.txt:2; submission/crack3d.py:298-309",
    "explains_groups": [
      "Part A: cosine-front coefficients",
      "Part B: loading level and front amplitude",
      "Part B: front shape f(s)",
      "Part B: kappa1(s)",
      "Part B: kappa2(s)",
      "Part B: kappa3(s)",
      "Part B: maximum of kappa3"
    ]
  },
  "notes": "There is no method.md and no results.json; score.md is 0/76 because every submitted value is None. The only theory in the folder is a spectral dislocation code (crack3d.py) that solves the half-plane crack and extracts K1, K2, K3 in the tilted front frame, so the six first-order coefficients are never written down. Its own outputs show the solve has not converged: sweep1.log has relative residuals of 0.6-0.8 and C1 swinging from -14 to +189, and the ridge-regularized sweep3.json still gives C1 about -0.6 instead of -pi. The comment in validate4.py expects the rotation-only values (C2=-0.6 pi, S2=+1.6 pi, C3=-0.8 pi, S3=-1.2 pi), which omit the nonlocal weights p2, p3 and the elastic slope couplings b2, b3; that expectation is not coded into the operator. Part B does not appear in any file. dbg1-dbg3 import a function crack3d.py does not define, and dbg4/dbg5 call an older interface, so those checks cannot run against the solver that is in the folder."
}
```
