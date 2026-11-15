```json
{
  "theory_source": {
    "answer": "own_computation",
    "evidence": [
      {"file": "submission/method.md", "lines": "78-110", "quote": "`solve.py` implements displacement finite elements... I fit ... sqrt(2 pi r) T(r) exp[-i epsilon log(r/h)] = K_h + c1 r + c2 r^2 over 1e-4 <= r/h <= 1e-2. The constant term gives the phase at r = h"},
      {"file": "submission/solve.py", "lines": "195-201", "quote": "kh = t[select] * (np.sqrt(2 * np.pi * rr) * np.exp(-1j * EPSILON * np.log(rr)))[:, None]"},
      {"file": "submission/fem_fine.json", "lines": "214-223", "quote": "\"phase\": [52.24897442996617, 4.287414384679429]"},
      {"file": "submission/results.json", "lines": "5-13", "quote": "\"phase_at_h_deg\": 52.25 ... \"phase_at_h_deg\": 4.29 ... \"phase_at_h_both_deg\": 30.44 ... \"phase_change_per_decade_deg\": -8.513574487346496"}
    ],
    "pieces": {
      "phase_change": {"used": "Analytic bimaterial formula: epsilon = (1/(2 pi)) log((1-beta)/(1+beta)) with beta = 0.2, then Delta psi = epsilon * ln(10) in degrees. Not fitted from the FEM; the FEM slope of arg T vs ln r is used only as a check.", "equals_key": "yes"},
      "energy_release_rate": {"used": "Plane-strain beam energy of both layers far behind the tip, G = P^2/(2 Ebar1 h) + 6 M^2/(Ebar1 h^3) + P^2/(2 Ebar2 H) + 6 [M + P(h+H)/2]^2/(Ebar2 H^3), giving 287/162 and 34/27.", "equals_key": "yes"},
      "amplitude_from_G": {"used": "G = ((1-beta^2)/2) (1/Ebar1 + 1/Ebar2) |K_h|^2 = (64/135) |K_h|^2, so A = |K_h|. Includes the bimaterial cosh factor (1-beta^2).", "equals_key": "yes"},
      "phase_force": {"used": "Own Q2 FEM of the plane-strain problem (case P=1, M=0). Traction T = sigma22 + i sigma12 sampled on both sides of the bonded interface, averaged, then least-squares fit of sqrt(2 pi r) T exp(-i epsilon ln r) = K_h + c1 r + c2 r^2 on 1e-4 <= r <= 1e-2. Phase of the constant term is the phase at r = h. Reported 52.25 deg is that fit (52.24897) rounded to 0.01 deg.", "equals_key": "yes"},
      "phase_moment": {"used": "Same FEM and fit, separate load P=0, M=1, not inferred from case A. Reported 4.29 deg is the fit (4.28741) rounded to 0.01 deg.", "equals_key": "yes"},
      "both_loads": {"used": "Complex superposition K = |K_A| exp(i psi_A) + |K_B| exp(i psi_B) with exact energy amplitudes and the two fitted phases, for M = P h (h=1). Reported 30.44 deg is 30.43694 rounded to 0.01 deg.", "equals_key": "yes"}
    }
  },
  "claimed_checks": {
    "answer": "some_unsupported",
    "checks": [
      {"claim": "Direct Gauss integration of the remote stress energy densities agrees with the analytic energies to floating-point precision for cases A, B, and their sum.", "claimed_at": "submission/method.md:132-134", "status": "supported", "evidence": "submission/postprocess.py:16-23 and 63-66 (check_beam_energy, rtol 1e-13); submission/checks.json:5-9 matches the analytic G values."},
      {"claim": "Final FEM amplitudes differ from the exact amplitudes by -0.00251% and -0.02544%.", "claimed_at": "submission/method.md:135-136", "status": "supported", "evidence": "submission/postprocess.py:74-75; submission/checks.json:17-20; amplitudes in submission/fem_fine.json:215-218 versus exact amplitudes in solve.py via G_COEFFICIENT."},
      {"claim": "Final upper-side and lower-side phase estimates differ by only 0.00112 and 0.00151 degrees.", "claimed_at": "submission/method.md:137-138", "status": "supported", "evidence": "submission/postprocess.py:77-79; submission/checks.json:21-24; upper and lower phases in submission/fem_fine.json:35-38 and 127-130."},
      {"claim": "Changing L from 8 to 12 on the refined meshes changes the two phases by 0.00388 and 0.01092 degrees; refining to the final mesh changes them by -0.00406 and -0.01896 degrees relative to the L=12 run.", "claimed_at": "submission/method.md:139-141", "status": "supported", "evidence": "submission/postprocess.py:81-88; submission/checks.json:25-32; phases in fem_refined.json:219-222, fem_length8.json:219-222, fem_fine.json:219-222."},
      {"claim": "Moving the final fit window to [1e-5, 1e-3] changes the phases by -0.00242 and -0.00010 degrees.", "claimed_at": "submission/method.md:142-143", "status": "supported", "evidence": "submission/postprocess.py:90-91; submission/checks.json:33-36; fem_fine.json average_0.0001_0.01 vs average_1e-05_0.001 phases at lines 219-222 and 265-268."},
      {"claim": "Fitting the logarithmic phase slope without imposing epsilon gives epsilon = -0.06451097 and -0.06453256.", "claimed_at": "submission/method.md:144-146", "status": "cannot_tell", "evidence": "Code is in submission/postprocess.py:93-105 and the numbers are stored in submission/checks.json:37-40, but the fit reads fem_fine.npz (sampled interface stresses), which MANIFEST.txt left out as not a text file. The JSON outputs do not contain the raw r and T samples, so the slope fit cannot be confirmed from the copied files."},
      {"claim": "The energy cross term requires cos(phase_A - phase_B) = 1/sqrt(G_A G_B), a difference of 47.97051474 degrees; the numerical difference is 47.96156005 degrees, error -0.00895 degrees.", "claimed_at": "submission/method.md:147-151", "status": "supported", "evidence": "submission/postprocess.py:68-72; submission/checks.json:14-16; phases used are the fem_fine average fit."},
      {"claim": "Superposing the fitted phases with the exact amplitudes recovers the combined-load energy within 0.00690%.", "claimed_at": "submission/method.md:152-153", "status": "supported", "evidence": "submission/postprocess.py:107-110; submission/checks.json:41."},
      {"claim": "The final unconstrained linear-system relative residual is 1.96e-11, and remote force and moment resultants cancel to about 1e-15.", "claimed_at": "submission/method.md:154-156", "status": "unsupported", "evidence": "solve.py:213-217 prints the residual and the force/moment sums, but no saved output contains 1.96e-11 or the resultant values. fem_*.json store work_matrix and elapsed time only. The .npz files that might have held a log were left out, and they are stress samples, not solver logs."}
    ]
  },
  "error_location": {
    "answer": "none",
    "first_wrong_step": "",
    "where": "",
    "explains_groups": []
  },
  "notes": "All 8 score.md checks pass. Phases in results.json are the solver's own Q2 finite-element fits (fem_fine.json, average of upper and lower tractions, window 1e-4 to 1e-2), rounded to 0.01 degree by postprocess.py:42-49; G and the amplitudes are the closed forms, not the FEM. The oscillation epsilon is derived analytically (beta=0.2, material 1 above) and then imposed in the K extraction; a free fit of the slope is only a check. Two claimed checks are not fully backed by the copied files: the free epsilon fit depends on fem_fine.npz, which was left out, and the stated residual 1.96e-11 and resultant cancellation of 1e-15 are printed by solve.py but never stored. The reported combined phase 30.44 comes from exact amplitudes plus fitted phases (30.43694), not from the slightly different both_phase stored inside fem_fine.json (30.43984), which uses the FEM amplitudes."
}
```
