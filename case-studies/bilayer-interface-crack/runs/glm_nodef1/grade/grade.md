```json
{
  "theory_source": {
    "answer": "own_computation",
    "evidence": [
      {"file": "submission/method.md", "lines": "51-80", "quote": "Phases and eps (numerical, Q8 FEM)... phase A ... 52.27 deg; phase B ... 4.25 deg; phase change per decade = eps ln 10 = -8.5136 deg"},
      {"file": "submission/code/final_run.py", "lines": "58-76", "quote": "y = T[m]*sqrt(2*pi*rs)*exp(-1j*eps*log(rs)); ... outD ... np.conj(Ad/Q)"},
      {"file": "submission/code/final_audit.txt", "lines": "36-41", "quote": "phase A (consensus median) = +52.2718 deg; phase B (consensus median) = +4.2506 deg; combined phase (exact |K|s) = +30.4330 deg"},
      {"file": "submission/code/fem.py", "lines": "148-229", "quote": "Structured Q8 mesh; interface nodes duplicated for x1 < 0 (crack)... End x1 = -L loaded with the Saint-Venant stress distribution"}
    ],
    "pieces": {
      "phase_change": {"used": "Analytic bimaterial index eps = ln[(1-beta)/(1+beta)]/(2 pi) = ln(2/3)/(2 pi) = -0.06453178, so the change per decade is eps*ln(10) degrees. The FEM log-slope of arg[T sqrt(2 pi r)] is used only as a check that the slope is about -0.0645, not as the reported number.", "equals_key": "yes"},
      "energy_release_rate": {"used": "Exact far-field energy of both layers with plane-strain moduli Ebar = E/(1-nu^2): G = sum [N_i^2/(2 Ebar_i t_i) + 6 M_i^2/(Ebar_i t_i^3)], with (N1,M1)=(-P,-M) and (N2,M2)=(P, M+P(h+H)/2). Case A 1.771604938, case B 1.259259259.", "equals_key": "yes"},
      "amplitude_from_G": {"used": "G = |K|^2 (1/Ebar1 + 1/Ebar2) / (2 cosh^2(pi eps)) = |K|^2 * (64/135) / E2, so A = sqrt(G/c). The cosh factor was checked by a numerical J-integral of the analytic crack-tip eigenpair.", "equals_key": "yes"},
      "phase_force": {"used": "Own Q8 FEM of load case A; K from a least-squares fit of T(r) sqrt(2 pi r) r^{-i eps} = K + c sqrt(r) + d r (and an independent crack-opening fit), then the median of the finest-mesh windows. With r measured in units of h, phase at r=h is arg K: 52.2718 deg.", "equals_key": "yes"},
      "phase_moment": {"used": "Same FEM pipeline run independently on load case B (P=0, M=1), not inferred from case A. Consensus median 4.2506 deg.", "equals_key": "yes"},
      "both_loads": {"used": "Superposition of the complex intensities with the exact amplitudes and the two FEM phases: K_tot = |K_A| exp(i psi_A) + |K_B| exp(i psi_B) at M = P h, giving 30.433 deg. A direct (P=M=1) FEM run is coded as a check of superposition.", "equals_key": "yes"}
    }
  },
  "claimed_checks": {
    "answer": "some_unsupported",
    "checks": [
      {"claim": "Q8 shape-function derivatives match finite differences to 1e-10, a uniform-strain patch test holds to 1e-10, and K*1 = 0.", "claimed_at": "submission/method.md:84-85", "status": "supported", "evidence": "submission/code/selftest.py:9-45 implements all three and asserts the tolerances; no saved output."},
      {"claim": "Homogeneous equal-arm moment benchmark: FEM K_I = 3.4634 against exact 2*sqrt(3) (-0.02%), phase -0.014 deg against exact 0.", "claimed_at": "submission/method.md:86-89", "status": "unsupported", "evidence": "submission/code/bench_homog.py:42-47 and submission/code/validations.py:48-55 run this comparison, but neither prints the claimed 3.4634 or -0.014 deg, and no output file contains them."},
      {"claim": "Homogeneous equal-arm force benchmark: FEM |K| = 2.6460 against exact sqrt(7) (+0.01%).", "claimed_at": "submission/method.md:90-91", "status": "unsupported", "evidence": "submission/code/bench_homog.py:49-52 expects |K|=1 and phase +-90, which is the wrong exact answer for h=H (the moment arm P(h+H)/2 is omitted). validations.py:40-46 has the same wrong target. No code computes sqrt(7) or prints 2.6460."},
      {"claim": "The G-K factor matches the closed form S/(2 cosh^2(pi eps)) to 5e-7 for two material pairs.", "claimed_at": "submission/method.md:92", "status": "supported", "evidence": "submission/code/check_c.py:16-63 numerically integrates the J-integral of the eigenpair and prints the ratio to S/2*sech^2 for (E,nu)=(9,1/3) and for E=(3,5) at two Poisson ratios. No output of that run is saved."},
      {"claim": "Traction-fit and crack-opening-fit K agree to 0.1% in magnitude and 0.05 deg in phase.", "claimed_at": "submission/method.md:93", "status": "contradicted", "evidence": "submission/code/final_audit.txt:9-20,22-34: on the finest mesh the two fits differ by about 0.04 deg (case A, window 0.001-0.02: 52.2490 vs 52.2850) but by about 0.16% in |K| (1.930090 vs the delta-fit magnitude), and coarser levels differ by up to 0.3 deg. The claimed 0.1% is not met."},
      {"claim": "FEM |K| agrees with the exact-relation |K| to 0.2%.", "claimed_at": "submission/method.md:94", "status": "supported", "evidence": "submission/code/final_audit.txt:4-5,17,31: exact 1.93312678 and 1.62980060; finest T-fit magnitudes 1.930090 and 1.627269, i.e. -0.16% and -0.15%, matching method.md:72-73."},
      {"claim": "Superposition of cases A and B was checked against a direct (P=M=1) FEM run, agreeing to 6 digits.", "claimed_at": "submission/method.md:78-79,95", "status": "unsupported", "evidence": "submission/code/validations.py:57-67 codes the direct-versus-sum comparison and would print both complex K values, but it does not check or print a 6-digit residual, and no output of that run is in the folder."},
      {"claim": "Domain length L = 15, 30, 45 changes |K| by less than 0.02% and the phase by less than 0.001 deg.", "claimed_at": "submission/method.md:59-60,95", "status": "unsupported", "evidence": "submission/code/validations.py:78-81 runs case B at L=15/30/45 and prints K, but nothing in the folder records the three values or the claimed differences."},
      {"claim": "The linear solver residual is below 1e-4.", "claimed_at": "submission/method.md:95", "status": "unsupported", "evidence": "submission/code/stab1.py:28-35,53 computes and prints ||Ku-f||/||f|| per mesh, but no printed residual is saved."},
      {"claim": "The displacement-jump amplitude matches the eigenpair constant Q to 0.1%.", "claimed_at": "submission/method.md:96", "status": "contradicted", "evidence": "The opening fit is exactly K_D = conj(A_delta/Q) in submission/code/final_run.py:69-75, so agreement with Q is the same comparison as the T-fit versus delta-fit. final_audit.txt:17 and :31 show |K_T| and |K_D| differing by about 0.16% on the finest mesh, beyond the claimed 0.1%."},
      {"claim": "eps measured from the FEM log-slope of arg[T sqrt(2 pi r)] is about -0.0645 and from arg[delta/sqrt(r)] about +0.0645, confirming the analytic value.", "claimed_at": "submission/method.md:74-76", "status": "unsupported", "evidence": "submission/code/fits.py:25-35,51-74 and submission/code/run_bimaterial.py:54-65 implement the scans and would print them, but bimaterial_raw.json stores eps_T/eps_d of -0.15 and -0.05 (a failed coarse scan), and no file records the claimed +/-0.0645 slopes."}
    ]
  },
  "error_location": {
    "answer": "none",
    "first_wrong_step": "",
    "where": "",
    "explains_groups": []
  },
  "notes": "All eight score.md checks pass. The reported G and amplitudes are the closed forms (both layers, plane strain, bimaterial cosh factor), and the three phases are the solver's own Q8 FEM results: consensus medians 52.2718 and 4.2506 deg, within 0.03 deg of the reference, combined by superposition. The phase change is the analytic eps*ln(10), with the FEM slope used only as a confirmation. Several checks named in method.md have code but no saved output, so they are unsupported; two (T-fit versus opening-fit to 0.1%, and the jump amplitude versus Q to 0.1%) are contradicted by final_audit.txt, where the finest-mesh magnitudes differ by about 0.16%. The homogeneous force benchmark in method.md quotes |K|=sqrt(7) and 2.6460, but bench_homog.py still tests the incorrect target |K|=1. The abandoned central-crack files and the truncated stab2_out.json are not used for the reported numbers. Large binary files under tmp/ were left out of the copy and are irrelevant."
}
```
