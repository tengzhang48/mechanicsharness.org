```json
{
  "theory_source": {
    "answer": "own_computation",
    "evidence": [
      {"file": "submission/method.md", "lines": "67-96", "quote": "The phases were obtained from a finite-element solution... The ratio phase was extrapolated linearly in r using the smallest resolved points"},
      {"file": "submission/work/fem.py", "lines": "60-221", "quote": "def build_mesh... def assemble... def apply_loads... def solve"},
      {"file": "submission/work/run_better.py", "lines": "3-10", "quote": "for name,P,M in [('A',1,0),('B',0,1)]: m=fem.build_mesh(x0=1e-6,q=1.07,Lx=3.0) ... np.savez('better%s.npz'%name"},
      {"file": "submission/work/analyze_better.py", "lines": "21-27", "quote": "ratio=(DY+1j*DX)/(DYs+1j*DXs) ... ph_h"},
      {"file": "submission/results.json", "lines": "5-12", "quote": "\"phase_at_h_deg\": 51.6 ... \"phase_at_h_deg\": 5.4 ... \"phase_at_h_both_deg\": 30.6"}
    ],
    "pieces": {
      "phase_change": {"used": "Exact bimaterial formula eps*ln(10), with eps=(1/(2*pi))*ln((1-beta)/(1+beta)), beta=0.2, eps=-0.0645318; not fitted from the FE field. Verified by an 8x8 Airy eigenfunction determinant vanishing at s=-1/2 +/- i*eps.", "equals_key": "yes"},
      "energy_release_rate": {"used": "Path-independent J of the far-behind beam stresses only (ahead stress-free), plane-strain moduli E*=E/(1-nu^2), both layers: case A G E2 h/P^2=287/162, case B G E2 h^3/M^2=34/27.", "equals_key": "yes"},
      "amplitude_from_G": {"used": "Two-material relation G=[(kappa+1)/mu1+(kappa+1)/mu2]|K|^2/(16 cosh^2(pi eps)) with C_G=64/135, A=|K|.", "equals_key": "yes"},
      "phase_force": {"used": "Own Q2 FEM of the bilayer; crack-opening jump delta=(uy+-uy-)+i(ux+-ux-) ratioed to a manufactured analytic eigenfunction, phase extrapolated in r and shifted by eps*ln(h) with h=0.5 in their units (H=1). Reported 51.6 deg. The npz outputs were left out of the copy, so the rounding to 51.6 is not shown as a number in a saved log.", "equals_key": "no"},
      "phase_moment": {"used": "Same FEM and jump-ratio extrapolation on its own case-B solve (P=0, M=1), not inferred from case A. Reported 5.4 deg.", "equals_key": "no"},
      "both_loads": {"used": "Superposition arg(K_A + h K_B) with h=0.5, claimed to match a direct combined-load jump phase at 30.6 deg. No driver applies both loads together.", "equals_key": "yes"}
    }
  },
  "claimed_checks": {
    "answer": "some_unsupported",
    "checks": [
      {"claim": "Dundurs exponent independently verified: the 8x8 Airy eigenfunction determinant vanishes at s=-1/2 +/- i*eps", "claimed_at": "submission/method.md:103-105", "status": "supported", "evidence": "submission/work/eigen.py:88-90 builds the matrix and prints det and min singular value at those s; analytic.py:28-50 and eigen2.py:48-53 do the same SVD nullspace."},
      {"claim": "Analytic eigenfunction interface traction is exactly a single power r^s at theta=0", "claimed_at": "submission/method.md:106-107", "status": "supported", "evidence": "submission/work/eigen2.py:71-77 evaluates T/r^s at r=0.1,1,3 and prints it; analytic.py:69-73 prints KA, KB from theta=0 stress."},
      {"claim": "Far-end beam tractions are self-equilibrated (net force and moment vanish to machine precision)", "claimed_at": "submission/method.md:108-109", "status": "unsupported", "evidence": "fem.py:160-204 applies the beam tractions but never sums force or moment. subproblem.py:72 prints net force only for the manufactured eigenfunction box, not the beam loads, and does not check moment."},
      {"claim": "FE far-field stresses reproduce the applied beam stresses (layer 1 to 1e-6; layer 2 to a few percent)", "claimed_at": "submission/method.md:110-111", "status": "unsupported", "evidence": "No code samples stresses at the loaded end or compares them with the beam formula. post.py extracts only the interface ahead of the tip and the crack-face jump."},
      {"claim": "Domain J of the FE solution is 3.6-3.75 for case A vs exact 3.543, and 10.078 for case B on r=0.05-0.2 vs exact 10.074", "claimed_at": "submission/method.md:112-115", "status": "cannot_tell", "evidence": "The check is coded (jdomain.py:3-48; analyze_better.py:13; checkB.py:5-6; analyzeA.py:15) and prints J, but the printed numbers are not saved as text. The npz inputs (betterA.npz, betterB.npz, fine_data.npz) were left out as too large, so the quoted 3.6-3.75 and 10.078 cannot be confirmed or contradicted."},
      {"claim": "Crack-opening amplitude ratios match energy-derived |K| ratios: case A expects 1.946, measured 1.94-1.99; case B expects 3.282, measured 3.35-3.6", "claimed_at": "submission/method.md:116-118", "status": "cannot_tell", "evidence": "analyze_better.py:22-27 and analyzeA.py:24-29 compute the jump-ratio magnitude and print it, but those prints were not copied and the npz files were left out, so the stated measured ranges cannot be verified."},
      {"claim": "Several meshes compared; extrapolated phases agree within about one degree", "claimed_at": "submission/method.md:119-120", "status": "cannot_tell", "evidence": "run_better.py (x0=1e-6, q=1.07) and run_fine.py (x0=1e-6, q=1.15) build two meshes and analyze_better.py/analyzeA.py print phases, but no saved output of the comparison remains and the npz results were left out."},
      {"claim": "Combined-load phase from K_A + h K_B agrees with a directly measured combined-load jump phase (30.6 deg), to better than 0.1 deg", "claimed_at": "submission/method.md:120-121", "status": "unsupported", "evidence": "No script superposes the two complex K's or solves M=P*h. run_better.py and run_fine.py run only (P,M)=(1,0) and (0,1). The 30.6 appears only as a written result in method.md:134 and results.json:12."}
    ]
  },
  "error_location": {
    "answer": "numerics",
    "first_wrong_step": "The closed forms (eps, G, |K| from G, and the (r/h)^{i eps} normalisation) are right and the phases come from the solver's own FEM, but the crack-opening phase is only averaged over r in [3e-5, 1e-3] and then shifted by eps*ln(h), not extrapolated as claimed, on a model only Lx=3H long whose case-A domain J is itself several percent off. That leaves the reported phases about 0.65 deg and 1.1 deg from the reference.",
    "where": "submission/work/analyze_better.py:24-27 (mean phase over 3e-5 to 1e-3, plus eps*ln(0.5)); submission/method.md:94-96 and 138-140 (the +/-1 deg extrapolation that was written into results.json:5 and :10)",
    "explains_groups": ["case A: phase at r = h (deg)", "case B: phase at r = h (deg)"]
  },
  "notes": "Only the two single-load phases fail (51.6 vs 52.249 and 5.4 vs 4.279); their difference is 46.2 deg against the exact 47.97 deg, so the error is not a common phase offset. The both-loads phase 30.6 passes the 0.5 deg tolerance even though the inputs are those two wrong phases (reference 30.433). G and the amplitudes are the exact fractions and do not depend on the FEM. notes.md:41-62 is an abandoned wrong intermediate (G_A E2 h/P^2 = 4.140625 and amplitudes 2.955 and 1.834) that was not what was reported; method.md:43-61 and results.json have the correct closed forms. The phase numbers 51.6, 5.4 and 30.6 are rounded claims: no copied text file prints them, and every npz the analysis scripts read was left out of the copy. Quarter-point elements (quarter.py) are written but never called by the run scripts."
}
```
