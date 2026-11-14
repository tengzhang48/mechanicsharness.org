```json
{
  "theory_source": {
    "answer": "own_computation",
    "evidence": [
      {"file": "submission/method.md", "lines": "32-43", "quote": "Wiener–Hopf analysis (exact, first_order.py)… The remaining √X-edge part is solved by weighted least squares… LS residual 1.7e-4"},
      {"file": "submission/first_order.py", "lines": "83-135", "quote": "def solve_shear(ML=6, Ncol=48…): … sols.append((s, res)) … R[i, j] = S2PI*(DVEC[j]/2.0)*A12m[i, j]"},
      {"file": "submission/method.md", "lines": "51-62", "quote": "δK₂^g = |k|δâ [(−1/2+γ)K₂⁰ + iσγ K₃⁰] … K₂^prob = K₂^g − a'K₃⁰ … verified against the exactly solvable inclined-front anti-plane problem"},
      {"file": "submission/solve_AB.py", "lines": "27-32", "quote": "GAM = NU/(2 - NU); D2 = -0.5 + GAM; D3 = -0.5 - GAM; T2p = D2*R2 + 1j*GAM*R3 - 1j*R3"}
    ],
    "coefficients": {
      "p1": {"used": "1", "equals_key": "yes"},
      "p2": {"used": "(2-3ν)/(2-ν) = 1-2γ with γ=ν/(2-ν)", "equals_key": "yes"},
      "p3": {"used": "(2+ν)/(2-ν) = 1+2γ", "equals_key": "yes"},
      "b2": {"used": "2(1-ν)/(2-ν) = 1-γ, from γ plus a kinematic rotation of 1", "equals_key": "no"},
      "b3": {"used": "2/(2-ν) = 1+γ, from γ plus a kinematic rotation of 1", "equals_key": "no"},
      "m3": {"used": "1/(1-ν), via G=(1-ν²)(K1²+K2²)/E+(1+ν)K3²/E", "equals_key": "yes"}
    }
  },
  "claimed_checks": {
    "answer": "some_unsupported",
    "checks": [
      {"claim": "m_stiffness.py: closed forms match the 3×3 sector solves, symmetry, isotropy, Navier residual and the 2D limits, all to ≤1e-15", "claimed_at": "submission/method.md:66", "status": "supported", "evidence": "submission/m_stiffness.py:104-134 carries out every one of these comparisons; no saved output."},
      {"claim": "check_J.py and check_I2.py: finite-part closed forms (Ei) to 12 digits", "claimed_at": "submission/method.md:67", "status": "supported", "evidence": "submission/check_J.py:6-11 and submission/check_I2.py:9-24; no saved output."},
      {"claim": "check_hadamard.py, check_arbiter.py, check_plus.py, check_b1ahead.py: X-domain kernels agree with contour quadrature to ~1e-8 ahead and b1 closed forms to 1e-12", "claimed_at": "submission/method.md:68-70", "status": "cannot_tell", "evidence": "The four scripts exist and implement the comparisons, but no output is saved, so the stated tolerances cannot be confirmed."},
      {"claim": "check_b1on.py: numerical on-crack traction of the shear b1 columns reproduces the residue closed form", "claimed_at": "submission/method.md:71-72", "status": "supported", "evidence": "submission/check_b1on.py:8-22; no saved output."},
      {"claim": "A12 of X^{1/2}e^{-X} equals −A12 of the b1 column for xx and zz to 1e-5; A32 = 1/(√(2π)D_i) exactly and off-diagonal A32 is 0", "claimed_at": "submission/method.md:73-74", "status": "supported", "evidence": "submission/check_identity.py:13-19, check_identity2.py:10-15, first_order.py:56-60 and 137-139; the 1e-5 figure is not in a saved output."},
      {"claim": "√x-series fit of the total ahead traction and the algebraic large-p tail of MΦ̂ agree to 1e-5; the slip-route extraction was rejected by an exact inclined-front anti-plane test", "claimed_at": "submission/method.md:75-79", "status": "unsupported", "evidence": "check_direct.py, check_tail.py, check_tail2.py, check_tail3.py and check_sliproute.py compare the two routes, but no file contains the inclined-front anti-plane test cited as the reason for rejection (also claimed at method.md:59-62 and solve_AB.py:10)."},
      {"claim": "check_nu.py: the closed form R = −½I + (ν/(2−ν))[[1,i],[i,−1]] is reproduced at ν = 0.2, 0.25, 0.35 to ~5e-6", "claimed_at": "submission/method.md:80-81", "status": "cannot_tell", "evidence": "submission/check_nu.py:6-19 reruns the pipeline against that formula; no output is saved, so the tolerance cannot be confirmed."},
      {"claim": "Basis size and collocation convergence (ML = 4…8, Ncol, X-range): R stable to ~1e-5", "claimed_at": "submission/method.md:82", "status": "cannot_tell", "evidence": "first_order.py:160, check_direct.py:11 and check_sliproute.py:13 loop over ML, but no script varies Ncol or the X-range and no output is saved."},
      {"claim": "Part B: G(z) = Γ(z) holds pointwise to machine precision, f is real and zero-mean, and the Fourier series (n ≤ 60) has converged", "claimed_at": "submission/method.md:83-84", "status": "supported", "evidence": "submission/solve_AB.py:59-69, 79, 104-109."}
    ]
  },
  "error_location": {
    "answer": "theory",
    "first_wrong_step": "The shear response matrix is right, but the rotation from the global (x,z) basis into the local (n,t) basis is written with the mode-2 and mode-3 weights swapped: the solver subtracts a'K3⁰ from K2 and adds a'K2⁰ to K3, each with weight 1, instead of the weights 2/(2−ν) and 2(1−ν)/(2−ν). That gives b2 = 1−γ = 14/17 and b3 = 1+γ = 20/17, interchanged relative to the key, so the slope terms no longer cancel in G.",
    "where": "submission/solve_AB.py:8-9 and 31-32; also submission/method.md:59-60 and 88-91",
    "explains_groups": ["Part A: cosine-front coefficients", "Part B: loading level and front amplitude", "Part B: front shape f(s)", "Part B: kappa1(s)", "Part B: kappa2(s)", "Part B: kappa3(s)", "Part B: maximum of kappa3"]
  },
  "notes": "Only S2, S3 and SG fail in Part A; C1, S1, C2, C3, CG and K10_scaled pass, which is exactly what swapped b2 and b3 predict (p1, p2, p3 and m3 are right, and m3·b3 = b2 no longer holds, so G gains an odd part). The nonlocal weights are not taken from Gao & Rice: mode 1 is re-derived by Wiener–Hopf and the shear matrix comes from the solver's own least-squares solve in first_order.py, checked across Poisson ratios in check_nu.py. The 1989 citation is only for the reproduced mode-1 result. The inclined-front test that is supposed to justify the rotation is described but never written down, so that one check is unsupported; three other checks have code but no saved output, hence cannot_tell. Nothing in the folder was left out that would change this: the omitted files are large binaries."
}
```
