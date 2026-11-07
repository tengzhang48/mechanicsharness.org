# Process log

- start: 16:04:00
- 16:04–16:20: read TASK.md, PROBLEM.md, KEY.md, paper text (Eq. 29, p. 1655 numerics, Fig. 6/9), MP-37/38/39 ANSWER.md + results.json
- 16:20–16:50: independent WH factorization scratch/factor_independent.py (grids 4096/16384/65536, scale 0.1) → n11=1.07679, n21=0.41696, K1_0=0.807593, A=2.44022, energy residual ~1e-7
- 16:50–17:05: independent spectral finite-crack scratch/spectral_check.py (n=16/24/32, q=20..320); 1/√q extrapolation K0≈0.8074–0.8075, A0≈2.439; 1/q fits do not go to 0.75
- 17:05: copy-run MP-37 factor_check.py at count=16384, scale=0.1 → K1_0=0.807780
- 17:05–17:10: scratch/tip_collocation.py failed (residuals O(1)); abandoned
- 17:12: scratch/identities.py — energy identity, weak-coupling K=1-(3/4)η, A_face=1.354055 ≡ paper Eq. 18, A_tips=1.08616, A_tot=2.44022; mesh √(ct) vs L/3 at t*=1e-4..0.1; Rice–Simons correspondence
- 17:13: scratch/kernel_check.py — D vs dmech share N(0)=diag(1+η,1) and N(∞)=I; energy identity uses only those
- 17:14: wrote VERDICT.md
- stop: 17:15:00
