# Visible messages and tool steps of glm_1, in order

Cut from the event stream of the run (opencode.jsonl, kept outside git). The stream holds the agent's messages and its
tool calls; it holds no private reasoning. Long tool inputs are shortened here.

- `read` (completed): `/work/PUBLIC_TASK.md`
- `read` (completed): `/work/PUBLIC_TASK.md`

## Message 1

The task is clear (I've read it fully). My plan: build the 2.5D crack-plane stiffness matrix analytically (sympy), then solve the scale-invariant first-order wavy-front problem per Fourier mode numerically, validate against the known Gao–Rice mode-I result, and assemble Parts A/B. Step 1: derive the Fourier-domain DD–traction matrix.

- `write` (completed): `/work/m_stiffness.py`
- `shell` (completed): `cd /work && python m_stiffness.py`
- `shell` (completed): `cd /work && (python3 --version; python3 -c "import sympy, numpy, scipy; print('sympy', sympy.__version__, 'numpy', numpy.__version__, 'scipy', scipy.__version__`
- `shell` (completed): `cd /work && python3 -c "import numpy, scipy; print('numpy', numpy.__version__, 'scipy', scipy.__version__)"; pip --version 2>/dev/null || pip3 --version 2>/dev/`

## Message 2

No sympy available — I'll do the algebra numerically in NumPy (exact linear algebra, no approximation): build the half-space solution, impose the sector mixed BCs, and verify my hand-derived closed forms to machine precision.

- `write` (completed): `/work/m_stiffness.py`
- `shell` (completed): `cd /work && python3 m_stiffness.py`

## Message 3

Two bugs in my *test* code (not the math): I passed an arbitrary `c` instead of `c_of_a` for the Navier check, and botched the 3×3 rotation matrix layout. Fixing the checks:

- `edit` (completed): `/work/m_stiffness.py`
- `shell` (completed): `cd /work && python3 m_stiffness.py`

## Message 4

M(p,q) is now verified exactly (closed forms match the exact sector-mixed solves to 1e-16; 2D limits, isotropy, Navier residual all pass). Step 2: the branch-cut contour machinery for the 2.5D kernels, with its constants pinned against known transform pairs.

- `write` (completed): `/work/kernels25d.py`

## Message 5

I made a folding error in the Gauss–Hermite helper. Let me rewrite the file cleanly with the correct quadrature plus the two exact verification tests:

- `write` (completed): `/work/kernels25d.py`
- `edit` (completed): `/work/kernels25d.py`
- `shell` (completed): `cd /work && python3 kernels25d.py`

## Message 6

The contour-constant bookkeeping needs to be settled definitively — my two reference tests implied inconsistent constants, so one of my cut-side value assignments has a bug. Let me write a diagnostic that checks cut values by direct near-cut sampling and the transform constants by direct real-line quadrature:

- `write` (completed): `/work/diag_constants.py`
- `shell` (completed): `cd /work && python3 diag_constants.py`

## Message 7

The two reference transforms force different contour constants, which is impossible for a pure Cauchy argument — something subtle (arc-at-infinity contributions) is being missed. Let me test the closed-contour identity numerically for the problematic case:

- `write` (completed): `/work/diag_contour.py`
- `shell` (completed): `cd /work && python3 diag_contour.py`

## Message 8

The contour-constant mystery is resolved (direct cut sampling confirms **C_upper = +i, C_lower = −i**; my earlier "i/√π" was a bookkeeping slip). However, the upper-cut integrals for the actual basis functions need Hadamard finite parts, so I'm switching to the cleaner equivalent: **X-domain collocation with closed-form Bessel kernels** (K₀/K₁) and classical singularity-subtraction (Ei-based closed forms for the finite-part pieces). Writing the kernel module:

- `write` (completed): `/work/kernels_x.py`
- `shell` (completed): `cd /work && python3 kernels_x.py`

## Message 9

Two bugs found: J has the wrong sign (J = +e^{−X}Ei(X), verified against the standard Ei definition), and C_xx's ν-term sign. Also cleanups. Fixing:

- `edit` (completed): `/work/kernels_x.py`
- `edit` (completed): `/work/kernels_x.py`
- `edit` (error): ``
- `edit` (completed): `/work/kernels_x.py`
- `edit` (error): ``
- `edit` (completed): `/work/kernels_x.py`
- `edit` (error): ``

## Message 10

The edit tool keeps dropping `path` on long payloads. Simpler to rewrite the file cleanly:

- `write` (completed): `/work/kernels_x.py`
- `shell` (completed): `cd /work && python3 kernels_x.py`

## Message 11

K_analytic is now smooth to 0 at t→0 for all entries. The J test was badly designed (asymmetric PV sampling); the closed form rests on the standard Ei definition, which I'll verify once more cleanly, then write the main first-order solver:

- `write` (completed): `/work/check_J.py`
- `shell` (completed): `cd /work && python3 check_J.py`

## Message 12

J is verified to 12 digits. Now the main first-order solver — null-space collocation, edge extraction, and all the built-in physical checks (Gao–Rice anchor, D-constants, slip-vs-traction route consistency):

- `write` (completed): `/work/first_order.py`
- `shell` (completed): `cd /work && timeout 1200 python3 first_order.py`

## Message 13

The adaptive-quadrature approach is drowning in cancellation noise and bad scaling. Restructuring: fixed-panel Gauss–Legendre quadrature (panels split at Xⱼ, deterministic nodes, no noise-chasing), proper row/column scaling, and — importantly — an independent cross-check of the finite-part machinery against the (already validated) contour method for the one case where the latter converges without finite parts:

- `write` (completed): `/work/kernels_x.py`
- `write` (completed): `/work/check_cross.py`

## Message 14

Before writing the full check, a striking analytic discovery: for b = X^{−1/2}e^{−X} in the y-sector, the two cut-side jumps cancel exactly (f_R − f_L = 0), suggesting zero on-crack traction. Let me verify numerically with the independent X-domain machinery — this affects the null-space structure:

- `write` (completed): `/work/check_cross.py`
- `shell` (completed): `cd /work && python3 check_cross.py`

## Message 15

The X-domain result contradicts my contour "cancellation" — correctly so: that integrand *grows* at infinity, where the naive arc argument fails; the Hadamard value from the X-domain machinery is the right one. The anchors (D-recovery, Gao–Rice) will be the decisive validation. Fixing the `a=0` geometric-panel bug and rewriting the solver with the vectorized API and proper scaling:

- `edit` (error): ``
- `edit` (completed): `/work/kernels_x.py`
- `write` (completed): `/work/first_order.py`
- `shell` (completed): `cd /work && timeout 1800 python3 first_order.py`

## Message 16

Key insight from the run: the **slip route gives exactly −1/2** across all basis sizes — but the traction route is garbage. The bug: I scale basis functions for conditioning but compute edge coefficients as raw coefficient sums, ignoring the scales. Quick verification of the fix and the ahead-evaluation magnitudes:

- `write` (completed): `/work/check_ahead.py`
- `shell` (completed): `cd /work && python3 check_ahead.py`

## Message 17

Ahead evaluation works correctly (trend → 0.357 as expected). Fixing the edge-coefficient normalization bug (divide by basis scales):

- `edit` (completed): `/work/kernels_x.py`
- `edit` (completed): `/work/first_order.py`
- `edit` (completed): `/work/first_order.py`
- `edit` (completed): `/work/first_order.py`
- `shell` (completed): `cd /work && timeout 1800 python3 first_order.py 2>&1 | tail -40`

## Message 18

The SVD reveals the real problem: the discretized null space is numerically 2-dimensional (two singular values at ~1e-16), so SVD returns arbitrary garbage combinations. Much better formulation: **solve the inhomogeneous problem directly** — prescribe the X^{−1/2}-edge input via a fixed function b₁ = X^{−1/2}e^{−X}, and least-squares solve for the √X-family remainder. Unique solution, no null-space ambiguity. Rewriting:

- `write` (completed): `/work/first_order.py`
- `shell` (completed): `cd /work && timeout 1800 python3 first_order.py 2>&1 | tail -60`

## Message 19

The LS residuals (~0.5) say the equations can't be satisfied — either the formulation or the kernel application is wrong. Time for a decisive test of the Hadamard machinery: for basis b = X^{3/2}e^{−X} the p-domain product *decays*, so the real-line Fourier integral converges and can be computed independently by oscillatory quadrature:

- `write` (completed): `/work/check_hadamard.py`
- `shell` (completed): `cd /work && python3 check_hadamard.py`

## Message 20

The ahead values disagree between two independent methods — one must have a bug. The arbiter: a third evaluation via the (already validated) lower-cut contour formula, which for this decaying integrand is a plain smooth integral:

- `write` (completed): `/work/check_arbiter.py`
- `shell` (completed): `cd /work && python3 check_arbiter.py`

## Message 21

The X-domain machinery is actually **validated** — the contour script had a cut-side sign slip (corrected: contour = −(my script's value) = X-domain value exactly). QAWO is simply unreliable for these integrals. So the kernels are right; the failure is in the solve. Diagnosing the LS residual structure directly:

- `write` (completed): `/work/diag_ls.py`
- `shell` (completed): `cd /work && python3 diag_ls.py`

## Message 22

Major breakthrough — the Wiener–Hopf structure cracks the problem analytically:

- **y-sector solved exactly**: b = c₋·e^{−X}X^{−1/2}, traction-free on the crack (Abel prescription), ahead traction = [1/(4(1−ν))]·e^{−|x|}|x|^{−3/2}, giving **r_y = −1/2 exactly** (Gao–Rice reproduced analytically) and exact D-consistency (0.357142857 = 1/(√(2π)D₁)).
- My Hadamard machinery used the wrong prescription for the X^{−1/2}-edge function (where M·Φ̂ grows); all other columns are fine.
- Shear sector: the b₁-columns have closed-form Abel values (e^{−X}·residue matrix on-crack; erfc/Γ forms ahead), with **exact** A₃₂ = −Q̂/2 → D_x, D_z exact, off-diagonal A₃₂ ≡ 0.

Rewriting the solver with these exact b₁-columns:

- `write` (completed): `/work/first_order.py`
- `shell` (completed): `cd /work && timeout 1800 python3 first_order.py 2>&1 | tail -45`

## Message 23

Two bugs found by re-deriving the cut-side algebra: P̂_xx = **+ccν** (not −ccν), and the on-crack residue matrix is missing the e^{−iπ/4} phase from i^{−1/2} (correct prefactor is the *real* √π·ccν/√2). Fixing and adding a direct-quadrature validation of the ahead closed forms:

- `edit` (completed): `/work/first_order.py`
- `edit` (completed): `/work/first_order.py`
- `write` (completed): `/work/check_b1ahead.py`
- `shell` (completed): `cd /work && python3 check_b1ahead.py`
- `edit` (completed): `/work/check_b1ahead.py`
- `shell` (completed): `cd /work && python3 check_b1ahead.py`
- `edit` (completed): `/work/check_b1ahead.py`
- `shell` (completed): `cd /work && python3 check_b1ahead.py`

## Message 24

The zz decomposition has a sign slip (1/(iρ) = −i/ρ was mishandled): P̂_zz = **−**ccν. The x=0.01 mismatches are just quadrature truncation (tail ~34 ≈ observed diff). Fixing:

- `edit` (error): `/work/first_order.py`
- `edit` (completed): `/work/first_order.py`
- `edit` (completed): `/work/check_b1ahead.py`
- `shell` (completed): `cd /work && python3 check_b1ahead.py && echo ===== && timeout 1800 python3 first_order.py 2>&1 | tail -32`

## Message 25

Structure is now clean (off-diagonals pure imaginary ✓), but the ~10% LS residual and ML-drift remain. Two suspects: (a) the m=0 plus-function (ν=3/2, conditionally-decaying product) may also need Abel correction; (b) residual rhs error. Diagnosing both:

- `write` (completed): `/work/check_plus.py`
- `shell` (completed): `cd /work && python3 check_plus.py`

## Message 26

Critical discovery: the *ahead* values are all validated (1e-8 agreement with contour), but the **on-crack** Hadamard values were never validated — and they look orders of magnitude too large. That explains the drifting LS. Testing the I₂ machinery on the one case with a closed form (b = e^{−X} → fp-integral = J − 1/Xⱼ):

- `write` (completed): `/work/check_I2.py`
- `write` (completed): `/work/check_I2.py`
- `shell` (completed): `cd /work && python3 check_I2.py`

## Message 27

Found the core bug: my I₁/I₂ subtraction models `b(Xj)e^{−X}` don't vanish at X=Xj, so the "regularized" integrands retain the pole. The correct model is `φ = e^{−(X−Xj)}` (equals 1 at Xj). Fixing kernels_x with proper closed forms (A₂ = e^{Xj}(J−1/Xj), etc.):

- `edit` (error): `/work/kernels_x.py`

## Message 28

Given the repeated edit-matching friction, rewriting the file whole:

- `write` (completed): `/work/kernels_x.py`
- `write` (completed): `/work/check_I2.py`
- `shell` (completed): `cd /work && python3 check_I2.py && python3 check_cross.py 2>/dev/null | head -8`

## Message 29

I₁ is now exact; I₂'s subtraction model must match b's value AND first two derivatives at Xj. Implementing proper Hermite-type models φ₀=e^{−u}(1+u+u²/2), φ₁=e^{−u}(u+u²), φ₂=e^{−u}u²/2 with closed-form finite parts:

- `edit` (completed): `/work/kernels_x.py`
- `write` (completed): `/work/check_I2.py`
- `shell` (completed): `cd /work && python3 check_I2.py && python3 check_cross.py 2>/dev/null | head -14`

## Message 30

I₁/I₂ now exact to 10 digits, and the b₁ y-column on-crack traction is ~10⁻⁴–10⁻³ smaller than the ahead values — confirming the WH prediction (exactly zero, up to quadrature noise). Now check whether the fixed machinery reproduces the shear residue closed forms:

- `write` (completed): `/work/check_b1on.py`
- `shell` (completed): `cd /work && python3 check_b1on.py`

## Message 31

The fixed machinery now reproduces the shear residue closed forms to ~4 digits (deviations only where fp-cancellation noise is worst). The closed forms and the machinery mutually validate. Running the full solve:

- `shell` (completed): `cd /work && timeout 1800 python3 first_order.py 2>&1 | tail -36`

## Message 32

LS residual is now 1.7e-4 and results are ML-converged. The response matrix R = [[−0.3699, −0.1859i],[−0.0876i, −0.6252]] — non-universal diagonals are actually *expected* theoretically (the ν/r terms add x^{−1/2}-only contributions to the ahead field, modifying A₁₂/(2A₃₂) away from −1/2 in the shear sector). Now the crucial independent validation: the slip route (completely independent of the ahead fits):

- `write` (completed): `/work/check_sliproute.py`
- `shell` (completed): `cd /work && timeout 1800 python3 check_sliproute.py`

## Message 33

The two routes disagree systematically — and I see why: my 3-term ahead fit omits the **constant** terms that the erfc-class closed forms produce, biasing A₁₂. The clean fix: the expansion is a regular power series in u=√x, so evaluate at very small x and fit f(u) = A₁₂ + A₀u + ... . Fixing:

- `edit` (completed): `/work/first_order.py`
- `shell` (completed): `cd /work && timeout 1800 python3 check_sliproute.py 2>&1 | tail -20`

## Message 34

The very-small-x evaluation fails because the ahead quadrature panels don't resolve down to |x|~1e-6. Making the ahead panels x-adaptive and simplifying the fit:

- `edit` (completed): `/work/kernels_x.py`
- `write` (completed): `/work/check_fit.py`
- `shell` (completed): `cd /work && python3 check_fit.py`

## Message 35

The √x-extrapolation is now clean (A₁₂ = 0.357138, A₃₂ → −5e-9). Re-running both extraction routes:

- `shell` (completed): `cd /work && timeout 1800 python3 check_sliproute.py 2>&1 | tail -21`

## Message 36

My weighted-lstsq formulation was structurally wrong (treated T as polynomial in u). Correcting: T·u² = A₁₂ + A₀u² + ... and T·u⁶ = A₃₂ + A₁₂u⁴ + ...:

- `edit` (completed): `/work/first_order.py`
- `shell` (completed): `cd /work && timeout 1800 python3 check_sliproute.py 2>&1 | tail -21`

## Message 37

Both routes can't be right. Key observation: the traction-route coupling ratio R12/R21 = 1.4279 matches D₃/D₁ exactly (clean structure), suggesting traction is right. Testing the y-sector identity A₁₂[φ₀⁺] = −A₁₂[b₁] per kernel entry — if it holds everywhere, the routes must agree and there's an arithmetic bug:

- `write` (completed): `/work/check_identity.py`
- `shell` (completed): `cd /work && python3 check_identity.py`

## Message 38

The identity **holds** (A₁₂[φ₀⁺] = −A₁₂[b₁] per entry), which means the two routes must agree for any solution — my last "fix" used a wrong power model (T·u not T·u²), and the b₁-xz A₁₂ (=+0.2143i) was missing. Correcting all three pieces:

- `edit` (completed): `/work/first_order.py`
- `edit` (completed): `/work/first_order.py`
- `edit` (completed): `/work/first_order.py`
- `shell` (completed): `cd /work && timeout 1800 python3 check_sliproute.py 2>&1 | tail -22`

## Message 39

Diagonals now agree to 3e-5 ✓. The off-diagonal assembly is missing the cross-kernel reg-terms (T_x gets K_xx·b_x **and** K_xz·b_z). Fixing the extraction assembly:

- `edit` (completed): `/work/first_order.py`

## Message 40

I wrote a mangled block — rewriting it cleanly:

- `edit` (error): `/work/first_order.py`
- `read` (completed): `/work/first_order.py`
- `edit` (error): `/work/first_order.py`
- `write` (completed): `/work/first_order.py`
- `shell` (completed): `cd /work && timeout 1800 python3 first_order.py 2>&1 | tail -40`

## Message 41

A₃₂/D are now exact to all digits and diagonals agree between routes (1e-5). The off-diagonal discrepancy must lie in the identity for the xz kernel or in the rotation cross-terms. Testing the identity for xz specifically:

- `write` (completed): `/work/check_identity2.py`
- `shell` (completed): `cd /work && python3 check_identity2.py`

## Message 42

The identity fails only for xz (A₁₂[φ₀⁺,xz] = 0 exactly), so the routes genuinely differ at the coupling — one has a sign error. To settle it without hand-phase algebra, I'll extract A₁₂ **algebraically from the large-p tail of M·Φ̂** (Watson pairing), first calibrating the tail-formula against the known b₁ closed forms:

- `write` (completed): `/work/check_tail.py`
- `write` (completed): `/work/check_tail.py`
- `shell` (completed): `cd /work && python3 check_tail.py`

## Message 43

Calibration constants are entry-independent (pure Watson phases) ✓. Now the fully-algebraic traction route applied to the solved coefficients:

- `write` (completed): `/work/check_tail2.py`
- `shell` (completed): `cd /work && timeout 1800 python3 check_tail2.py 2>&1 | tail -24`

## Message 44

A third answer! The routes disagree because of extraction-layer issues. Cleanest possible arbiter: evaluate the **total** ahead traction T(x) from validated components only (b₁ closed forms + validated ahead machinery), subtract the exactly-known A₃₂x^{−3/2} term, and Richardson-extrapolate the remainder in u=√x:

- `write` (completed): `/work/check_direct.py`
- `shell` (completed): `cd /work && timeout 1800 python3 check_direct.py 2>&1 | tail -24`

## Message 45

Major resolution: the exact inclined-plane test **falsifies the slip route** — displacement-jump asymptotics don't pair with the traction-defined K at O(a′) in coupled 3D situations (the u_y,n-term essential to the D₂-pairing is absent). The traction route is correct. Now to arbitrate its two variants, testing the p-tail algebra on the known φ₀⁺ value:

- `write` (completed): `/work/check_tail3.py`
- `shell` (completed): `cd /work && python3 check_tail3.py`

## Message 46

The algebraic route passed all calibration tests — but my Φ̂_sol assembly in check_tail2 ignored the basis **scales** and the multi-term Laguerre structure. Fixing:

- `write` (completed): `/work/check_tail2.py`
- `shell` (completed): `cd /work && timeout 1800 python3 check_tail2.py 2>&1 | tail -14`
- `edit` (completed): `/work/check_tail2.py`
- `shell` (completed): `cd /work && timeout 1800 python3 check_tail2.py 2>&1 | tail -12`

## Message 47

**Breakthrough — exact structure found.** The corrected algebraic route now matches the x-fit route perfectly:

**R = −½I + γ·[[1, i],[i, −1]]**, with γ = 0.17647 = **3/17 = ν/(2−ν)** exactly (numerics match to 1e-5). The slip route was falsified by the exact inclined-plane test. Strong validation: re-run at other ν to confirm γ = ν/(2−ν):

- `write` (completed): `/work/check_nu.py`
- `shell` (completed): `cd /work && timeout 1800 python3 check_nu.py 2>&1 | tail -24`

## Message 48

**Exact operator confirmed across three ν values** (agreement to ~5e-6):

**R = −½I + γ[[1, i],[i, −1]] with γ = ν/(2−ν)** — plus the exact Gao–Rice mode-I diagonal. Now assembling the final answer (Part A single mode + Part B per-mode equilibrium) with the basis-rotation conversion to the problem's local K definitions:

- `write` (completed): `/work/solve_AB.py`
- `shell` (completed): `cd /work && python3 solve_AB.py`

## Message 49

The fine-grid series evaluation is too slow in pure Python. Vectorizing:

- `edit` (completed): `/work/solve_AB.py`
- `edit` (completed): `/work/solve_AB.py`
- `shell` (completed): `cd /work && timeout 600 python3 solve_AB.py`

## Message 50

All values computed — just a JSON serialization fix needed, plus a pointwise equilibrium check (G(z) must equal Γ(z)):

- `edit` (completed): `/work/solve_AB.py`
- `shell` (completed): `cd /work && timeout 600 python3 solve_AB.py`

## Message 51

Equilibrium holds to machine precision. Verifying the sharp κ₃ peak is genuine with a targeted fine evaluation:

- `shell` (completed): `cd /work && python3 -c " import numpy as np, json exec(open('solve_AB.py').read().split('# ---------------- Part B')[0]) q = 0.5 nvals = np.arange(1, 61) fhat =`

## Message 52

The peak is genuine and smooth (max ≈ 1.4650 at s ≈ 0.0324). Checking the final results file:

- `shell` (completed): `cd /work && cat results.json`

## Message 53

Everything is verified. Writing the method documentation:

- `write` (completed): `/work/method.md`

## Message 54

Final end-to-end reproducibility run of the production chain:

- `shell` (completed): `cd /work && rm -rf __pycache__ && python3 m_stiffness.py > /dev/null && echo "m_stiffness OK" && timeout 1200 python3 first_order.py 2>&1 | tail -14 && timeout `

## Message 55

Clean reproduction. Final tidy-up and folder inventory:

- `shell` (completed): `cd /work && rm -f tmp_y.npy first_order_results.json && ls -la && python3 -c "import json; d=json.load(open('results.json')); print('results.json OK:', list(d.k`
- `shell` (completed): `cd /work && rm -rf __pycache__ kernels25d.py && ls`

## Message 56

Done. `results.json`, `method.md`, and all code are in the folder.

## Summary

**What I reached**

- **Exact first-order wavy-front operator** (per z-Fourier mode, ν = 0.3, γ = ν/(2−ν) = 3/17):
  - δK₁/K₁⁰ = −(|k|/2)δâ — the classical Gao–Rice result, reproduced **analytically** via Wiener–Hopf (the mode-I correction is exactly b = c₋e^{−X}X^{−1/2}, traction-free on the crack);
  - δK₂^global = |k|δâ[(−½+γ)K₂⁰ + iσγK₃⁰], δK₃^global = |k|δâ[(−½−γ)K₃⁰ + iσγK₂⁰], converted to the problem's local (n,t)-basis K definitions by the rotation K₂^prob = K₂^g − a′K₃⁰, K₃^prob = K₃^g + a′K₂⁰ (verified against the exactly-solvable inclined-front anti-plane problem).

- **Part A**: C₁ = −π, S₁ = 0, C₂ = −1.2197, S₂ = +4.1395, C₃ = −3.4003, S₃ = −4.4352, C_G = −6.8237, S_G = −2.2733.

- **Part B**: K10_scaled = 0.66310 (= [(1−ν²)(1.36)+(1+ν)0.64]^(−1/2)·(1−ν²)^(1/2), exact); per-mode algebraic solution of the trapping equilibrium (Poisson-kernel Γ profile); front_peak_to_peak = 0.30412; the 16-point tables for f, κ₁, κ₂, κ₃; κ₃max = 1.46499 at s = 0.03242. Equilibrium G(z)=Γ(z) verified pointwise to machine precision.

**Validation**: exact stiffness matrix (1e-16 checks, 2D limits, isotropy); closed-form Bessel kernels with finite-part application cross-checked against independent contour quadrature (1e-8–1e-12); exact D-consistency of the A₃₂ coefficients; two independent K-extraction routes agreeing to 1e-5 (a third, slip-based route was deliberately rejected using an exact test); the closed form R = −½I + (ν/(2−ν))[[1,i],[i,−1]] re-derived by the full pipeline at ν = 0.2/0.25/0.35 to ~5e-6.

**What remains unverified**

- No comparison against literature values for the mixed-mode (II+III) coupling coefficient γ = ν/(2−ν) — I had no access to references; it rests on the internal cross-validation above (the mode-I anchor is the one externally known result).
- The far field was taken as the pure straight-front K-field; possible non-singular far-field terms (T-stress etc.) are unspecified by the task and assumed absent — they don't enter the first-order K-perturbation for this loading class, but a different interpretation of "loaded from infinity" could in principle add O(A/λ) terms.
- Results are first-order accurate in A/λ and ε by construction; the reported numbers should be exact within that theory (numerical error ≲ 10⁻⁵), but second-order corrections were not assessed.

