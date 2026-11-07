# Verdict

Parameters: ν=0.2, ν_u=0.4, α=0.89. Length unit L, stress unit the face traction of loading 1. Report K in units of √(πL) and the small-time volume as V = V_undrained + A (L²/G) √(ct/L²) + o(√t).

Order of limits used for K(0+): **inner first**,

    K(0+)  :=  lim_{t→0+}  lim_{r→0+}  √(2π r) σ_yy(r, θ=0, t).

For any t>0 the r→0 limit is a drained r^{-1/2} field. The outer undrained intensity (r of order L, then t→0) is 1. These two numbers differ because a diffusion layer of thickness √(ct) sits at each tip.

Coupling that enters the inner problem is the **modulus ratio**

    η  = (ν_u − ν)/(1 − ν_u)  = 1/3,     1+η = (1−ν)/(1−ν_u) = 4/3.

This is independent of α. The paper’s η = α(1−2ν)/(2(1−ν)) = 0.33375 is a different object: the long-time loading-2 backstress K_{2,∞}.

---

## 1. K(0+)

**Value:** 0.80759 ± 0.00005 (in units √(πL)).

**Supports:** GPT (MP-37/38/39), not the paper’s 0.75.

### Decisive argument

The inner Laplace problem (length unit √(c/z), z→∞) is a Wiener–Hopf Riemann–Hilbert problem on the real line for a 2×2 Hermitian symbol N(ξ) with N(∞)=I and

    N(0) = diag(1+η, 1).

A plus/minus factorization N = N_-^* N_+ with N(∞)=I implies, at ξ=0,

    n_{11}^2 + n_{21}^2  =  N_{11}(0)  =  1+η.          (energy / J identity)

The drained-tip intensity in units of the outer undrained K is the (1,1) matching factor

    K(0+) / √(πL)  =  n_{11} / (1+η).

Independent numerical factorization of the same kernel (scratch/factor_independent.py; FFT projector, grids 4096 / 16384 / 65536 at scale 0.1, linear extrapolation in 1/count) gives

    n_{11} = 1.07679,   n_{21} = 0.41696,
    n_{11}^2 + n_{21}^2 − (1+η)  ∼  10^{-7},
    K(0+)  =  0.807593.

Copy-run of MP-37 factor_check.py at count=16384, scale=0.1: K(0+)=0.807780 (one-grid; their extrapolated 0.807593). Independent finite-crack Chebyshev–Fourier solve (scratch/spectral_check.py, n=16/24/32, q=√(s/c) up to 320) produces drained-COD intensities falling 0.848 → 0.817; quadratic/linear fit in 1/√q extrapolates to 0.8074–0.8075. Fits in 1/q do **not** go to 0.75.

Paper Eq. 29 is the substitution of the **undrained** elliptical COD into the **drained** Irwin formula,

    K_paper  =  (1−ν_u)/(1−ν)  =  1/(1+η)  =  0.75.

That is n_{11}=1, i.e. the inner WH factor is taken as the identity. It is a heuristic, as the paper itself says (p. 1650). It would be the answer only if there were no hydraulic column in N_-. Pure elastic energy matching (no Darcy dissipation, n_{21}=0) would give K=1/√(1+η)≈0.866. The true value lies between those two bounds, as required by n_{11}^2 + n_{21}^2 = 1+η with both n_{11}, n_{21} nonzero.

Weak-coupling expansion of the true WH result is K = 1 − (3/4)η + O(η²). At this particular η=1/3 that first-order number **coincides** with 1/(1+η)=0.75. The O(η²) remainder is +0.057, which is the whole discrepancy.

Rice–Simons (steadily moving permeable crack, high speed) is the same inner problem: diffusion length c/V → 0, drained tip inside an undrained outer K-field. Their K_tip/K_outer is a nontrivial function of speed, **not** the modulus ratio (1−ν_u)/(1−ν). The paper invoked Rice–Simons only to keep an r^{-1/2} tip, then replaced the whole COD by the undrained ellipse.

### Paper numerics cannot resolve this K

Six equal quadratic elements over 2L (p. 1655) ⇒ element length L/3. Fig. 9’s smallest times are t* = ct/L² ∼ 10^{-4}–10^{-3} (cf. Fig. 7 axis and the t*<0.001 statement for Fig. 6). Diffusion length √(ct)/L = √t* is 0.01–0.032, i.e. 10–33 times smaller than one element. Six Stehfest terms sample Laplace parameters n ln 2 / t; they do not create spatial resolution the mesh does not have. The code reports

    K  =  G √π  φ(L) / (2(1−ν)),

the COD-based intensity of the **outer** quadratic shape function. At small t* that φ(L) is the undrained ellipse, so the plotted K is **identically** the Eq. 29 heuristic. Fig. 9’s approach to 0.75 is therefore not an independent confirmation of a resolved inner tip.

### Check that could have failed

- Energy residual n_{11}^2+n_{21}^2−(1+η) remaining O(1) under grid refinement (it fell to 10^{-7}).
- Spectral K(q) extrapolating in 1/√q to 0.75 rather than 0.807, or depending on polynomial degree n.
- Using the paper’s η=0.33375 in the inner kernel (K would then depend on α; it does not).

### GPT code (K)

MP-37 factorization is the right inner problem. Minor notes, none of which move K at the 10^{-4} level: (i) their D(ξ)=1+η−2η k²/(m(m+k)) is the (1,1) entry after hydraulic elimination, not the naive mechanical factor 1+η(m−k)/(m+k); both share the same N(0) and N(∞), and the energy identity uses only those; (ii) extrapolation in 1/count is the leading projector bias. MP-38 spectral Schur system is the finite-crack version of the same kernel; reported 0.808 is a truncated-q rounding of 0.8076. MP-39 leaves K approximate 0.807 and is consistent.

Not settled to more than ±5×10^{-5}: a fully resolved inner collocation (scratch/tip_collocation.py) did not converge (truncation/Fourier). WH + spectral extrapolation already agree to that tolerance.

---

## 2. Small-time volume coefficient A

**Value:** 2.4402 ± 0.001.

**Supports:** GPT MP-37 and MP-38 (A≈2.44), **not** the paper and **not** MP-39 (A=1.354).

### Decisive argument

Reciprocity (paper Eq. 18 / loading-1 ↔ loading-2) gives

    V(t) − V_undrained  =  (const) × ∫_body ζ s^a dA

with s^a the elastic weight function. Two regions contribute at the **same** order √t:

1. **Face layers** (thickness √(ct) along |x|<L). ζ = O(1), area ∼ L √(ct). This is the paper’s calculation and MP-39’s 1.354:

       A_face  =  16 (ν_u−ν) / ((1+η) √π)
               =  16 (1−ν_u)(ν_u−ν) / ((1−ν) √π)
               =  1.354055.

2. **Tip cells** (two of them, size √(ct)). The inner stress-trace / fluid-content field Z is K-like, ∼ r^{-1/2}, so

       ∫ ζ s^a dA  ∼  (1/z) √L · √(c/z)   →   Laplace z^{-3/2}   →   √t,

the same order as the faces. The WH inner solution supplies the coefficient through n_{21}:

       A_tips  =  4 (ν_u−ν) π √(2/η) n_{21} / ((1+η) √π)  =  1.08616,
       A       =  A_face + A_tips
               =  4 (ν_u−ν) / ((1+η) √π)  [ 4 + π √(2/η) n_{21} ]
               =  2.44022.

Independent spectral volume (scratch/spectral_check.py) A_est(q) rises 1.83 → 2.27 as q=20→320; 1/√q extrapolation → 2.439, matching the WH formula, not 1.354. MP-37’s separate mixed-boundary FFT finite-crack check (finite_check.md; N up to 2^{21}, s up to 10^7) extrapolates A→2.44 from below.

The paper (and MP-39) computed only the face integral. That piece is correct as far as it goes; it is not the full √t coefficient. A_tips/A_face ≈ 0.80, so omitting tips is an O(1) error, not a higher-order correction.

### Paper mesh vs volume

Fig. 6: they report excellent small-time volume for t*<0.001. Face layers along the crack **are** resolved (the faces are the element surfaces; the 1-D diffusion into the body is in the Green’s function, not the crack mesh). Tip cells are not. A method that resolves faces and smears tips will converge to 1.354, which is what they plot. That is consistent with both their analysis and their mesh, and it does not test A_tips.

### Check that could have failed

- A_face failing to reproduce paper Eq. 18 algebraically (it matches to all printed digits).
- Spectral A(q) extrapolating to 1.354, or WH A depending on the FFT scale.
- Tip contribution scaling as t^{3/4} rather than √t (would happen if Z were bounded rather than r^{-1/2}).

### GPT code (volume)

MP-37: formula and finite_check are consistent; reported 2.44022 is the WH number above. MP-38: 2.44 from the same spectral u_0(q) used for K; coarser rounding. MP-39: agrees with the paper because it keeps only the face layer — internally consistent, incomplete.

---

## Summary

| quantity | paper | GPT MP-37/38 | GPT MP-39 | this review |
|---|---|---|---|---|
| K(0+)/√(πL) | 0.75 (heuristic + unresolved mesh) | 0.8076 / 0.808 | ≈0.807 | **0.80759 ± 0.00005** |
| A | 1.354 (faces only) | 2.44 | 1.354 (faces only) | **2.4402 ± 0.001** |

Uncontested values (K_{1,∞}=1, V_{1,0}=π(1−ν_u), V_{1,∞}=π(1−ν), p_{faces,0}=2(ν_u−ν)/[α(1−2ν)]≈0.74906, K_{2,0}=0, K_{2,∞}=−0.33375, V_{2,∞}=−πα(1−2ν)/2) are not in dispute.

If the problem had asked for the COD-based / outer intensity lim_{t→0+} lim_{r/√(ct)→∞, r/L→0} √(2πr) σ_yy, the paper’s 0.75 would be the answer to a **different** question. The stated definition is the true-tip limit, for which 0.75 is wrong.
