# Independent FE check of problem B: is K(0⁺) = 0.75 √(πL)? (log P365–P366)

**Why.** For loading 1, the key follows Detournay & Cheng 1991, Eq. 29: K(0⁺) = [(1 − ν_u)/(1 − ν)]
√(πL), which is 0.75 √(πL) for ν = 0.2 and ν_u = 0.4. The paper calls it a "heuristic argument"
(p. 1650).
- **GPT-6.1-sol challenged it in all three runs** (MP-37, 38, 39), giving 0.807–0.808 by a
  Wiener–Hopf-type analysis: a drained tip layer inside the undrained outer field.
- **Teng:** "can you check the solution of B yourself?" and "I have a feeling 0.75 has some
  approxmination".
- **This check** solves the stated problem directly, independently of both the paper's
  boundary-element method and GPT's factorization.

## Method (`fe_check_K0.py`; FEniCSx, dolfinx 0.10)

- **Model.** Biot poroelasticity exactly as in problem B:
  - total stress σ = 2Gε + λ tr ε I − αp I;
  - fluid content ζ = α div u + p/M, with M from λ_u = λ + α²M;
  - Darcy flow with κ from the paper's Eq. 8 so that c = 1. Checked: the standard formula
    κM(λ + 2G)/(λ_u + 2G) also gives 1.000000;
  - G = 1, ν = 0.2, ν_u = 0.4, α = 0.89.
- **Geometry.** A quarter domain [0, 20L]², L = 1, with the crack on y = 0, 0 ≤ x ≤ 1. The crack
  face carries traction σ_yy = −1, zero shear and p = 0 (loading 1). Ahead of the tip: u_y = 0,
  zero shear, no flux. x = 0: u_x = 0, no flux. Far sides: traction-free, p = 0.
- **Discretization.**
  - Taylor–Hood elements (u P2, p P1).
  - Backward Euler from rest. The first step approximates the undrained t = 0⁺ response; steps grow
    geometrically.
  - The mesh is refined NREF times in shrinking disks around the tip (radius 2 × 0.5^k). NREF = 16
    gives h_min = 1.1 × 10⁻⁵ (5,226 cells); NREF = 18 gives 2.7 × 10⁻⁶ (5,482 cells).
- **Extraction from the crack-face opening u_y(r),** r = distance from the tip:
  - **inner, drained:** K = u_y G √(2π)/[2(1 − ν)√r] for 5h_min < r < 0.05√(ct);
  - **outer, undrained:** the same with ν_u, for 5√(ct) < r < 0.03;
  - the reported ratio is the median K_inner/K_outer at the same time.

## Results (all values in units of √(πL))

| Run | h_min | time steps | resolved from t = | K_inner | K_outer | ratio (earliest steps) |
|---|---|---|---|---:|---:|---:|
| NREF 16 | 1.1e-5 | growth 1.25 from 1e-7 | 2.6e-6 | 0.8129 | 1.0180 | 0.7978–0.7986 |
| NREF 16 | 1.1e-5 | growth 1.1 from 1e-8 | 2.5e-6 | 0.8128 | 1.0178 | 0.7980–0.7988 |
| NREF 18 | 2.7e-6 | growth 1.25 from 1e-7 | 2.3e-7 | 0.8097 | 1.0113 | 0.8003–0.8013 |
| NREF 14 | 4.3e-5 | growth 1.25 | — (the tip layer is not resolved before the run ends) | — | — | — |

- **The profile across the layer** (NREF 16, t = 1.1 × 10⁻⁵, √(ct) = 3.3 × 10⁻³):
  - the drained-relation K has a plateau of 0.818–0.819 at r/√(ct) ≈ 0.03–0.15;
  - it reads low (0.805) only within two elements of the tip node;
  - the undrained-relation K returns to about 1.0 at r/√(ct) ≈ 10–30.
- **Physically consistent:** the near-tip field is the drained one, inside a thin layer, while the
  outer field is undrained.
- **Converged in time:** the time-step refinement changes the ratio by less than 0.001.
- **Mesh refinement raises it slightly,** 0.7985 → 0.801, with earlier resolved times.
- **Finite domain:** K_outer slightly above 1 (1.011–1.018) comes from the finite domain (W = 20L)
  and the outer window. The ratio removes it.

## Conclusion

**K(0⁺)/√(πL) ≈ 0.80 ± 0.01** for the stated problem.
- **The paper's heuristic 0.75 is low by about 7%.** It is an approximation, as the paper itself
  says, and as Teng suspected.
- **GPT's 0.807 is within about 1%** of the FE value. GPT's challenge is supported.
- **Opus's reading:** a plausible reason for the paper's value is that its numerical solution (six
  quadratic elements over 2L, p. 1655) cannot resolve a tip layer of width √(ct) ≈ 0.01L at
  t* = 10⁻⁴. Its Fig. 9 curve therefore shows the heuristic value at small times, not the true
  limit. This is not verified.
- **Limits of this check:**
  - the quarter-domain model assumes symmetry and shear-free faces, as GPT assumed;
  - the extraction uses windowed medians;
  - the finest mesh resolves t ≥ 2 × 10⁻⁷, and the t → 0⁺ value is read from the earliest steps.

## Second contested number: the small-time volume coefficient (1.354 paper against 2.44 GPT)

**FE with only the tip refined** (NREF 16, time-step growth 1.1 from 10⁻⁸; volume
V = 4∫u_y dx over the quarter-model face):
- V at the first step = 1.920, against the paper's 0.6π = 1.885 (+1.9%, from the finite domain and
  the coarse face mesh away from the tip);
- fits of V = V₀ + A√(ct/L²):

| t range | A |
|---|---:|
| 10⁻⁷ – 10⁻⁶ | 1.072 |
| 10⁻⁶ – 10⁻⁵ | 1.108 |
| 3×10⁻⁶ – 3×10⁻⁵ | 1.152 |

- The local slope is 1.07–1.09 for t ≤ 2×10⁻⁶.

**Not reliable as a total.** Away from the tip, the mesh does not resolve the diffusion layer along
the faces (√(ct) = 3×10⁻⁴ – 5×10⁻³ against an element size of about 0.06–0.25 there). That layer
carries the flat-face part of the volume growth. The FE number therefore contains mainly the tip
contribution.

**Codex's review** (`review_B_K10/codex/VERDICT.md`, below) splits A into 1.3541 from the flat faces
and 1.0862 from the tips. The FE's 1.07–1.09 is consistent with Codex's tip part. That is
consistency, not proof.

**Two failed attempts to resolve the face layer, both Opus's errors:**
- Refining a strip along the face down to 10⁻³ made the direct LU too heavy (11 GB, no step after
  12 min).
- A lighter strip re-refined the already-fine tip cells at every level (the strip included the tip
  region), producing 22.7 million cells and h_min = 8×10⁻⁸.

Both runs were stopped (logs kept as `*.attempt*`). The model is 2D plane strain throughout (Teng:
"Are you doing 3D simulations / 22 million is too large"). The fix is size-field refinement: split
an edge only if it is longer than the target size at its location, a few ×10⁴ cells.

**FE with the face layer resolved** (`fe_check_volume.py`, P367):
- **Mesh:** a size field, 2.5×10⁻³ within 0.05 of the face, graded to 1.1×10⁻⁵ at the tip; 39,926
  cells, 2D.
- **Time stepping:** growth 1.15 from 10⁻⁶ to 1.5×10⁻³.
- **V at the first step:** 1.8935, against the paper's 1.885 (+0.5%).
- **Local slope dV/d√(ct/L²):**

| t | 1.1×10⁻⁵ | 1.0×10⁻⁴ | 2.8×10⁻⁴ | 1.0×10⁻³ |
|---|---|---|---|---|
| slope | 2.21 | 2.02 | 1.92 | 1.75 |

- **Fits over a decade:**

| t range | A |
|---|---:|
| 10⁻⁵ – 10⁻⁴ | 2.11 |
| 3×10⁻⁵ – 3×10⁻⁴ | 2.00 |
| 10⁻⁴ – 10⁻³ | 1.88 |

- **Reading:**
  - The slope rises as t → 0. The decrease at later times is the O(ct/L²) term.
  - At the earliest time the face layer is only about one element thick, so that point is still
    under-resolved.
  - **Extrapolated, A ≈ 2.2–2.4: far above the paper's 1.354, consistent with Codex's 2.4402.**
    The precision is not enough to confirm the digits.
- **K(0⁺) on this mesh:** 0.8137 at the earliest step, the same as with the tip-only mesh (0.8128).
  Resolving the face layer does not change the tip quantity.

## The cross-family review (Codex, `review_B_K10/codex/VERDICT.md`, exit 0 at 16:30)

**Its own derivation from PROBLEM.md:**
- a Laplace–Fourier reduction to exact surface relations;
- a Wiener–Hopf factorization of the tip problem, solved by more than one method;
- a separate finite-crack calculation to check the model.

**Results:**
- **k₀ = K(0⁺)/√(πL) = 0.8075930** (numerical uncertainty about ±2×10⁻⁶). This supports GPT.
  - "The drained inner tip opening has an amplitude 1.07679067 times its matched undrained outer
    amplitude. The paper's Eq. (29) effectively sets this matching factor to one."
  - "The reversed stress limit, taking time to zero at a fixed positive distance first, is
    exactly 1. The paper's 0.75 is neither of these two stress limits."
- **A = 2.4402196 = 1.3540550 (flat faces) + 1.0861646 (two tips).**
  - This supports GPT's MP-37 and MP-38; MP-39 and the paper miss the tip part.
  - The tip part "survives at order √τ, even though the tip region shrinks".
- **Grok's review** (`review_B_K10/grok/VERDICT.md`, P372): K(0⁺) = 0.80759 ± 0.00005, A = 2.4402 ± 0.001. It supports GPT.

**Where this leaves B:**
- **Three independent routes now agree on K(0⁺) ≈ 0.807:**
  - GPT (three runs, Wiener–Hopf);
  - Codex (an independent derivation);
  - this FE (0.80 ± 0.01, direct simulation, slightly low as the windowed extraction suggests).
- **The paper's Eq. 29 is a heuristic** that is about 7% low. Its small-time volume Eq. 18 misses a
  tip contribution of the same order: FE A ≈ 2.2–2.4 against 1.354; Codex 2.4402.
- The key of B must be corrected before grading (section "Key correction" in `KEY.md`, after Grok's
  verdict).

## Gel parameters: reproducing Yu et al. 2020, Eq. 12 (P371)

- **Same script, same mesh** (NREF 18), at ν = 0.24, ν_u = 0.5 (1/M = 0), α = 1. It is selected by
  environment variables: `NU=0.24 NUU=0.5 ALPHA=1 TAGX=_gel`. The defaults are still problem B.
- **Ratio K_inner/K_outer at the earliest steps: 0.7256–0.7264.** GPT's factorization at the same
  values gives 0.7337, and Yu et al.'s Eq. 12 (α ≈ 0.735) gives 0.7351.
- **The FE is about 1% low here, as on B** (0.8003–0.8013 against 0.8076). The comparison is in
  `CASE_STUDY_GPT_CHALLENGE.md`, section "Reproduction of Yu et al. 2020".
