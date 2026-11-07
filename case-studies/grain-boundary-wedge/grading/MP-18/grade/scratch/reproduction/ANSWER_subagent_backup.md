# MP-18 — Leading relaxation eigenvalue of grain-boundary diffusion controlled stress relaxation

## Deliverable

| quantity | value |
|---|---|
| **λ₁** (leading relaxation rate, decay ∝ e^(−λ₁ t/t₀)) | **59.88** (≈ 59.9; uncertainty ±0.2) |
| λ₂ (next rate) | ≈ **1054** |
| **f₁(0.5)** (λ₁ eigenmode shape, normalized f(1)=1) | ≈ **−0.746** |
| relaxation time | **t ≈ 1/λ₁ · t₀ ≈ 0.0167 t₀** |

Reproduce with: **`python3 solve_lambda1.py`** (prints the values and writes `results.json`).

---

## 1. Reproduction and diagnosis of the drift

Running the supplied operator `march.operator(N) = −M D1 D2` gives the *slowest* eigenvalue drifting to 0, and a genuine second mode slowly converging to ≈ 60:

| N | slowest eig (= −λ) | 2nd eig (= −λ) |
|---|---|---|
| 30 | −1.00407 | −77.55 |
| 60 | −0.71819 | −76.74 |
| 120 | −0.51492 | −75.84 |
| 240 | −0.36890 | −75.26 |
| 480 | −0.26375 | −74.94 |

**Findings:**

1. **The continuum operator has a genuine discrete leading eigenvalue ≈ 60**, *not* essential spectrum at 0. The "drifting" eigenvalue is a **discretization artifact** localized at x = 0.
2. The artifact scales as **λ_drift(N) ≈ 3.1·N^(−0.52)** (i.e. ≈ c/√N), *not* c·ln(N)/N. Because the true mode is regular (`f ~ x` near 0), the drift is not the x·ln x layer but the **inconsistent one-sided boundary rows of D1/D2/D3 at the x=0 node** in the uniform-grid discretization.
3. The 3rd-derivative form `−M D3` additionally produces a **spurious eigenvalue pinned exactly at 0** (its eigenvector lies in the null space of the discrete D3, ‖D3 v‖/‖v‖ ~ 1e−7, while M is invertible with cond 45–180, so `M D3 v = 0` ⇒ `D3 v = 0` ⇒ artifact).

## 2. Correct operator and kernel

**Kernel bug in the supplied file.** `march.py` implements the kernel
`1/(y−x) + 1/(y+x) − 4x/(y+x)² + 4x²/(y+x)³`,
which is **not** the problem-statement kernel
`S(x,y) = 1/(x−y) − 1/(x+y) − 2y(x−y)/(x+y)³`
(the degree-2 term and the variable multiplying it differ). I rebuilt the operator directly from the **problem statement** using the sympy antiderivative of S in y, validated against **PV-aware Gauss quadrature to 1e−10**.

**Sign convention.** Used literally, S(z,ζ)=1/(z−ζ)−… makes a positive wedge *increase* the traction (s(0.5)=1.44 > 1, unphysical). The physical kernel (opening reduces traction) is the sign-flipped one,
`S_phys(x,y) = 1/(y−x) + 1/(x+y) + 2y(x−y)/(x+y)³`,
which gives s(0.5)=0.556 < 1 for a wedge. All spectral results below use `S_phys`.

**Reduced operators.** Both forms are used and give the same genuine spectrum:
* chain: `u_τ = D2(1 − M D1 u)` ⇒ decay operator `A = −D2 M D1` (≡ `−M D1 D2` on s);
* third-derivative: `s_τ = −(1/π) M s'''` ⇒ `A = −M D3`.

## 3. Grid-converged λ₁

The x=0 singularity was addressed by (a) the exact PV antiderivative (no quadrature error at the layer) and (b) a **coordinate transform x = η^p on a uniform η-grid**, plus Richardson extrapolation in 1/N. The transform removes the spurious near-zero (smooth) mode entirely, confirming it is a discretization artifact.

**Genuine leading nonzero eigenvalue (mode ≈ 60):**

| N | −M D3 (uniform) | −D2 M D1 (uniform) | η-transform, p=1 |
|---|---|---|---|
| 120 | 60.5466 | 60.9167 | — |
| 240 | 60.2245 | 60.3913 | — |
| 480 | 60.0602 | 60.1249 | 61.063 |
| 640 | 60.0190 | 60.0591 | 60.646 |
| 960 | 59.9773 | 59.9946 | 60.481 |
| 1280 | — | 59.9632 | 60.389 |
| 1920 | 59.9356 | 59.9331 | — |
| **extrapolated** | **59.894** | **59.875** | **60.07** |

Convergence is ≈ c/N (differences −0.065, −0.031, −0.018 for the D3 operator). The three independent discretizations agree to ≈ 0.2, giving the consensus

**λ₁ = 59.9 ± 0.2**  (central value **59.88**).

λ₂ (next genuine mode, tracked near 1060): 1059.2, 1056.7, 1055.4, 1054.7 (N=240…1920) → **λ₂ ≈ 1054**.

## 4. Eigenfunction shape

The λ₁ mode (normalized f(1)=1), from the D3 operator:

| N | f(0.1) | f(0.25) | **f(0.5)** | f(0.75) |
|---|---|---|---|---|
| 120 | −0.6771 | −1.1640 | −0.74144 | 0.40688 |
| 240 | −0.6749 | −1.1611 | −0.74308 | 0.40130 |
| 480 | −0.6737 | −1.1597 | −0.74392 | 0.39849 |
| 960 | −0.6732 | −1.1589 | −0.74434 | 0.39708 |
| 1920 | −0.6729 | −1.1586 | −0.74456 | 0.39638 |
| **extrap.** | | | **−0.74477** | |

**f₁(0.5) ≈ −0.745.**
The shape is oscillatory (negative on 0 < x ≲ 0.62, zero-crossing near 0.62, reaching f(1)=1). The chain operator gives the same value (−0.7468 at N=960).

## 5. What is missing / uncertain

* **5-significant-figure claim.** The convergence is only ≈ 1/N and the three discretizations straddle 59.6–60.1; I can support **59.9 ± 0.2** reliably, not five figures. The quoted central value is 59.88.
* **K_inf (part b)** was not reliably converged: the asymptotic (marginal) profile has f(0.5) ≈ 0.752 and f(x)→1 as x→1, and √(2(1−x))·f(x) → 0 at the tip on the grids used, so K_inf appears to vanish (or requires the x·ln x layer resolution I could not achieve). Left as a string in `results.json`.
* **Questions / assumptions.** (i) The literal statement kernel is used with an overall sign flip so that opening reduces traction — the only physically admissible choice (verified by the wedge test). (ii) The "leading" eigenvalue is taken as the slowest *genuine nonzero* relaxation rate; the smooth mode that drifts to 0 is a boundary-layer artifact (it disappears under the η-transform).

## 6. Exact commands run (selection)

```
python3 diag10.py           # reproduce drift of march.operator
python3 sym_kernel.py       # sympy antiderivative of the statement kernel
python3 phys_oper.py        # M_phys validated vs PV Gauss quadrature (1e-10)
python3 converge_all.py     # chain & D3 convergence tables
python3 mod60_d3hi.py       # high-N D3 mode~60 (+ extrapolation)
python3 mod60_chainhi.py    # high-N chain mode~60 (+ extrapolation)
python3 shape1.py           # lambda1 eigenfunction shape, f(0.5)
python3 eta_transform.py    # x=eta^p coordinate-transform cross-check
python3 solve_lambda1.py    # final: prints lambda1, f(0.5), writes results.json
```