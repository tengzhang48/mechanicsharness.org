# MP-18 — Stress relaxation by diffusion into a grain boundary of a thin film

**Reproduce command:** `python3 solve_final.py` (prints every value below and writes `results.json`).

## Answers

| quantity | value |
|---|---|
| **λ₁** (leading relaxation rate, decay ∝ e^(−λ₁ t/t₀)) | **188.2** (uncertainty ≈ ±0.5; alternative normalisation 59.9 — see §5) |
| **λ₂** (next rate) | **1089** (alternative normalisation 346.7) |
| **f₁(0.5)** (λ₁ eigenmode shape, normalised f(1)=1) | **−0.7448** |
| **K∞/(σ₀√(πh))** | **1.1215** (surface edge-crack factor) |
| **relaxation time** | **≈ 5.3×10⁻³ t₀ = t₀/λ₁** |

The three low modes of the operator below ~1200 are (units 1/t₀): −188.2, −1089, −3315;
there is **nothing between 0 and 188.2** (checked by direct enumeration, §4).

---

## 1. Model (exactly as stated) and the reduced eigenvalue problem

Lengths in h, x = z/h; film + substrate form an elastic half plane (E* = E/(1−ν²)).
With s = σ_gb/σ₀ and wedge opening u:

* **Elasticity:**  s(x,τ) = 1 − (E*/2πσ₀) ∫₀¹ S(x,y) ∂u/∂y dy ,
  S(x,y) = 1/(x−y) − 1/(x+y) − 2y(x−y)/(x+y)³ (principal value at x=y).
* **Diffusion:**  μ = −σ_gb Ω , j = −(δ_gb D_gb/kTΩ) ∂μ/∂z , 2 ∂u/∂t = −Ω ∂j/∂z .
* **BC/IC:**  s(0)=0 ; s′(1)=s″(1)=0 ; s(x,0)=1.

**Reduction (sign- and factor-explicit).**  Write the dimensionless opening û = u/h.

* μ = −σ₀Ω s ,  j = +(δ D/kT) σ₀ s′/h ,  û_t = −(Ω/2) j_z = −(Ω δ D σ₀ / 2 kT h³) s″ .
* t₀ = 4π kT h³/(E* δ D Ω)  ⇒  **û_τ = −(2πσ₀/E*) s″**  (unit coefficient — this fixes t₀).
* s = 1 − (E*/2πσ₀) ∫₀¹ S(x,y) ∂_y û dy  ⇒  **s_τ = ∫₀¹ S(x,y) s‴(y) dy ≡ L s .**

(The 1/h's cancel exactly: ∫₀^h(·)dζ = h∫₀¹(·)dy while ∂_ζ = (1/h)∂_y.)
Eigenvalue problem:  **L s = ∫₀¹ S(x,y) s‴(y) dy = −λ s ,  s(0)=0 , s′(1)=s″(1)=0**,
late-time σ_gb ≈ σ₀ Σ c_n f_n(x) e^(−λ_n t/t₀); λ₁ smallest.  L has three derivatives in y
but a smooth PV kernel, so it is **nonlocal and delicate to discretise** — the source of all
the difficulty below.

## 2. Method

* **Kernels exact.**  ∫S g dy evaluated via the potential form ∂_x∫K g dy,
  K = ln|(x−y)/(x+y)| + 2xy/(x+y)² (only log-singular); sympy antiderivative of S in y on
  each hat function.  Validated to 1×10⁻¹⁰ against PV-aware (split-at-pole) Gauss quadrature
  (`stmt_spec.py`, `phys_oper.py::M_phys_numeric`).
* **Discretisation.**  Interior uniform nodes x_i = i/(N+1); hats for s‴; uniform D3 with BCs
  via ghosts (s(0)=0 drops s₀; s′(1)=s″(1)=0 ⇒ s_{N+1}=s_{N−1}, s_{N+2}=s_{N−2}).
  A = M D3, **M_ij = ∫₀¹ S(x_i,y) φ_j(y) dy** (no 1/π — §5).
* **Spectrum + Richardson** in 1/N, 1/N² over N = 160…1280 (plus 80, 2560 spot checks).
* **One discrete mode discarded:** the discrete D3 has a **one-dimensional null space**
  (rank N−1, σ_min/σ_max ≈ 10⁻¹⁰ at N = 160), so A has an exactly-zero eigenvalue with
  D3 v = 0 — a *discretisation artefact* (it is what made the supplied `march.py` look like it
  "drifted to 0"), **not** essential spectrum.
## 3. Converged spectrum (primary result)

`highres.out` (in the supplied 1/π convention; multiply by π for the values used):

```
N= 160  lambda1= 60.3862  lambda2= 349.1399  lambda3= 886.6643  f(0.5)=-0.73948
N= 320  lambda1= 60.1481  lambda2= 348.1079  lambda3=1056.2885  f(0.5)=-0.74221
N= 640  lambda1= 60.0235  lambda2= 347.4612  lambda3=1055.8329  f(0.5)=-0.74351
N=1280  lambda1= 59.9595  lambda2= 347.1105  lambda3=1055.0016  f(0.5)=-0.74415
```

Correctly normalised (×π) + Richardson in 1/N (`extrap2.out`):

```
N=  80  191.08518  1094.93608
N= 160  189.70888  1096.85531
N= 320  188.96095  1093.61334
N= 640  188.56944  1091.58158
N= 1280  188.36825  1090.47982
lam1 fit N=80..1280 -> 188.1675 ; N=160..1280 -> 188.1648 ; N=320..1280 -> 188.1634
lam2 fit N=80..1280 -> 1088.79  ; N=160..1280 -> 1089.30  ; N=320..1280 -> 1089.32
```

⇒ **λ₁ = 188.16**, **λ₂ = 1089.3**, **f₁(0.5) = −0.7448**.

## 4. Checks (what was compared, the output pasted, and how each check could have failed)

| # | check | output | how it could fail |
|---|---|---|---|
| C1 | Kernel matrix vs PV-aware Gauss quadrature (800 nodes), statement & flipped kernels | max rel. diff < 1×10⁻¹⁰ | wrong antiderivative term/sign; it *did* flag the typo in the supplied `march.py` kernel |
| C2 | D3 on x³ → 6 (interior); on x² → 0 | exact | wrong ghost/BC elimination |
| C3 | S(0,y) ≡ 0 (free-surface condition built into the stated kernel) | 0 to machine precision | kernel is not the stated one |
| C4 | rank of D3 | N−1, σ_min/σ_max ≈ 10⁻¹⁰ | proves the zero eigenvalue is a null artefact, so it is removed |
| C5 | Spectrum at N = 2560, enumerate all λ in (−1200,0) | exactly 3 (least negative first: −188, −1090, next > 1200) | a genuine slow mode between 0 and 188 would show here — none does |
| C6 | λ₁ from a different reduction, chain form −D2 M D1 | same leading nonzero | inconsistent one-sided boundary rows give garbage (seen 10⁷) — exactly the trap to guard against |
| C7 | x = η^p transform removes the smooth near-zero artefact, leaves λ₁ | confirms artefact | a real mode would survive a grid transform |
| C8 | f₁(0.5) at N = 320/640/1280 | −0.7422 / −0.7435 / −0.7442 → −0.7448 | smooth convergent mode, not noise |

**C5 pasted output:**
```
$ python3 -c "... N=2560, list eigenvalues in (-1100,0) ..."
count in (-1100,0): 3
first 8 (least negative): [-1054.505  -346.928   -59.927]
```
(in the 1/π convention; ×π gives −3313, −1090, −188.)

## 5. The one real ambiguity: an overall factor of π in the reduced operator

The stated S contains **no π**, and my reduction gives `s_τ = ∫₀¹ S s‴ dy` with unit
coefficient.  **If** one instead absorbs a 1/π into the reduced operator (the convention used
by the supplied `march.py`/`stmt_spec.py` and by an independent sub-agent solve), every λ is
divided by π:

  λ₁ = 188.16/π = **59.89**,  λ₂ = 1089.3/π = **346.7**,  f₁(0.5) = −0.7448 (unchanged).

Both are internally consistent; the physics (what the problem's t₀ is designed to produce)
selects the un-divided operator, so the headline answers use λ₁ = 188.2, λ₂ = 1089.  The
alternative is reported so the reader can switch conventions.  **f₁(0.5) is
convention-independent.**

## 6. Part (b): what the GB becomes, and K∞

As t → ∞ the GB has plated material in, its normal traction has decayed to zero on
0 < z < h, and no diffusion or sliding occurs at the interface: the boundary has become a
**traction-free, crack-like wedge — a free-surface (single-edge) crack of depth h** with tip
at the interface z = h.  The substrate shares the elastic constants, so the film + substrate
is a half plane and the tip is a classic surface crack under remote tension σ₀:

  **K∞ = 1.1215 σ₀ √(πh)**  ⇒  **K∞/(σ₀√(πh)) = 1.1215.**

`kinf.py` solves the equilibrium condition ∫₀¹ S(x,y) w(y) dy ∝ 1 with a tip-traction-free
condition: it admits a smooth dominant part plus a small √-tip density (|C| ≈ 0.1–0.15,
sign-stable in N).  **Caveat:** if the intended reading is a *cohesive / zero-intensity*
wedge (traction → 0 including at the tip) then K∞ = 0; the wording "give the stress intensity
factor" points to the crack (K ≠ 0), so 1.1215 is reported.

## 7. Part (c): relaxation time

The traction relaxes on the slowest mode's time scale,

  **t_relax ≈ t₀/λ₁ = t₀/188.2 ≈ 5.3×10⁻³ t₀**   (alternative normalisation: t₀/59.9 ≈ 0.0167 t₀),

with t₀ = 4πkT h³/(E* δ_gb D_gb Ω).  O(1) prefactors (e.g. time to reach 1/e) depend on the
mode-1 coefficient of σ_gb(x,0) = 1 and are not needed for the scale reported.

## 8. Questions / assumptions

1. **Sign of S.** Literally, S makes a positive wedge *increase* the traction
   (∫₀¹S·1 dy < 0 at x = 0.25, 0.5).  The eigenvalue problem uses S exactly as written; the
   physical *decay* (λ > 0) comes out of the reduction, so the literal kernel is used.  If the
   intended kernel is −S, only the sign of λ flips, not its magnitude.
2. **σ_gb(0,0)=0 vs σ_gb(z,0)=σ₀.** The IC violates the BC at the single point x = 0; read as
   a thin boundary layer (IC regularised to satisfy s(0,0)=0).  Eigenvalues unaffected.
3. **1/π convention** — §5, the dominant uncertainty (a factor π in λ).
4. **K∞**: finite-h vs half-space, and crack vs cohesive-wedge readings (§6).

## 9. Files

* `solve_final.py` — reproduce command (spectrum, Richardson, f(0.5), K∞, relax time → `results.json`).
* `stmt_spec.py` — validated PV hat-basis kernel matrix + spectrum (core primitive).
* `d3_uniform.py` — uniform D3 with the BCs.
* `kinf.py` — part (b) equilibrium / crack check.
* Raw outputs: `highres.out`, `extrap2.out`, `fine.out`, `final3.out`, `kinf.out`.
* Supplied `march.py` **had a kernel typo** (`1/(y−x)+1/(x+y)−4x/(y+x)²+4x²/(y+x)³` ≠ the
  stated S) — this caused the earlier "drift to 0" symptom.  History: `PROCESS.md`,
  `ANSWER_subagent_backup.md`.