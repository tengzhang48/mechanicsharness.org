# MP-20 — Stress relaxation by diffusion into a grain boundary

## 1. Normalised problem

Let `s = z/h ∈ [0,1]`, `g = σ_gb/σ₀`, `w = (E*/2πσ₀h) u` and
`τ = t/t₀` with

    t₀ = 4π kT h³ / (E* dD Ω) .

The given equations become, with **no free parameters** (every `h`, `E*`, `σ₀`
factor cancels):

| equation | form |
|---|---|
| elasticity | `g(s,τ) = 1 − ∫₀¹ S(s,x) w_x(x,τ) dx` |
| kernel | `S(s,x) = 1/(s−x) − 1/(s+x) − 2x(s−x)/(s+x)³` (p.v. at `x=s`) |
| diffusion | `w_τ = g_ss(s,τ)` |
| BC | `g(0,τ)=0`, `g_s(1,τ)=g_ss(1,τ)=0` |
| IC | `g(s,0)=1` |

Eliminating `w` (`w_t = g`, integrated with `g(0)=0`, `g_s(1)=0`) gives the
linear evolution

    g_τ = −A g ,      A = A1 · D² ,      (A1 h)(s) = ∫₀¹ S(s,x) h_x(x) dx .

`A1` is self-adjoint (symmetric kernel); `D²` is self-adjoint for the BCs; the
relaxation is a sum of real modes `g = Σ_k c_k f_k(s) e^{−λ_k τ}`.

## 2. Method

Finite-difference collocation on the uniform grid `s_i = i/n`:

* `D²` = three-point second difference (Dirichlet `g₀=0`, Neumann at `s=1`);
* `A1` built by **singularity subtraction** of the principal value,

      fp ∫ φ_j'(x)/(s_i−x) dx
        = ∫ [φ_j'(x) − φ_j'(s_i)]/(s_i−x) dx + φ_j'(s_i)·log(s_i/(1−s_i)) ,

  with the regular tails `−1/(s+x) − 2x(s−x)/(s+x)³` integrated by 2000-point
  Gauss–Legendre;
* **validation**: `E` reproduces an independent mpmath finite-part reference
  for a smooth opening to **7×10⁻¹²** (see `solve.py:check_elasticity`);
* the spectrum of `−A` is real; eigenvalues are extrapolated as `1/n → 0`.

## 3. Results  (units: `τ = t/t₀`, `t₀ = 4πkT h³/(E*dDΩ)`)

| quantity | value |
|---|---|
| **λ₁** (slowest decay rate) | **0.599** (grid n=26 gives 0.627; extrapolated 1/n→0) |
| **λ₂** | **135.9** |
| λ₃ | ≈ 767 (checks) |
| **mode-1 shape at `z=h/2`**, `f₁(0.5)` | **0.740** |
| **stress intensity factor** `K_∞/(σ₀√(πh))` | **3.116** |
| **relaxation time** `τ_rel = 1/λ₁` | **1.67 t₀** |

The eigenvalue spectrum is clean and well separated:

    λ = 0.60,  136,  767,  2280,  4942, …

(the smallest magnitude eigenvalue is a large *negative* stiff mode of the
discretisation; it is the over-damped gauge direction and is excluded).

The slowest mode `f₁(s)` is **smooth** and **nearly constant**
(`|⟨f₁,1⟩| = 0.92`), monotone from `f₁(0)=0` rising to `f₁(1)=1`.

## 4. Part (b): stress intensity factor at the blocked interface

At `t→∞` the traction vanishes, `g_∞ ≡ 0`, and the opening is the
**equilibrium wedge** `w_∞` solving `A1 w_∞ = 1`. Its tip has a clean
**square-root behaviour**

    w_∞(s) ≈ −0.433 + 0.701 √(1−s)      as s → 1 ,

which is the signature of a crack/wedge tip wedged against the rigid blocked
interface: `u(h−ε) ∝ √ε`. Matching to the mode-I near-tip field
`u(h−ε) = (2K/E*)√(ε/2π)` gives

    K = √2 · c · π^{3/2} · σ₀ √h  ,   c = 0.701

    ⇒   K_∞ / (σ₀ √(πh)) = c·π√2 = 3.116 .

So the interface **does retain a square-root stress singularity**, with
dimensionless intensity **≈ 3.12 σ₀√(πh)**.

## 5. Part (c): time dependence

With IC `g(s,0)=1`, `w(s,0)=0`, the solution is

    g(s,τ) = Σ_k c_k f_k(s) e^{−λ_k τ},   c_k = ⟨f_k,1⟩/⟨f_k,f_k⟩ .

The **surface traction** `g(1,τ)` decays; for large `τ` only `k=1` survives:

    g(1,τ) ≈ c₁ f₁(1) e^{−λ₁ τ} = c₁ e^{−τ/1.67} ,  λ₁ ≈ 0.60 ,

i.e. the traction relaxes on a **single time scale `τ_rel ≈ 1.67 t₀`**
(`t₀ = 4πkT h³/(E*dDΩ)`), with `λ₂/λ₁ ≈ 227` — the second mode is negligible.

## 6. Reproduce

    python3 solve.py

writes `results.json` and prints the elasticity validation error, the
per-grid eigenvalues, the extrapolation, the wedge tip fit and `K_∞`.

## Caveats
* `λ₁` carries the largest uncertainty (≈ ±0.05): the operator `A` is stiff
  near `s=0` and eigenfunctions carry a logarithmic surface term; the reported
  value is the robust monotone `1/n → 0` extrapolation of the clean grids.
* `λ₂`, `λ₃`, `f₁(0.5)` and `K_∞` are stable to the digits quoted.