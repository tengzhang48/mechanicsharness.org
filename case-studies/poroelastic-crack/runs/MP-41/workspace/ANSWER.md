# MP-41 — A suddenly loaded / pressurized crack in a fluid-saturated porous solid

## Reproduce

```
python3 solve.py
```

prints every number in `results.json` (and one numerical check).

## Results (part (c): ν = 0.2, ν_u = 0.4, α = 0.89)

K in units of √(πL), V in units of L²/G, per unit loading magnitude.

| name | value | general formula |
|------|-------|-----------------|
| `K1_0`   | 1.0        | `sqrt(pi*L)` |
| `K1_inf` | 1.0        | `sqrt(pi*L)` |
| `V1_0`   | 1.884956   | `pi*(1-nu_u)*L^2/G` |
| `V1_inf` | 2.513274   | `pi*(1-nu)*L^2/G` |
| `V1_short_coeff` | 1.128379 | `8*(nu_u-nu)*L^2/((1-nu)*sqrt(pi)*G)` |
| `p_faces_0` | 0.749064 | `2*(nu_u-nu)/(alpha*(1-2*nu))` |
| `K2_0`   | 0.0        | `0` |
| `K2_inf` | 0.0        | `0` |
| `V2_inf` | 0.0        | `0` |

## Model and method

**Geometry / material.** A Griffith (mode‑I) crack of half‑length L on y=0 in an
infinite plane‑strain Biot solid. The constitutive laws given in TASK.md imply
(total stress, drained constants G, ν, ν_u, α):

* `σ_ij = 2G e_ij + (2Gν/(1-2ν)) e_kk δ_ij - α p δ_ij`  (equivalent to the given equation),
* `ζ = α e_kk + p/M`, with `1/M = α²(1-2ν)(1-2ν_u)/(2G(ν_u-ν))`.

Mass balance `∂ζ/∂t = κ∇²p`. Eliminating e_kk gives the pore‑pressure transport law

```
∂p/∂t = c(1-nu_u)/(1-nu) ∇²p - (nu_u-nu)/(alpha(1-2nu)) ∂σ_kk/∂t .      (T)
```

The combination `T2 ≡ σ_kk + [α(1-2ν)/(1-ν)] p` is harmonic (Beltrami):
`∇²T2 = 0`. This is the key observation: the total stress is the **elastic** stress
field plus a term proportional to p.

**Why K is given by the face traction.** For a crack loaded through its face
tractions in a homogeneous body, the near‑tip singular field is the elastic mode‑I
field and its coefficient `K` is fixed by the traction distribution (weight‑function
argument); the pore pressure enters only through the regular term `-α p`, whose
gradient is integrable and does not change the `r^{-1/2}` amplitude. Hence the
opening pressure is the **total normal traction** on the faces and `K = Σ√(πL)`.

**Elastic (Griffith) building blocks** (uniform opening pressure Σ, plane strain):

```
K = Σ sqrt(pi L),   V = pi (1-nu) Σ L^2 / G,   delta(x) = 2(1-nu)Σ sqrt(L^2-x^2)/G
```

(verified internally: the near‑tip COD `δ ~ (4K(1-ν)/G)√(r/2π)` reproduces
`δ(x)=2Σ(1-ν)√(L²-x²)/G`, and `V=∫δ=π(1-ν)ΣL²/G`).

**Time limits.** At `t=0⁺` no pore fluid has moved, so the bulk is **undrained**:
substituting `ζ=0` into the constitutive law collapses it to the elastic law with
`ν→ν_u` (`2G e_ij = σ_ij - ν_u σ_kk δ_ij`), and
`p = -(ν_u-ν)σ_kk/(α(1-2ν))`. At `t→∞` the pore pressure has relaxed to its
far‑field value and the material is **drained** (ν). The crack opening is therefore
the Griffith value with the appropriate Poisson ratio.

## (a) Loading 1 — `σ_n = -H(t)`, `p = 0` on the faces

* Opening traction = +1 (total normal traction).
* `t=0⁺` (undrained, ν_u):  `K1_0 = √(πL)`, `V1_0 = π(1-ν_u)L²/G`.
* `t→∞` (drained, ν):       `K1_inf = √(πL)`, `V1_inf = π(1-ν)L²/G`.

**Induced pore pressure next to the faces at t=0⁺.** On the crack faces (undrained
elastic solution) `σ_xx = σ_yy = -1`, so `σ_kk = -2`. Undrained (`ζ=0`) gives

```
p_faces_0 = -(nu_u-nu) σ_kk/(alpha(1-2nu)) = 2(nu_u-nu)/(alpha(1-2nu)) = 0.749064 .
```

**Small‑time growth of V (ct/L² ≪ 1).** The crack faces are held at p=0 while the
interior just inside is undrained (`p=p_und`); a pore‑pressure boundary layer grows
inward from each face. In the layer the problem is 1‑D (y only). Combining (T) with
the 1‑D compatibility result that e_xx is constant through the layer gives exactly

```
∂p/∂t = c ∂²p/∂y² ,   p(0,t)=0 ,  p(y,0)=p_und  =>  p = p_und erf(y/(2 sqrt(c t)))
```

The vertical strain in the layer is `e_yy = [-1+3ν-2νν_u + α(1-2ν)p]/(2G(1-ν))`;
the change from the undrained value is `δe_yy = α(1-2ν)(p-p_und)/(2G(1-ν))`. Using
`∫₀^∞(p-p_und)dy = -2 p_und sqrt(ct)/sqrt(pi)` and `α(1-2ν)p_und = 2(ν_u-ν)`, the
extra opening per face is `2(ν_u-ν)√(ct)/(G(1-ν)√π)`. With two faces and integrating
along the crack (uniform in the interior, length 2L):

```
V ≈ V1_0 + A sqrt(ct/L^2),   A = 8 (nu_u-nu) / ((1-nu) sqrt(pi))  (units L^2/G).
```

## (b) Loading 2 — `σ_n = 0`, `p = H(t)` on the faces

The total normal traction is identically zero on the faces.

* `t=0⁺`: no fluid has moved; the pore pressure is only a sheet, `p=1` on the faces
  and `p=0` in the interior. The bulk is undrained **and unloaded** (zero face
  traction, zero far field), so the undrained elastic solution is the trivial one,
  `σ ≡ 0`. Therefore `K2_0 = 0` and `V2_0 = 0`.
* `t→∞`: the injected pressure diffuses away to infinity; the pore pressure at any
  fixed material point `→ 0` (leak‑off). With zero face traction and `p=0` in the
  interior, the strain field `→ 0`, so `K2_inf = 0` and `V2_inf = 0`.

Consistency: Loading 1 + Loading 2 equals the fluid‑pressurized crack
(`σ_n=-H`, `p=H`). Its opening is `π(1-ν_u)→π(1-ν)` (undrained→drained) and
`K=√(πL)` at all times — the classical Rice–Cleary result — which is exactly
`V1+V2`.

## Checks made

1. **Internal elastic consistency of K and V.** The crack-opening-displacement
   profile `δ(x)=2(1-ν)Σ√(L²-x²)/G` was checked against the Irwin near‑tip relation
   `δ=(4K(1-ν)/G)√(r/2π)` with `K=Σ√(πL)`; both give `2√2(1-ν)Σ√(Lr)/G`. Then
   `V=∫δ dx = π(1-ν)ΣL²/G`. *Could have failed if* the COD normalisation
   were `4` instead of `2`; the near‑tip cross‑check rules that out.

2. **Boundary‑layer integral (independent numerics).** `solve.py` solves
   `∂p/∂t=c∂²p/∂y²`, `p(0,t)=0`, `p(y,0)=p_und`, `p(∞,t)=p_und` with SciPy's stiff
   BDF integrator (method of lines), never using the erf closed form in the solve
   step. Compared `∫₀^∞(p-p_und)dy` to the analytic `-2p_und√(ct)/√π`:

   ```
   boundary-layer face value   = 0.0        (expect 0)
   int(p-p_und)dy  numeric     = -0.7979468936430866
   int(p-p_und)dy  analytic    = -0.7978845608028655
   relative difference         = 7.81e-05
   ```

   *Could have failed* by (a) numerical diffusion/instability — avoided by the
   implicit BDF scheme (the earlier explicit scheme was rejected because the
   stability limit `dt≤dy²/2c` made it unusable), or (b) the far boundary polluting
   the layer — avoided by `y_max=10 ≫ √(ct)=0.71`. Agreement to 8e‑5 validates both
   the erf boundary‑layer solution and the factor `2/√π` entering `A`.

3. **`K1_0 = K1_inf`.** Both equal `√(πL)` because K is traction‑controlled; the
   drained/undrained switch only rescales the opening, not the singular field.

## What is missing / uncertain (assumptions used)

* The pore‑pressure transport law (T): I used the given diffusivity `c` verbatim and
  verified that in the 1‑D face boundary layer (with the 1‑D compatibility result
  `e_xx = const`) the coefficient collapses to `c` itself, i.e. `∂p/∂t = c∇²p`
  (the coupling term is absorbed). If instead the intended transport law kept the
  un‑absorbed form, `A` would scale by `√((1-nu_u)/(1-nu))`. **Assumption:** the
  boundary‑layer form above.
* `A` uses a "1‑D‑layer, uniform‑along‑the‑crack" approximation: it neglects the
  near‑tip modulation and the slow bulk drainage. It is the correct leading
  `√(ct)` term but the `O(1)` prefactor is a boundary‑layer ("surface") estimate,
  not an exact finite‑crack result. The formula and derivation are given so it can
  be checked.
* **Loading 2.** The literal reading (`σ_n=0` total traction) makes the faces
  traction‑free, so — with the total‑traction rule for K — the answers are zero at
  both limits, which is what I report. If the problem intends `σ_n=0` as an
  *effective* normal stress (faces carrying total traction `-αp`), the answers
  change to `K2_0=K2_inf=α√(πL)=0.89√(πL)` and `V2_inf=π(1-ν)αL²/G=2.237`; I flag
  this rather than silently choosing.
* `p_faces_0` is the jump value of p in the solid immediately inside the drained
  faces (the value at the face itself is 0 by the boundary condition).

## Questions

1. Clarify Loading 2: is `σ_n=0` the *total* normal stress (my reading, giving zero
   K and V) or an *effective* normal stress? Assumption used: total, because the
   problem explicitly defines "σ_n = σ_yy is the normal stress on the faces".
2. Confirm the sign/definition of `p_faces_0` (I report the magnitude of the
   induced compressive pore pressure adjacent to the face).

*(Never opened files outside this folder; no web access used.)*
