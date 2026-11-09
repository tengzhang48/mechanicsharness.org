# ENTRY.md — Knowledge unit: HRR near-tip field at a tensile crack in a hardening material

Source: J. W. Hutchinson, *Singular behaviour at the end of a tensile crack in a
hardening material*, J. Mech. Phys. Solids **16** (1968) 13–31 (P015). This is the
plane-strain/plane-stress companion of Rice–Rosengren (1968); the near-tip field is
now known as the **HRR field**.

## Problem

A **central finite tensile crack of half-length `L`** in an infinite plate under a
far uniaxial tension normal to the crack (traction-free faces); the near-tip region
is treated as locally semi-infinite. Lengths are non-dimensionalised by the crack
half-length `L`, stresses by the yield stress `σ_y`, and strains by the yield strain
`ε_y = σ_y/E` (paper §2, Fig. 2). The material is elastic–plastic, described by
**total deformation theory** with power-law (Ramberg–Osgood) hardening
`ε = σ + α σⁿ`. Find the dominant stress/strain singularity as the crack tip is
approached, where plastic strains dominate, and determine its amplitude under
**small-scale yielding** by means of Rice's path-independent `J`-integral. Solve
for generalized **plane stress** and **plane strain**.

## Assumptions and range of validity

- Total (deformation) plasticity, valid only for **monotonic, nearly proportional
  loading and no unloading**. The paper examines this only a posteriori and does
  **not** establish it: §8 states that the stress proportions in the elastic and
  fully-plastic zones are "quite different" on the same ray and calls for a check
  against Budiansky's acceptability criterion (pp. 14, 29). Nearly-proportional
  loading therefore remains an assumption to be validated for each use case.
- The loss of J's energetic interpretation applies to **unloading during crack
  extension** in a real metal (p. 20); it is a separate issue from the validity of
  the local singularity, which (26) extends even when the plastic zone is not small
  provided the outer integral does not vanish (p. 21).
- Small elastic strains in the singular zone; high-`n` hardening, so the plastic
  `n`-power term controls the tip.
- Two dimensions, plane stress or plane strain; linear kinematics and equilibrium.
- **Small-scale yielding is needed only for the far-field amplitude `K`**: the
  plastic zone must be small compared with `L` and embedded in the elastic singular
  field, so that J from the outer path is the elastic `πσ^{∞2}` (plane stress) or
  `π(1-ν²)σ^{∞2}` (plane strain). Rice (1966) suggested validity for applied stress
  below ~1/3, at most ~1/2, of yield. Loads near fracture violate this amplitude
  relation; they do not by themselves invalidate the local `r`-powers.
- Symmetry about the crack line (mode I); traction-free crack faces.

## Main results

### 1. Singularity exponent  — `derived` and `reproduced`

`φ = K rˢ f(θ)` with the exact, `n`-independent-of-details result

    s = (2n+1)/(n+1)      (eqs. 15, 26)

Consequences: stresses `σ ~ r^{-1/(n+1)}`, strains `ε ~ r^{-n/(n+1)}`, and the
**strain-energy density `W ~ 1/r`** for every `n`.

- `reproduced`: after substituting `φ = K rˢ f(θ)` in the governing equations, the
  fourth-order angular boundary-value problem (both plane stress eq. 13 and plane
  strain eq. 31) with `f'(0)=f'''(0)=0`, `f(π)=f'(π)=0` was solved by shooting on
  `(s, f''(0))`. The converged `s` equals `(2n+1)/(n+1)` to `<10^{-6}` for
  `n=3,5,9,13` in both cases (test `test_exponent_from_eigenvalue`).
- `derived`: at `s=(2n+1)/(n+1)`, `(n+1)(s-2) = -1`, so the **energy density scales
  as `W ~ 1/r`**, not as a radius-independent quantity. The J integral has one extra
  factor `r` from the contour element `ds`, so the radial power of its integrand is
  `(n+1)(s-2)+1 = 0` and J is radius-independent. (The earlier statement that `W` is
  `r`-independent and that the equation-(22) exponent vanishes at the eigenvalue was
  wrong; only `J`, not `W`, is radius-independent.)

### 2. Angular stress fields  — `reproduced` (qualitative) and `derived`

The normalised angular functions `σ̃_r, σ̃_θ, σ̃_rθ, σ̃_e(θ)` were obtained from the
same ODE. They match the shapes of the paper's Figs. 3a/b (plane stress) and 5a/b
(plane strain): e.g. plane stress `n=3` gives `σ̃_r(0)=0.77`, `σ̃_θ(0)=1.11`,
`σ̃_e(0)=0.98` (figure reads ≈0.75, ≈1.05, ≈1.0); plane strain `n=3` gives
`σ̃_θ(0)=1.94` (computed 1.9407), `σ̃_r(0)=1.65`, `σ̃_e` peaking near `π/2`
(Fig. 5a reads ≈1.94 and ≈1.65 at `θ=0`, peak near `π/2`). Figure comparison is
qualitative (read by eye); the earlier "≈1.4" for `σ̃_r(0)` was a misreading of
Fig. 5a and "1.95" rounds `σ̃_θ(0)` too high.

### 3. J-integral normalisation `I` and energy coefficient `E_n`  — `reproduced`

Using the paper's eq. (24) with the printed similarity displacements and
`σ̃_e^{max}=1`:

| n | case | `I` (ours) | `I` (paper) | rel. err | `E_n` (ours) | `E_n` (paper) | rel. err |
|---|------|-----------|-------------|----------|--------------|---------------|----------|
| 3 | plane stress | 3.855 | 3.86 | −0.13% | 2.162 | 2.16 | +0.11% |
| 5 | plane stress | 3.408 | 3.41 | −0.07% | 2.447 | 2.45 | −0.13% |
| 9 | plane stress | 3.030 | 3.03 | −0.01% | 2.671 | 2.67 | +0.05% |
| 13 | plane stress | 2.871 | 2.87 | +0.04% | 2.755 | 2.76 | −0.17% |
| 3 | plane strain | 5.507 | 5.51 | −0.05% | 1.227 | 1.23 | −0.23% |
| 5 | plane strain | 5.024 | 5.01 | +0.27% | 1.382 | 1.37 | +0.85% |
| 9 | plane strain | 4.602 | 4.60 | +0.04% | 1.485 | 1.48 | +0.34% |
| 13 | plane strain | 4.401 | 4.40 | +0.02% | 1.518 | 1.52 | −0.10% |

Every `I` is reproduced to ≤0.3% and every `E_n` to ≤0.85%. The associated
column `(π/I)^{1/(n+1)}` (Table 1) is reproduced to `<3×10^{-3}`. The test
tolerances are 0.5% for `I` and 1.0% for `E_n`, matching these measured
agreements; they are bounds on **agreement with the printed table**, not error
estimates for the present computation. Refining the quadrature/ODE changes
`E_n` by only ~1.1e-8 (e.g. plane strain `n=5`: `1.3816415028` at 4001 points vs
`1.3816415137` at 16001 points), so the 0.85% gap from the printed `1.37` is
**not** a current integration error; the cause of the historical table difference
is unestablished.

### 4. Small-scale-yielding amplitude  — `taken` (formula) + `reproduced` (column)

    K = (1/α)^{1/(n+1)} (π/I)^{1/(n+1)} (σ^∞)^{2/(n+1)}          (plane stress, eq. 27)
    K = ((1-ν²)/α)^{1/(n+1)} (π/I)^{1/(n+1)} (σ^∞)^{2/(n+1)}    (plane strain, eq. 32)

`taken` from J-matching; internally consistent with the reproduced `I`.

### 5. Plane strain has the higher tensile stress and the lower fracture toughness — `taken`

The ratio `σ_θ(0)_{p.strain}/σ_θ(0)_{p.stress}` (eq. 50, Fig. 7) exceeds one and
grows with `n`; the paper reports the tensile singularity ahead of the crack is
larger in plane strain (paper §7, p. 27). The physical inference drawn there
(pp. 28) is the opposite of "tougher": the plane-strain **fracture toughness**, as
reflected by the tensile singularity ahead of the crack, is significantly **lower**
than the plane-stress toughness. This is a near-tip field comparison, not itself a
fracture criterion. (Not reproduced numerically here; the authors flag the
single-term piecewise-linear version of this as "suspicious".)

## Conditions and limits

- Valid only as `r→0` inside the plastic zone; the dominant term is only an
  accurate representation in a very small region.
- `n` finite and `>1`; `n=1` is linear elastic and the plastic/elastic split used
  here degenerates (the elastic biharmonic term reappears).
- `E_tan = 0` (perfect plasticity) is a singular limit **not** described by this
  solution; the stress becomes bounded and perfect-plasticity slip-line theory
  applies (Rice 1967).
- Applicable to smooth-yielding materials (e.g. aluminium), not to the
  Lüders-band mode of mild steel (paper §7).

## Relations to other results

- Specialises to the **linear-elastic** crack field at `n=1`/as the plastic zone
  shrinks: `φ ∝ r^{3/2}(cos θ/2 + ⅓ cos 3θ/2)` (eq. 40), with `s=3/2`. The exported
  `elastic_J_plane_stress` is the closed form `K_I²/E`; the exported
  `elastic_J_numeric` evaluates the circular `J` integral of that field from first
  principles (full Cartesian derivatives, including the polar-basis terms) and the
  two agree to quadrature accuracy (`test_elastic_J_plane_stress_export_matches_numeric`).
- The exponent is the same for plane stress and plane strain and equals the
  Rice–Rosengren (1968) result; the angular functions and `I_n` differ between the
  two.
- `I` and `K` are distinct objects: `I` is a dimensionless angular normalisation
  integral depending only on `n` and the constraint (plane stress/strain); the HRR
  **amplitude `K`** is the load intensity, fixed by J (and `I`) through
  eqs. (25)/(27)/(32). `K` is the HRR analogue of the elastic `K_I`; `I` is its
  normalisation. The plane-strain `I_n`, `E_n` are the standard normalisations used
  in later HRR/J-integral estimation schemes.
- The biharmonic control of an elastic singularity and the `J`-integral are due to
  Rice (1967/1968); McClintock & Irwin (1965) give the background.

## Errata / inconsistencies found

1. **`paper.txt` OCR error**: plane-stress `I(n=13)` is printed `2.87` in the
   image (and is required by the tabulated `(π/I)^{1/14}=1.006`); `paper.txt`
   gives `2.37`.
2. **Digit ambiguity**: plane-strain `E_n(n=9)` is read as `1.48` in the image
   (our computed `1.485`); `paper.txt` OCR gives `1.43`. Similarly the OCR of the
   plane-strain `E_n` row is unreliable; our computed `1.227,1.382,1.485,1.518`
   match `1.23,1.37,1.48,1.52`.
3. **Eq. (22), traction-work line — printed prefactor error.** The image prints
   `σ_{ij}n_ju_{i,x} = K^{n+1} r^{(n-1)(s-2)}{...}`, i.e. it **omits `α`** and uses
   the exponent `(n-1)(s-2)`. The W line of the same equation correctly carries
   `αK^{n+1} r^{(n+1)(s-2)}`. Dimensionally (displacements carry `αKⁿ`) and by
   comparison with (23)/(24) the traction-work line must be
   `αK^{n+1} r^{(n+1)(s-2)}{...}`. This is a genuine printing/derivation slip; eqs.
   (23)–(27) and the implemented `I` use the correct scaling. See `EXTRACT.md`.
4. **Eq. (22) and a first-principles J integral — not a discrepancy.** The earlier
   "open question" (a direct J integral giving `I≈5.65` for plane stress `n=3`) was
   an artifact of omitting the derivatives `ũ_r'`, `ũ_θ'` of the polar basis /
   displacement angles. Re-doing the Cartesian derivatives carefully, or a
   high-resolution independent quadrature, gives `I=3.8549663167`, reproducing Table
   1. The printed eq. (24) bracket (with the `-\sinθ(...)` term) is correct and is
   what `integral_I` implements. The direct route is validated by
   `test_elastic_mode_I_J_equals_KI2_over_E_numeric` and
   `test_elastic_J_plane_stress_export_matches_numeric` on the elastic special case.
5. **Eq. (37) printed without a square.** The image prints the hydrostatic term as
   `((1-2ν)/6) σ_KK`; the complementary-energy form requires `σ_KK²` (from
   `σ_{pp}²`). Kept as printed in `EXTRACT.md`, noted here.
6. **Eq. (44) printed differential.** The image prints `∫ r dr ∫ W dσ`, i.e. the
   angular differential is printed `dσ` rather than `dθ`. Kept as printed in
   `EXTRACT.md`, noted here.

## Using this unit

**Questions it answers.** (i) What is the dominant stress/strain singularity at a
mode-I crack tip in a power-law hardening material? (ii) What are the angular
stress distributions and the dimensionless normalisations `I`, `E_n`? (iii) What is
the small-scale-yielding amplitude `K` in terms of J (or `σ^∞`)? (iv) How does the
plane-strain tensile stress ahead of the tip compare with plane stress?

**Input conventions and units.** The HRR quantities are dimensionless: stresses
in units of the yield stress `σ_y`, strains in units of `ε_y = σ_y/E_ref`,
lengths in units of the crack half-length `L` (`r = r_phys/L`). The remote load
entering `amplitude_K` is `sigma_inf = σ_remote/σ_y` (positive). The
Ramberg–Osgood law is `ε = σ + α σⁿ` with `α > 0` the dimensionless hardening
coefficient and `n > 1` the hardening exponent. `n=1` is linear elastic (the
plastic split used here degenerates, so it is outside the intended range);
`α > 0` and the normalisation `I > 0` are required. `nu` is Poisson's ratio
(default `0.3`). `plane_strain` is a boolean selecting the constraint. `theta`
is measured from the crack line, in radians, `0 <= theta <= pi` by symmetry.
The physical stress field is
`sigma_phys(r,theta) = sigma_y * K * (r_phys/L)^(s-2) * sigma_tilde(theta)`,
where `sigma_tilde` is the normalised angular function returned by
`solve_angular` and `K` is the dimensionless HRR amplitude. Note `σ_y` is the
unit of stress; the reference modulus used to make `ε_y` is independent, so no
extra factor of `E` appears in the HRR formulas.

**Elastic exceptions to the unit convention.** The two elastic functions
`elastic_J_plane_stress` and `elastic_J_numeric` take the *conventional*
linear-elastic stress-intensity factor `K_I` (units stress·sqrt(length)) and
Young's modulus `E` (units stress), and return `J` in stress·length; they are
not dimensionless. If instead one normalises as in the paper,
`K_I_dimless = K_I_phys/(σ_y sqrt(L))`, `E` is measured in units of the reference
Young modulus (`E=1` when normalised), and
`J_phys = σ_y ε_y L * J_dimless`. In the HRR formulas `K` (the dimensionless
amplitude fixed by `J` through eqs. 25/27/32) must be distinguished from this
conventional elastic `K_I`: `K_I` is an input to the elastic helpers, `K` is an
output of `amplitude_K`.

The remaining elastic helpers are **plane-stress, unit-amplitude** helpers:
`elastic_angular_stress_function`, `elastic_angular_displacements`,
`elastic_J_integrand` and `elastic_displacement_cartesian` use the mode-I
angular field of eq. (40) normalised so that its unit stress amplitude
corresponds to the conventional `K_I = sqrt(2*pi)` (so that the dimensionless
circular integral equals `2*pi`). `theta` is in radians; `nu` is Poisson's
ratio; the modulus `E` must be positive (the displacement helper scales as
`1/E`); `r > 0` is required for `elastic_displacement_cartesian` (it uses
`r = hypot(x,y)` and `theta = atan2(y,x)`).

**Functions to call** (module `hrr`, all in `__all__`):

```
singular_exponent(n) -> s = (2n+1)/(n+1)
effective_stress_coeff(n, plane_strain) -> 3x3 M with se_tilde^2 = z^T M z, z=(f,f',f'')
solve_angular(n, plane_strain=True, s=None, normalize=True, n_theta=4001) -> dict
    {theta,f,fp,fpp,fppp,se,sr,st,srt,c,p,s}
integral_I(n, plane_strain=True, sol=None) -> I               (eq. 24)
energy_coefficient_En(n, plane_strain=True, sol=None, I=None) -> E_n   (eq. 34)
amplitude_K(n, plane_strain=True, alpha=1.0, sigma_inf=1.0, nu=0.3, I=None) -> K
elastic_angular_stress_function(theta) -> (f,f',f'',f''')     (eq. 40)
elastic_angular_displacements(theta, nu=0.3) -> (ur,ut,ur',ut')  (eq. 40)
elastic_J_integrand(theta, nu=0.3) -> (integrand, terms)      (E=1 integrand)
elastic_displacement_cartesian(x, y, nu=0.3, E=1.0) -> (ux, uy)  (unit amplitude)
elastic_J_plane_stress(K_I, E) -> K_I^2/E                     (closed form)
elastic_J_numeric(K_I, E=1.0, nu=0.3, n_theta=200001) -> J    (first-principles integral)
```

**Meaning of the `solve_angular` keys.** `theta` is the angular grid;
`f,fp,fpp,fppp` are the angular function and its first three derivatives;
`se,sr,st,srt` are the dimensionless effective, radial, hoop and shear stresses
`sigma_tilde`; `s = (2n+1)/(n+1)` is the eigenvalue; `p = f''(0)` is the
shooting value. `c` is the normalisation factor applied to the raw solution:
with `normalize=True` (default) the returned field has `max_theta se = 1`
(`c = 1/max(se_raw)`); with `normalize=False` `c = 1` and the raw field is
returned, whose n=3 plane-stress `max(se)` is `1.1866376456455583`. The raw
fields equal the normalised fields divided by `c`; `p` is unaffected by the
choice.

**Optional arguments, exponent override and angular grid.** Only `integral_I`
and `energy_coefficient_En` accept a precomputed `sol`; `amplitude_K` does
**not** accept `sol` (calling `amplitude_K(..., sol=...)` raises `TypeError`).
Only `energy_coefficient_En` and `amplitude_K` accept a supplied `I`; `integral_I`
computes `I` and accepts a precomputed `sol` (it does not accept `I`), so an
independently obtained normalisation (for example the paper's Table 1 value)
can be supplied; `amplitude_K` also accepts `sigma_inf`, `alpha` and `nu`, but
**not** a `J` argument. When passing `sol`/`I`, `n`, `plane_strain` and the
`normalize` setting must be consistent with the object passed (the functions do
not re-check this). `s` is an optional fixed-exponent override: leaving it
`None` uses the physical HRR eigenvalue `s(n)=(2n+1)/(n+1)`. `n_theta` is the
number of uniform angular samples; it sets the output grid and the quadrature
accuracy, not the ODE tolerance.

**Normalisation scope.** A tabulated value of `I` is a property of `n` and the
constraint **only after the angular normalisation is fixed**. The standard
values (and Table 1) are for a field normalised to `max_theta se = 1`, i.e. the
`normalize=True` default. `I` is homogeneous of degree `n+1` in the angular
field, so it scales as the field scale to the power `n+1`. For the same `n=3`
plane-stress field the default normalised value is `I = 3.854966353`, whereas
the raw (`normalize=False`) field with `c = 1` gives `I = 7.643514258`; in
general `I_raw = I_norm / c_norm**(n+1)`, where `c_norm = 1/max(se_raw)` is the
`c` returned by the **normalized** solution. `I`, `K` and the angular field must
therefore share one common normalisation before they are combined.

**Amplitude reconstruction.** With the inner-contour identity
`J_dimless = alpha K^(n+1) I`, the amplitude follows from a known dimensionless
J as
`K = (J_dimless/(alpha*I))^(1/(n+1))` for **both** plane stress and plane
strain. The factor `1-nu^2` belongs only in computing the *outer* elastic J from
the remote load, `J_outer = pi sigma_inf^2` (plane stress) or
`J_outer = pi (1-nu^2) sigma_inf^2` (plane strain); it must **not** be applied
again when inverting J for K. `amplitude_K` implements the equivalent
remote-load form of eqs. (27)/(32). Worked check: `n=3`, `alpha=2`,
`sigma_inf=0.2`, `nu=0.3`, plane-strain `I=5.51` gives
`J_outer = pi(1-nu^2)sigma_inf^2 = 0.1143539726` and `K = 0.319166584`; the
round-2 documented extra division by `1-nu^2` gives `0.326781`, 2.39% high.
This is bound by `test_amplitude_K_documented_reconstruction`.

**Demonstrated range.** The numerical results here were demonstrated for
`n = 3, 5, 9, 13`, both constraints. The asymptotic solution is defined
theoretically for any finite `n > 1`; values outside the demonstrated set are
extrapolation.

**Where the model holds, and where it does not.**
- Holds: the asymptotic `r -> 0` field inside the plastic zone, for monotonic,
  nearly-proportional, non-unloading loading, `n` finite and `>1`, small elastic
  strains. The `n`-power plastic term controls the tip.
- Does not hold: `E_tan = 0` (perfect plasticity) is a singular limit — there the
  stress is bounded (slip-line theory, Rice 1967). `n=1` is elastic. The
  small-scale-yielding amplitude `K` fails at loads near fracture, where the plastic
  zone is not small (this is a limit on the amplitude relation, not on the local
  exponents). The energetic meaning of J fails under unloading during crack
  extension in a real metal. Nearly-proportional loading is **assumed**, not proven
  (paper §8). Mild-steel Lüders-band deformation is not modelled.

**Provenance of each main result.**
- Exponent `s=(2n+1)/(n+1)`: **re-derived** from the angular eigenvalue problem
  (shooting), and **independently checked** by the review.
- Angular stress fields: **reproduced** (solved from the ODE); figure comparison is
  qualitative only.
- `I` (eq. 24): **implemented** from the printed integral; **reproduced** against
  Table 1 to ≤0.3%; **independently checked** by a Cartesian finite-difference J
  route (`test_integral_I_matches_cartesian_J`, agreement < 1.3e-7 relative) and by
  the review's independent derivation.
- `E_n` (eq. 34): **implemented**; **reproduced** against Table 1 to ≤0.85%.
- Amplitude `K` (eqs. 27/32): **taken** from J-matching, consistent with the
  reproduced `I`; **checked** in `test_amplitude_K_satisfies_J_balance` against the
  paper's `I`.
- Elastic special case (eq. 40): **implemented**; the `J=K_I²/E` identity is
  **independently checked** by a first-principles circular J integral
  (`elastic_J_numeric`), whose **returned** integrand is compared pointwise to
  `W cos(theta)` minus the traction work evaluated from full Cartesian
  displacement gradients (`test_elastic_reduced_integrand_matches_cartesian_components`,
  residual `< 1e-8`; it fails, residual `1.4`, on the old derivative expression)
  and shown to integrate to `2π`
  (`test_elastic_dimensionless_integral_is_2pi`). The displacement helper's
  `1/E` scaling is checked at `E != 1`
  (`test_elastic_displacement_cartesian_modulus_scaling`).
- Eq. (50) / Fig. 7 plane-strain vs plane-stress ratio: **taken**, not reproduced.
