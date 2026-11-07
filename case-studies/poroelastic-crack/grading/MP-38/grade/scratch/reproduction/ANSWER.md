# Suddenly Pressurized Poroelastic Crack

## Results

Use dimensionless time `tau = ct/L^2`, intensity `k = K/sqrt(pi L)`, and volume `v = GV/L^2`. The crack is traction-prescribed, without a contact constraint. Define

`eta = alpha(1-2 nu)/(2(1-nu))`.

| Quantity | General expression | nu=0.2, nu_u=0.4, alpha=0.89 |
|---|---|---:|
| K1_0 | Actual-tip limit of the spectral system below | approximately 0.808 |
| K1_inf | 1 | 1 |
| V1_0 | pi(1-nu_u) | 1.8849555921538759 |
| V1_inf | pi(1-nu) | 2.5132741228718345 |
| V1_short_coeff | High-frequency volume limit below | approximately 2.44 |
| p_faces_0 | 2(nu_u-nu)/(alpha(1-2nu)) | 0.7490636704119851 |
| K2_0 | 0 | 0 |
| K2_inf | -eta | -0.33375 |
| V2_inf | -pi alpha(1-2nu)/2 | -0.8388052385084748 |

Thus, for loading 1,

`V(t) = (L^2/G)[pi(1-nu_u) + A sqrt(ct/L^2) + o(sqrt(ct/L^2))]`,

with `A approximately 2.44` for the specified parameters. Numerical uncertainty estimates are about `0.001` in `K1_0` and `0.01` in `A`; these are convergence estimates, not rigorous bounds.

**Important initial-tip convention:** the intensity of the initial *outer undrained elastic field* is exactly 1. The intensity at the actual tip, obtained by taking distance to the tip to zero at each positive time before taking time to zero, is approximately 0.808. These are different limits. The table and JSON use the latter, because `K` is defined at the crack tips. If `t=0+` instead means the strictly undrained initial field with diffusion suppressed, use `K1_0=1` and leave the other reported values unchanged.

## Model and Brief Derivation

The given two-index constitutive equations reduce to

`sigma_ij = 2G e_ij + lambda e_kk delta_ij - alpha p delta_ij`,

`zeta = alpha e_kk + S p`,

where `lambda = 2G nu/(1-2nu)` and `S = alpha^2(1-2nu)(1-2nu_u)/(2G(nu_u-nu))`. The undrained Lamé constant is `lambda_u = lambda + alpha^2/S = 2G nu_u/(1-2nu_u)`, and `kappa/(S+alpha^2/(lambda+2G))` reproduces the stated diffusivity. Thus the two-index equations are used as written, not as an unreduced three-dimensional compliance law.

For loading 1 the initial outer solution has zero fluid-content change and is undrained. The late solution is drained. The elastic crack solution under face traction `-P` is

`K=P sqrt(pi L)`, `w(x)=2(1-nu)P sqrt(L^2-x^2)/G`, `V=pi(1-nu)P L^2/G`.

Replacing `nu` by `nu_u` gives the initial volume. The initial outer face stress trace is `-2`; setting `zeta=0` gives the listed induced pressure. This replacement gives the *outer* initial intensity, not necessarily the actual intensity inside the shrinking drainage region.

For loading 2, pressure spreads out to a distance of order `sqrt(ct)`. At any fixed point its late limit is 1, but at every finite time it still decays at infinity. The asymptotically radial distant field satisfies `(lambda+2G)e_kk=alpha p`. Its local stress tends to `-eta delta_ij`. A drained crack correction with face stress `+eta` cancels that background traction. Its intensity and opening are negative, giving the listed loading-2 values. Initially the pressure occupies a vanishing layer and `K2_0=0`.

## Questions and Assumptions

- The crack is assumed frictionless, with zero shear traction, and symmetry across its plane. No unilateral contact constraint was specified. Loading 2 can therefore have negative mathematical opening.
- `p_faces_0` means the solid-side, undrained outer pressure before diffusion, not the actual Dirichlet boundary value. The latter is exactly zero for loading 1 at every positive time; these spatial and time limits do not commute.
- A decaying steady harmonic pressure equal to 1 on a bounded crack does not exist in the infinite two-dimensional exterior. The requested long-time results are limits of the specified transient problem, not of an inconsistent steady boundary-value problem.
- The incompressible endpoint `nu_u=1/2` is interpreted by continuity.
- The initial-tip limit is not specified explicitly. I use the actual-tip limit `lim(t->0+) lim(r->0+) sqrt(2 pi r) sigma_yy(L+r,0,t)`. Reversing these limits gives the outer undrained value 1. Equating drained and undrained elastic energy-release rates is not justified without accounting for Darcy dissipation in the tip layer.

## Transform Method

This method is chosen to resolve the hydraulic crack-tip layer as well as the long, flat face layers. A one-dimensional face-only diffusion approximation misses contributions of the same order as the requested volume growth. The two difficult general expressions are supplied as explicit spectral limits, not closed elementary functions.

Set `L=G=1` and use Laplace transform in `tau=ct/L^2`. Write `q=sqrt(s)`, `m=sqrt(k^2+q^2)`, `a=1-nu_u`, and `d=2(nu_u-nu)/(1-nu)`. All responses below are step-normalized, i.e. multiplied by `s` after Laplace transformation. The coefficients `u_j(q)` represent the upper-face displacement, so the volume is `pi u_0(q)`.

Use even crack-supported bases

`phi_j(x)=sqrt(1-x^2) U_(2j)(x)` and `psi_j(x)=T_(2j)(x)/sqrt(1-x^2)`.

Their Fourier transforms are `(-1)^j Phi_j` and `(-1)^j Psi_j`, where

`Phi_j(k)=pi(2j+1) J_(2j+1)(k)/k`, `Psi_j(k)=pi J_(2j)(k)`.

The alternating signs are absorbed into the unknown coefficients in the code. Let `D_ij=pi(2i+1)delta_ij/2`. The matrices are

`A_ij = D_ij - (d/pi) integral_0^infinity [k^3/(m(m+k))] Phi_i Phi_j dk`,

`B_ij = (d/pi) integral_0^infinity [q^2 k/(m(m+k))] Phi_i Psi_j dk`,

`C_ij = (d/pi) integral_0^infinity [q^2/m] Psi_i Psi_j dk`.

Solve

`(A - B C^(-1) B^T) u = a (pi/2) e_0`.

For `d=0`, discard the hydraulic system and use `A=D`. This avoids a singular but physically irrelevant zero hydraulic matrix.

These symbols follow by solving equilibrium and mass balance for decaying Fourier modes in the upper half-plane. The elastic harmonic modes decay as `exp(-k y)` and the diffusion mode as `exp(-m y)`. Eliminating their amplitudes using zero shear traction gives the above coupled normal-displacement/flux boundary system. Galerkin restriction to the crack enforces zero displacement and zero normal flux on its intact continuation, and the specified traction and pore pressure on the crack. In particular, the hydraulic inverse must be taken **after restriction to the crack**; replacing it with a full-line Fourier inverse gives an incorrect tip response.

The explicit general formulas for the remaining requested quantities are

`K1_0/sqrt(pi L) = lim(q->infinity) [sum_(j>=0) (-1)^j(2j+1) u_j(q)]/(1-nu)`,

`A = (2/sqrt(pi)) lim(q->infinity) q[pi u_0(q)-pi(1-nu_u)]`.

The infinite-basis limit is taken before the `q` limit. The factor `2/sqrt(pi)` follows from `s Laplace[sqrt(tau)] = sqrt(pi)/(2q)`. The intensity sum follows from `U_(2j)(1)=2j+1` and the drained displacement-intensity relation at any positive time. These expressions depend on `nu,nu_u` only; `alpha` enters dimensional time through `c`.

`solve.py` uses composite Simpson quadrature and stable low-argument Bessel evaluation with forward recurrence only at high arguments. The omitted Fourier tails are estimated by their oscillatory averages. If the cutoff is `R`, the corrections to `A` and `C` are respectively

`-d outer((2j+1)(-1)^j)/[sqrt(R^2+q^2)+R]`,

`d q asinh(q/R) outer((-1)^j)`.

The leading mixed tail is oscillatory and is not added. Cutoff doubling checks its residual influence. Production uses 80 basis functions, cutoff 160000, spacing 0.125, and `q=80,320,1280,5120`. A cubic polynomial in `q^(-1/2)` extrapolates volume, and a quadratic through the first three points extrapolates intensity. The last intensity point is excluded because fixed-basis tip resolution begins to deteriorate there.

## Missing or Uncertain

- No elementary closed form or rigorous numerical error bound has been established for the initial actual-tip intensity or short-time coefficient. Their general answers are the fully specified spectral limits above; numerical values are rounded to the supported precision.
- Tip convergence is more sensitive than volume convergence. Higher basis order at a fixed cutoff can amplify inaccurate oscillatory tails. Convergence checks below show this explicitly; they should not be mistaken for a proof of the limits.
- The long-time loading-2 derivation uses the asymptotically isotropic distant diffusion field. The nonexistence of a decaying two-dimensional steady pressure field is essential to its interpretation.

## Reproduction and Checks

Run **`python3 solve.py`** in this folder. It recomputes the spectral responses, writes `results.json`, and prints the same JSON, including every requested number and general formula. It requires Python, NumPy, and SciPy. No input outside this folder or network access is used.

### Constitutive and Uncoupled Checks

Compared the reduced diffusivity with the supplied expression, and the reduced undrained Lamé constant with its elastic expression. Also removed coupling (`nu_u=nu`) and compared the spectral solution with the exact elastic volume and intensity. Actual output:

```text
no coupling: expected volume 2.5132741228718345 K 1
[[  1.           2.51327412   0.           1.        ]
 [100.           2.51327412   0.           1.        ]]
diffusivity reduced, supplied, residual 1.8703212744369169 1.8703212744369166 2.220446049250313e-16
undrained lambda reduced, expected 4.000000000000001 4.000000000000001
```

These checks could fail through a wrong plane-strain reduction, incorrect trace convention, Fourier normalization, or volume/intensity conversion. The first attempt at the uncoupled check raised `numpy.linalg.LinAlgError: A singular matrix detected: slice(s) [0] are singular.` because `C=0` at zero coupling. The solver now bypasses the irrelevant hydraulic block at exactly zero coupling, and the check above is the successful rerun.

### Drained Endpoint Check

Compared the small-`q` spectral solution with the exact drained volume `pi(1-nu)` and intensity 1. Columns are `q`, volume, finite-`q` short-coefficient estimator, intensity. Actual output:

```text
drained endpoint check: target volume 2.5132741228718345 target K 1
[[1.00000000e-03 2.51327405e+00 7.08981459e-04 1.00000902e+00]
 [1.00000000e-02 2.51320449e+00 7.08902971e-03 9.99982379e-01]
 [1.00000000e-01 2.50955574e+00 7.04785799e-02 9.98632704e-01]
 [1.00000000e+00 2.40115554e+00 5.82469271e-01 9.62898516e-01]
 [4.00000000e+00 2.16264129e+00 1.25333901e+00 8.98053262e-01]]
```

This could fail through incorrect symbols or step normalization. The last few digits of tip intensity are affected by finite cutoff and basis order; the third column is not an estimate of short-time growth in this small-`q` regime.

### Basis and Quadrature Checks

Compared 48, 96, and 160 basis functions at the same cutoff and spacing. Columns here are `q`, volume, short-coefficient estimator. Actual output:

```text
N 48 cutoff 40000 spacing 0.25
[[1.00000000e+02 1.90398691e+00 2.14745457e+00]
 [2.00000000e+02 1.89483305e+00 2.22910313e+00]
 [4.00000000e+02 1.89002641e+00 2.28872417e+00]
 [8.00000000e+02 1.88753885e+00 2.33191341e+00]
 [1.60000000e+03 1.88626444e+00 2.36300012e+00]
 [3.20000000e+03 1.88561608e+00 2.38489156e+00]]
N 96 cutoff 40000 spacing 0.25
[[1.00000000e+02 1.90398691e+00 2.14745459e+00]
 [2.00000000e+02 1.89483305e+00 2.22910318e+00]
 [4.00000000e+02 1.89002641e+00 2.28872428e+00]
 [8.00000000e+02 1.88753885e+00 2.33191366e+00]
 [1.60000000e+03 1.88626444e+00 2.36300922e+00]
 [3.20000000e+03 1.88561619e+00 2.38528074e+00]]
N 160 cutoff 40000 spacing 0.25
[[1.00000000e+02 1.90398694e+00 2.14745816e+00]
 [2.00000000e+02 1.89483308e+00 2.22911064e+00]
 [4.00000000e+02 1.89002645e+00 2.28873984e+00]
 [8.00000000e+02 1.88753888e+00 2.33194633e+00]
 [1.60000000e+03 1.88626448e+00 2.36307926e+00]
 [3.20000000e+03 1.88561623e+00 2.38543819e+00]]
```

Compared increasing cutoff and halving spacing, including overlapping large-`q` extrapolations. Actual output:

```text
N96 cutoff120000 spacing0.25
[[4.00000000e+02 1.89002641e+00 2.28872427e+00]
 [8.00000000e+02 1.88753885e+00 2.33191386e+00]
 [1.60000000e+03 1.88626444e+00 2.36301146e+00]
 [3.20000000e+03 1.88561619e+00 2.38529989e+00]
 [6.40000000e+03 1.88528809e+00 2.40120801e+00]
 [1.28000000e+04 1.88512260e+00 2.41216531e+00]]
quadratic extrapolations [np.float64(2.440061382113711), np.float64(2.440124224292605), np.float64(2.438870923602375)]
N160 cutoff40000 spacing0.125
[[4.00000000e+02 1.89002635e+00 2.28869635e+00]
 [8.00000000e+02 1.88753883e+00 2.33190129e+00]
 [1.60000000e+03 1.88626446e+00 2.36303308e+00]
 [3.20000000e+03 1.88561622e+00 2.38539118e+00]
 [6.40000000e+03 1.88528812e+00 2.40140257e+00]
 [1.28000000e+04 1.88512262e+00 2.41243307e+00]]
quadratic extrapolations [np.float64(2.4403607284375415), np.float64(2.440645750314374), np.float64(2.43940756974985)]
```

These checks could fail through unresolved tip scales, Bessel recurrence instability, insufficient quadrature resolution, incorrect Fourier tails, or unjustified extrapolation. The small discrepancies justify rounding the coefficient to 2.44, not reporting all extrapolated digits.

### Initial Intensity Checks

Compared tip-amplitude convergence at multiple basis sizes and cutoffs. Actual outputs, now with a fourth intensity column:

```text
N160 cutoff80000 spacing0.125; q volume A_est K_est
[[2.00000000e+01 1.96612632e+00 1.83182723e+00 8.49316966e-01]
 [8.00000000e+01 1.90838893e+00 2.11533542e+00 8.28871280e-01]
 [3.20000000e+02 1.89124648e+00 2.27152346e+00 8.18786487e-01]
 [1.28000000e+03 1.88658546e+00 2.35405562e+00 8.13791696e-01]
 [5.12000000e+03 1.88537043e+00 2.39661756e+00 8.11290820e-01]]
N96 cutoff120000 spacing0.125; q volume A_est K_est
[[2.00000000e+01 1.96612632e+00 1.83182722e+00 8.48039148e-01]
 [8.00000000e+01 1.90838893e+00 2.11533538e+00 8.27622845e-01]
 [3.20000000e+02 1.89124648e+00 2.27152325e+00 8.17549797e-01]
 [1.28000000e+03 1.88658546e+00 2.35405480e+00 8.12561515e-01]
 [5.12000000e+03 1.88537043e+00 2.39661887e+00 8.09654966e-01]]
N48 cutoff120000 spacing0.25; q volume A_est K_est
[[2.00000000e+01 1.96612755e+00 1.83185486e+00 8.48023231e-01]
 [8.00000000e+01 1.90838935e+00 2.11537287e+00 8.27606294e-01]
 [3.20000000e+02 1.89124660e+00 2.27156613e+00 8.17532893e-01]
 [1.28000000e+03 1.88658549e+00 2.35409916e+00 8.12099148e-01]
 [5.12000000e+03 1.88537009e+00 2.39466854e+00 7.96049665e-01]]
K quadratic fits [np.float64(0.8075495387906676), np.float64(0.8064007065492166), np.float64(0.771111773897171)]
N96 cutoff300000 spacing0.5; q volume A_est K_est
[[2.00000000e+01 1.96614780e+00 1.83231188e+00 8.48057276e-01]
 [8.00000000e+01 1.90839633e+00 2.11600285e+00 8.27622935e-01]
 [3.20000000e+02 1.89124861e+00 2.27229015e+00 8.17541152e-01]
 [1.28000000e+03 1.88658603e+00 2.35487652e+00 8.12548956e-01]
 [5.12000000e+03 1.88537057e+00 2.39747269e+00 8.09641770e-01]]
K quadratic fits [np.float64(0.8075496261092631), np.float64(0.8075892226710333), np.float64(0.806460525912287)]
N64 cutoff200000 spacing0.125
[[2.00000000e+01 1.96612632e+00 1.83182722e+00 8.48021371e-01]
 [8.00000000e+01 1.90838893e+00 2.11533538e+00 8.27605448e-01]
 [3.20000000e+02 1.89124648e+00 2.27152325e+00 8.17532534e-01]
 [1.28000000e+03 1.88658546e+00 2.35405483e+00 8.12527418e-01]]
K1_0 fit 0.8075431965629399
V1_short fit 2.4398607431830546
```

This check can fail because intensity weights the high-order coefficients much more heavily than volume. Indeed the 48-basis high-`q` fits fail, and the 160-basis run at insufficient cutoff shows a systematic offset. Resolved lower-`q` fits agree near 0.8076; the quoted 0.808 uses only that precision.

### Production Cutoff Check

Compared the production cutoff 160000 against 320000 with everything else fixed. Actual output:

```text
production convergence data
[[8.00000000e+01 1.90838893e+00 2.11533537e+00 8.27607592e-01]
 [3.20000000e+02 1.89124648e+00 2.27152324e+00 8.17534662e-01]
 [1.28000000e+03 1.88658546e+00 2.35405480e+00 8.12546211e-01]
 [5.12000000e+03 1.88537042e+00 2.39660118e+00 8.08414627e-01]]
K1_0 quadratic fit first3 0.8075897671067087
V1_short cubic fit 2.4400664995666825
N80 cutoff320000 spacing0.125
[[8.00000000e+01 1.90838893e+00 2.11533537e+00 8.27605418e-01]
 [3.20000000e+02 1.89124648e+00 2.27152325e+00 8.17532504e-01]
 [1.28000000e+03 1.88658546e+00 2.35405484e+00 8.12544085e-01]
 [5.12000000e+03 1.88537042e+00 2.39660231e+00 8.08412812e-01]]
K1_0 quadratic fit first3 0.807587689880131
V1_short cubic fit 2.440069858483082
```

This could fail if averaged tails were inadequate or the mixed oscillatory tail were large. Cutoff doubling changes neither reported rounded value. It does not eliminate extrapolation uncertainty.

### Reproduction and Opening Integrals

Compared freshly recomputed JSON printed by the reproduction command with the saved JSON, checked required keys, and independently integrated the drained elliptic openings for loading 1 and the loading-2 crack correction. Actual output:

```text
reproduction stdout equals stored results.json: True
required keys present: True
integrated drained opening, formula, residual: 2.5132741228718385 2.5132741228718345 3.9968028886505635e-15
integrated loading2 correction, formula, residual: -0.8388052385084762 -0.8388052385084748 -1.4432899320127035e-15
reproduction output:
{
  "K1_0": 0.808,
  "K1_inf": 1.0,
  "V1_0": 1.8849555921538759,
  "V1_inf": 2.5132741228718345,
  "V1_short_coeff": 2.44,
  "p_faces_0": 0.7490636704119851,
  "K2_0": 0.0,
  "K2_inf": -0.33375,
  "V2_inf": -0.8388052385084748,
  "formulas": {
    "K1_0": "lim_{q->infinity} sum_{j>=0} (-1)^j*(2*j+1)*u_j(q)/(1-nu); u solves the spectral Schur system specified in ANSWER.md",
    "K1_inf": "1",
    "V1_0": "pi*(1-nu_u)",
    "V1_inf": "pi*(1-nu)",
    "V1_short_coeff": "(2/sqrt(pi))*lim_{q->infinity} q*(pi*u_0(q)-pi*(1-nu_u)); u solves the spectral Schur system specified in ANSWER.md",
    "p_faces_0": "2*(nu_u-nu)/(alpha*(1-2*nu))",
    "K2_0": "0",
    "K2_inf": "-alpha*(1-2*nu)/(2*(1-nu))",
    "V2_inf": "-pi*alpha*(1-2*nu)/2",
    "spectral_definitions": "q=sqrt(s), m=sqrt(k^2+q^2), d=2*(nu_u-nu)/(1-nu), Phi_j=pi*(2*j+1)*J_(2*j+1)(k)/k, Psi_j=pi*J_(2*j)(k), D_ij=pi*(2*i+1)*delta_ij/2",
    "spectral_A": "A_ij=D_ij-(d/pi)*integral_0^infinity k^3*Phi_i*Phi_j/[m*(m+k)] dk",
    "spectral_B": "B_ij=(d/pi)*integral_0^infinity q^2*k*Phi_i*Psi_j/[m*(m+k)] dk",
    "spectral_C": "C_ij=(d/pi)*integral_0^infinity q^2*Psi_i*Psi_j/m dk",
    "spectral_u": "(A-B*C^(-1)*B^T)*u=(1-nu_u)*(pi/2)*e_0; use infinite-basis limit before q->infinity; at d=0 use A=D without hydraulic block",
    "V1_short": "V1(t)=(L^2/G)*[pi*(1-nu_u)+V1_short_coeff*sqrt(c*t/L^2)+o(sqrt(c*t/L^2))]"
  },
  "reproduce_command": "python3 solve.py",
  "numerical_uncertainty": {
    "K1_0": 0.001,
    "V1_short_coeff": 0.01
  },
  "K1_0_outer_undrained": 1.0
}
```

The reproduction check could fail because of nondeterministic rounding, a stale JSON file, omitted keys, or mismatched printed and saved values. It verifies reproducibility, not physical correctness. The opening integrals could fail through a missing factor of two between upper-face displacement and full opening, a wrong sign, or inconsistent volume units; their residuals are at floating-point quadrature precision.

After adding the complete spectral definitions as formula strings in JSON, the final artifact check produced:

```text
final printed JSON equals saved JSON: True
final JSON output is pasted verbatim in ANSWER.md: True
deliverables exist: True
all reported formulas are strings: True
```

This compares the final command output to the final saved JSON and the verbatim block above. It could fail through an outdated pasted output, missing deliverable, invalid JSON type, or reproduction mismatch; it is an artifact-consistency check rather than a physical validation.
