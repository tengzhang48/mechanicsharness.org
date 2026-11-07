# Primary Check Audit

All output blocks below are pasted from terminal outputs, not reconstructed
numbers. Independent finite-crack runs, including every initial experiment and
the timeout, are in `finite_check.md`.

## 1. Initial Factorization Refinement

Compared four frequency grids, the canonical norm identity, imaginary parts at
zero frequency, and the lower-half-plane projection residual.

Command: `OPENBLAS_NUM_THREADS=1 python3 factor_check.py`.
At this stage the script used scale 1 and four grids; its current default
instead runs the scale-comparison check in section 3. The solver function is
unchanged except for its subsequently parameterized map scale.

```text
{"K1_0": 0.8087920541442544, "count": 4096, "factor_norm_error": 8.337032702732117e-05, "imaginary_zero_error": 7.101220300318127e-16, "iterations": 41, "n11": 1.0783894055256724, "n21": 0.4129079724470665, "projection_residual": 9.927991183889148e-15}
{"K1_0": 0.8081877744032346, "count": 16384, "factor_norm_error": 2.171795205585525e-05, "imaginary_zero_error": 2.1618173313190687e-15, "iterations": 41, "n11": 1.0775836992043129, "n21": 0.41493182873158596, "projection_residual": 1.0010621830902275e-14}
{"K1_0": 0.8078891162784099, "count": 65536, "factor_norm_error": 5.539811696486652e-06, "imaginary_zero_error": 2.1728615150922834e-15, "iterations": 41, "n11": 1.077185488371213, "n21": 0.41594506462693004, "projection_residual": 1.0005834971553701e-14}
{"K1_0": 0.8077407312274688, "count": 262144, "factor_norm_error": 1.3988069476056353e-06, "imaginary_zero_error": 3.815573907224162e-15, "iterations": 42, "n11": 1.076987641636625, "n21": 0.41645210036961205, "projection_residual": 5.006966149441424e-15}
```

The K and hydraulic factor entries converge only as the square root of grid
spacing. The norm identity converges faster. This check could fail through a
wrong spectral projection, factor normalization, inadequate grid resolution,
iteration nonconvergence, or a correct norm but wrong orientation of the
factor's first column. A small iterative residual does not bound grid error.

## 2. Extended Scale-1 Refinement

Compared the preceding grids with a further fourfold refinement.
Command: `OPENBLAS_NUM_THREADS=1 python3 factor_check.py`, after adding the
1048576 grid, before parameterizing its scale.

```text
{"K1_0": 0.8087920541442544, "count": 4096, "factor_norm_error": 8.337032702732117e-05, "imaginary_zero_error": 7.101220300318127e-16, "iterations": 41, "n11": 1.0783894055256724, "n21": 0.4129079724470665, "projection_residual": 9.927991183889148e-15}
{"K1_0": 0.8081877744032346, "count": 16384, "factor_norm_error": 2.171795205585525e-05, "imaginary_zero_error": 2.1618173313190687e-15, "iterations": 41, "n11": 1.0775836992043129, "n21": 0.41493182873158596, "projection_residual": 1.0010621830902275e-14}
{"K1_0": 0.8078891162784099, "count": 65536, "factor_norm_error": 5.539811696486652e-06, "imaginary_zero_error": 2.1728615150922834e-15, "iterations": 41, "n11": 1.077185488371213, "n21": 0.41594506462693004, "projection_residual": 1.0005834971553701e-14}
{"K1_0": 0.8077407312274688, "count": 262144, "factor_norm_error": 1.3988069476056353e-06, "imaginary_zero_error": 3.815573907224162e-15, "iterations": 42, "n11": 1.076987641636625, "n21": 0.41645210036961205, "projection_residual": 5.006966149441424e-15}
{"K1_0": 0.8076667837231959, "count": 1048576, "factor_norm_error": 3.514379700586545e-07, "imaginary_zero_error": 4.731371389486066e-15, "iterations": 42, "n11": 1.0768890449642612, "n21": 0.41670573503044667, "projection_residual": 5.022002395156206e-15}
```

This confirms the observed square-root convergence, but could fail for the
same reasons as section 1 or through an asymptotic extrapolation based on too
coarse a grid. It is not an independent validation of the governing symbol.

## 3. Frequency-Map Scale Comparison

Compared map scales 0.01, 0.1, and 1.0, each with three fourfold refinements.
The reported extrapolations remove one N^(-1/2) term. This checks that clustering
frequencies near zero does not change the extrapolated physical result.

Command: `OPENBLAS_NUM_THREADS=1 python3 factor_check.py` (current version).

```text
{"K1_0": 0.8076520040311014, "count": 16384, "factor_norm_error": 1.9616627322527336e-07, "imaginary_zero_error": 8.238148675181965e-16, "iterations": 42, "n11": 1.0768693387081352, "n21": 0.4167564718752548, "projection_residual": 7.420153846023195e-15, "scale": 0.01}
{"K1_0": 0.807622493811575, "count": 65536, "factor_norm_error": 5.2774377534703376e-08, "imaginary_zero_error": 8.275983584824932e-16, "iterations": 42, "n11": 1.0768299917487667, "n21": 0.41685795539735343, "projection_residual": 7.388665881698362e-15, "scale": 0.01}
{"K1_0": 0.8076077454627785, "count": 262144, "factor_norm_error": 1.3660317099706276e-08, "imaginary_zero_error": 9.650755099907097e-16, "iterations": 42, "n11": 1.0768103272837046, "n21": 0.4169087022944126, "projection_residual": 7.419914281809807e-15, "scale": 0.01}
{"extrapolated_A": 2.440219570440159, "extrapolated_K1_0": 0.8075929971139819, "extrapolated_norm_error": -1.952986616871044e-08, "scale": 0.01}
{"K1_0": 0.8077799479050168, "count": 16384, "factor_norm_error": 2.1489512707084657e-06, "imaginary_zero_error": 5.802580511269474e-16, "iterations": 42, "n11": 1.0770399305400224, "n21": 0.4163177516116117, "projection_residual": 5.203795106648773e-15, "scale": 0.1}
{"K1_0": 0.807686353472059, "count": 65536, "factor_norm_error": 5.511525302814846e-07, "imaginary_zero_error": 7.486998220606582e-16, "iterations": 42, "n11": 1.0769151379627453, "n21": 0.41663853651882077, "projection_residual": 5.236421857428607e-15, "scale": 0.1}
{"K1_0": 0.8076396460352899, "count": 262144, "factor_norm_error": 1.395290567529628e-07, "imaginary_zero_error": 1.2267101375420699e-15, "iterations": 42, "n11": 1.0768528613803865, "n21": 0.41679897768500346, "projection_residual": 5.31430825793628e-15, "scale": 0.1}
{"extrapolated_A": 2.4402194914047914, "extrapolated_K1_0": 0.8075929385985208, "extrapolated_norm_error": -2.1285493567724245e-07, "scale": 0.1}
{"K1_0": 0.8081877744032346, "count": 16384, "factor_norm_error": 2.171795205585525e-05, "imaginary_zero_error": 2.1618173313190687e-15, "iterations": 41, "n11": 1.0775836992043129, "n21": 0.41493182873158596, "projection_residual": 1.0010621830902275e-14, "scale": 1.0}
{"K1_0": 0.8078891162784099, "count": 65536, "factor_norm_error": 5.539811696486652e-06, "imaginary_zero_error": 2.1728615150922834e-15, "iterations": 41, "n11": 1.077185488371213, "n21": 0.41594506462693004, "projection_residual": 1.0005834971553701e-14, "scale": 1.0}
{"K1_0": 0.8077407312274688, "count": 262144, "factor_norm_error": 1.3988069476056353e-06, "imaginary_zero_error": 3.815573907224162e-15, "iterations": 42, "n11": 1.076987641636625, "n21": 0.41645210036961205, "projection_residual": 5.006966149441424e-15, "scale": 1.0}
{"extrapolated_A": 2.4402187548800103, "extrapolated_K1_0": 0.8075923461765279, "extrapolated_norm_error": -2.1497406514914275e-06, "scale": 1.0}
```

Failure modes: extreme map scales could underresolve either intermediate
frequencies or infinity; a mistaken frequency-zero prescription could survive
all scales; extrapolation could hide common model error. The largest scale's
coarser extrapolation is visibly less accurate, motivating scale 0.01 for the
final computation.

## 4. Final Calculation And Extrapolation

Compared raw factor grids, linear and quadratic extrapolations, and the
independent norm identity. Also printed every requested result and formula.

Command: `OPENBLAS_NUM_THREADS=1 python3 solve.py --diagnostics`.

```text
{
  "K1_0": 0.8075930027492326,
  "K1_inf": 1.0,
  "K2_0": 0.0,
  "K2_inf": -0.33375,
  "V1_0": 1.8849555921538759,
  "V1_inf": 2.5132741228718345,
  "V1_short_coeff": 2.4402195794425063,
  "V2_inf": -0.8388052385084748,
  "formulas": {
    "K1_0": "n11/(1+eta) (true-tip r->0 first, then t->0+)",
    "K1_0_bulk_alternative": "1 (t->0+ at fixed distance first, then r->0)",
    "K1_inf": "1",
    "K2_0": "0",
    "K2_inf": "-a/(2*(1-nu))",
    "V1_0": "pi*(1-nu_u)",
    "V1_inf": "pi*(1-nu)",
    "V1_short_coeff": "4*d/((1+eta)*sqrt(pi))*(4+pi*sqrt(2/eta)*n21)",
    "V2_inf": "-pi*a/2",
    "definitions": "d=nu_u-nu; a=alpha*(1-2*nu); eta=d/(1-nu_u); n11=N_minus(0;eta)[0,0]; n21=N_minus(0;eta)[1,0]",
    "factor_definition": "N=N_plus*N_minus; N_plus and inverse analytic in Im(xi)>0; N_minus and inverse analytic in Im(xi)<0; N_plus(infinity)=N_minus(infinity)=I; square roots have positive real part",
    "factor_kernel": "k=abs(xi); m=sqrt(1+xi^2); t=k/(m+k); D=1+eta-2*eta*k*t/m; k_plus=sqrt(-i*xi); k_minus=sqrt(i*xi); m_plus=sqrt(1-i*xi); m_minus=sqrt(1+i*xi); N=[[D,sqrt(2*eta)*t/(m_plus*k_plus)],[sqrt(2*eta)*t/(m_minus*k_minus),1]]; N(0)=diag(1+eta,1)",
    "normalization": "K_report=K/sqrt(pi*L); V_report=G*V/L^2; unit step loading",
    "p_faces_0": "2*d/a (undrained outer-limit pressure); actual prescribed face trace is 0"
  },
  "p_faces_0": 0.7490636704119851,
  "reproduce_command": "OPENBLAS_NUM_THREADS=1 python3 solve.py"
}
{
  "extrapolated_n11": 1.0767906703323102,
  "extrapolated_n21": 0.4169594526473146,
  "extrapolated_norm_error": -4.666793618213205e-10,
  "face_only_A": 1.354055000514615,
  "factor_grids": [
    {
      "K1_0": 0.8076224938115749,
      "count": 65536,
      "factor_norm_error": 5.277437731265877e-08,
      "imaginary_zero_error": 8.392822674359131e-16,
      "iterations": 42,
      "n11": 1.0768299917487667,
      "n21": 0.41685795539735343,
      "projection_residual": 7.378395846389504e-15,
      "scale": 0.01
    },
    {
      "K1_0": 0.8076077454627784,
      "count": 262144,
      "factor_norm_error": 1.3660316877661671e-08,
      "imaginary_zero_error": 8.713135287368024e-16,
      "iterations": 42,
      "n11": 1.0768103272837046,
      "n21": 0.41690870229441257,
      "projection_residual": 7.419293514385307e-15,
      "scale": 0.01
    },
    {
      "K1_0": 0.8076003734015992,
      "count": 1048576,
      "factor_norm_error": 3.473425547184661e-09,
      "imaginary_zero_error": 9.401426776692667e-16,
      "iterations": 42,
      "n11": 1.076800497868799,
      "n21": 0.41693407703888324,
      "projection_residual": 7.464662548365556e-15,
      "scale": 0.01
    }
  ],
  "linear_extrapolated_A": 2.4402195771919195,
  "linear_extrapolated_K1_0": 0.80759300134042
}
```

Failure modes: a wrong assumed expansion in N^(-1/2), a shared derivation error,
or factor norm agreement concealing orientation error. These checks establish
numerical stability rather than a rigorous error bound. Independent finite
crack numerics address the shared-factorization risk at lower precision.

## 5. Reproduction, Schema, Energy, And Weak Coupling

Compared stored JSON against stdout and regenerated JSON, verified required
number/string types, compared A with the energy-derived alternative using K,
and checked the first-order weak-coupling intensity.

Command: `OPENBLAS_NUM_THREADS=1 python3 check_results.py`.

```text
Reproduction: stored JSON == command stdout JSON == regenerated JSON: PASS
Schema: all 9 Report quantities are numbers; all formulas are strings: PASS
Volume from independently derived energy identity: 2.440219580900
Reported volume coefficient: 2.440219579443
Difference: 1.458e-09
Weak coupling eta=0.0001: K=0.999925013475; first-order K=0.999925000000
Weak-coupling scaled deficit (1-K)/eta: 0.749865248; expected limit 0.75
Energy relation and weak-coupling checks: PASS
```

Failure modes: nondeterminism or mismatched stdout would fail equality; missing
or wrongly typed quantities would fail the schema assertions; wrong hydraulic
normalization or tip signs would usually fail the energy/weak-coupling checks.
Equality alone cannot validate the mathematics. The weak-coupling comparison
has both finite eta and discretization error, and does not establish finite-eta
accuracy by itself.
