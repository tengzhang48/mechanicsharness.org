# Reference

| Quantity | Reference | Absolute tolerance |
| --- | ---: | ---: |
| p | 0.156175260286 | 0.0003 |
| R | 0.373262313312 | 0.002 |

Run the standalone [reference.py](reference.py) with `python T2/reference.py`. It writes all checks to [reference_results.json](reference_results.json). Its dependencies are NumPy and SciPy.

## Model adaptation

For a singular field and n > 1, the power-law strains grow faster than the linear elastic strains as r tends to zero. Thus the **inner** leading field uses the nonlinear part of the stated constitutive law. The elastic outer solution remains useful for matching, but its notch exponent does not determine this limiting field. Monotonic proportional loading permits the stipulated total-deformation description. This is an ideal asymptote; a finite notch radius and large deformations would cut it off physically.

Let Φ = A rˢf(θ), h = s−1, k = n(s−2), D = f″+s(2−s)f, and Q² = 3D²/4 + 3h²f′². Plastic plane strain gives σzz = (σrr+σθθ)/2. As obtained from equilibrium and strain compatibility,

\[
(Q^{n-1}D)''-k(k+2)Q^{n-1}D
 +4h(k+1)(Q^{n-1}f')'=0.
\]

**The geometry changes the eigenproblem:** impose f′(0) = f‴(0) = 0 and f(β) = f′(β) = 0 at β = 2π/3, not π. With f(0) = 1, the two unknowns are s and f″(0). Choose the leading symmetric tensile branch, with positive f in the interior. Finite energy requires p = 2−s < 2/(n+1).

The solution is s = 1.843824739714 and f″(0) = −1.263078996343. Therefore

\[
p=2-s,\qquad
R=\frac{s+f''(0)}{s(s-1)}.
\]

The radial power of energy inside radius r is 2−(n+1)p = 1.219123698570, rather than the crack value 1. The crack's nonzero constant-J scaling does not fix a notch exponent: inclined traction-free notch faces still contribute the W Nₓ term to a closed energy-flux contour. No absolute stress amplitude can be obtained from this local eigenproblem alone. E, ν, α, and σ₀ set matching/scales but cancel from the two requested inner-field quantities.

## Independent checks actually run

The primary solver uses **mixed variables** (f, v=f′, g=Qⁿ⁻¹D, g′), inverts the monotone constitutive relation for D at each evaluation, and applies adaptive boundary-value collocation. Compatibility becomes g″ = k(k+2)g−4h(k+1)(Qⁿ⁻¹v)′.

The second route explicitly differentiates the constitutive relation to form a fourth-order ODE in (f,f′,f″,f‴), then uses DOP853 shooting and a two-parameter root solve. It does not use the primary solver's constitutive-inversion code. The two methods differ by 2.58×10⁻¹⁴ in p, 3.14×10⁻¹⁴ in R, and at most 1.25×10⁻¹² in f over 5,001 angles. The BVP residual is below 10⁻⁹; the shooting face residual is below 1.8×10⁻¹⁵.

Additional independent limits check the equation, boundary-angle convention, and branch:

* **Linear limit n = 1:** the symmetric traction-free elastic wedge obeys sin(2βλ)+λ sin(2β) = 0, λ=s−1. The analytic root gives p = 0.384268940509 and R = 0.444272010788. The mixed solver agrees to 3.27×10⁻¹³ and 6.36×10⁻¹³ respectively. R is independently calculated from the analytic angular function A cos(sθ)+B cos[(s−2)θ]. This linear result is not the answer at n = 4.
* **Zero-opening limit β = π, n = 4:** the solver gives p = 0.200000000000098, agreeing with the crack's exact 1/(n+1) to 9.8×10⁻¹⁴.
* **Rounded source special case:** a 90° opening at n = 3 gives p = 0.225016183019, consistent with the supplied paper's reported 0.225. This is corroboration at another geometry/material parameter, not a source for the new reference values.

The agreement verifies the idealized eigenfield, not the physical extent of its dominance in a finite specimen.
