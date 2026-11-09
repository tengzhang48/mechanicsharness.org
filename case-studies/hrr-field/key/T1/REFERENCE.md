# Reference

The requested values are

| Quantity | Reference | Absolute tolerance |
| --- | ---: | ---: |
| p | 0.133333333333 | 0.00001 |
| H | 1.893569725213 | 0.005 |
| C | 0.457305783974 | 0.002 |

[reference.py](reference.py) is standalone NumPy/SciPy code. Run `python T1/reference.py` from the build folder; it writes [reference_results.json](reference_results.json). The numerical values above concern the stipulated limiting constitutive law, rather than stresses at a finite physical distance where elasticity or blunting might matter.

## Derivation and normalization

Write Φ = A rˢ f(θ), normalize f(0) = 1, and omit the arbitrary amplitude A during the dimensionless calculation. With primes denoting θ derivatives, put

\[
b_r=sf+f'',\quad b_\theta=s(s-1)f,\quad b_{r\theta}=(1-s)f',
\quad D=f''+s(2-s)f,\quad
Q^2=\tfrac34D^2+3(s-1)^2f'^2.
\]

The plane-strain closure is σzz = (σrr+σθθ)/2, so the effective stress involves the **difference** of the in-plane normal stresses, not their hydrostatic mean. The dimensionless strain amplitudes are eᵣ = 3Qⁿ⁻¹D/4, eθ = −eᵣ, and eᵣθ = 3Qⁿ⁻¹bᵣθ/2.

Compatibility gives the angular equation

\[
(Q^{n-1}D)''-k(k+2)Q^{n-1}D
 +4(s-1)(k+1)(Q^{n-1}f')'=0,\qquad k=n(s-2).
\]

Impose f′(0) = f‴(0) = 0 and f(π) = f′(π) = 0, plus the normalization. The code solves for s as an unknown eigenvalue. A nonzero, radius-independent J requires (n+1)(s−2)+1 = 0, independently predicting s = 1.866666666667 and p = 1/(n+1).

For displacement exponent a = k+1, write uᵣ = rᵃG and uθ = rᵃL. Then G = eᵣ/a, L = (2eᵣθ−G′)/(a−1), and L′ = eθ−G. At r = 1,

\[
(u_{,x})_r=aG\cos\theta-(G'-L)\sin\theta,\qquad
(u_{,x})_\theta=aL\cos\theta-(L'+G)\sin\theta.
\]

Integrating the stated contour with W = nQⁿ⁺¹/(n+1) gives J₀ = 0.307101353930 in the f(0) = 1 convention. Hence H = s(s−1)/J₀¹⁄⁽ⁿ⁺¹⁾. Since W ∝ 1/r, radial integration gives C = [2n/(n+1) ∫₀^π Qⁿ⁺¹dθ]/J₀. The arbitrary amplitude cancels in both quantities. The calculated initial curvature is f″(0) = −0.689256761930.

## Independent checks actually run

* **Different numerical method:** adaptive collocation solves the five-condition eigenproblem; DOP853 IVP shooting independently adjusts s and f″(0) to satisfy both free-face conditions. Absolute differences are 3.69×10⁻¹⁴ in p, 3.53×10⁻¹³ in H, and 4.45×10⁻¹³ in C. The shooting traction residual is below 5×10⁻¹⁵ in the chosen normalization. Both methods use the derived angular equation, so this comparison alone checks its numerical solution rather than its derivation.
* **Different evaluation of J:** reconstruct Cartesian displacements, take a centered numerical x derivative, and integrate the Cartesian traction flux on dimensionless radii 0.6, 1.0, and 1.7. The maximum relative difference from the analytic polar derivative is 3.24×10⁻¹². This checks basis rotation, energy normalization, and radius independence.
* **Analytic special solution:** at n = 1 the given law is incompressible linear elasticity. Its exact opening field gives p = 1/2, H = √[2/(3π)], and C = 1/4. The largest absolute error across these quantities is 1.63×10⁻¹³. This limit tests the constitutive, shear, and energy factors.
* **Source corroboration:** the supplied paper's plane-strain n = 3 table gives I ≈ 5.51 when max Q = 1 and its energy coefficient Eₙ ≈ 1.23. Our code obtains 5.507323232908 and πC = 1.227162152550, consistent with that rounding. The n = 6.5 reference was solved directly, with no interpolation of the table.

The very close solver agreement is an internal numerical check, not a claim of comparable physical accuracy for a real blunted elastic-plastic crack.
