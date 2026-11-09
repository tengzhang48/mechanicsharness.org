# Grade

| rubric item | points awarded / max | reason |
| --- | --- | --- |
| Stress exponent | 15 / 15 | **10**: reported p = 0.1333333333 differs from 0.133333333333 by 3.3×10⁻¹¹, below 10⁻⁵. **5**: `solve.py` derives W ∝ r^((n+1)(s−2)) = r⁻¹, gives the contour power (n+1)(s−2)+1 = 0, and explains that both contour terms scale as W r; `final.py` demonstrates a positive, radius-independent flux. This supplies the energy/flux scaling connection beyond simply quoting p. |
| Mechanics and boundary conditions | 25 / 25 | **10**: the Airy stress components satisfy equilibrium, and `make_rhs` implements and solves the fourth-order strain-compatibility equation. **10**: the effective stress √[¾(σrr−σθθ)²+3σrθ²] implements σzz = (σrr+σθθ)/2; the normal-strain factor ¾ and tensor-shear factor 3/2 are correct, with no added elastic term. **5**: f′(0)=f‴(0)=0, f(π)=f′(π)=0, correct even/odd continuation, and explicit normalization max σ̃e=1 are implemented. |
| Tensile coefficient H | 25 / 25 | **20**: reported H = 1.8935697252 differs from 1.893569725213 by 1.3×10⁻¹¹, below 0.005; rerunning `solve.py` gives 1.893569725217. **5**: the answer and `solve.py` correctly calibrate the dimensional stress amplitude K using J = αε₀ K^(n+1) I/σ₀ⁿ and obtain H = σ̃θθ(0) I⁻ᵖ. The dimensional reconstruction returns the same H. |
| Energy coefficient C | 20 / 20 | **15**: reported C = 0.4573057839 differs from 0.457305783974 by 7.4×10⁻¹¹, below 0.002; rerunning `solve.py` gives 0.457305783939. **5**: W ∝ 1/r makes the radial integral proportional to R, and the full angular integral gives C = [n/(n+1)] shape2/I. U is energy per unit thickness, with no extra thickness factor. |
| Verification | 15 / 15 | **10**: the Cartesian finite-difference contour calculation gives J ≈ 7.70940655×10⁴ J/m² at R = 3×10⁻⁴, 10⁻³, 3×10⁻³, and 10⁻² m; H also agrees across amplitudes 0.5, 1, and 2. The separate elastic contour check reproduces 4.550000×10⁻¹² at all three stated radii. **5**: reproduced endpoint residuals are f(π)=1.990×10⁻¹⁵ and f′(π)=8.251×10⁻¹⁷, before effective-stress normalization. These successful checks satisfy the rubric despite the additional failed checks described below. |

**Total: 100 / 100.**

**Scientific errors**

- `final.py:280` and `final.py:314` confuse its amplitude A, defined by σ = σ₀ A r⁻ᵖ σ̃, with the dimensional stress amplitude K = σ₀ A. Its closed expressions must use αε₀σ₀ A^(n+1), rather than αε₀ A^(n+1)/σ₀ⁿ. The missing factor is σ₀^(n+1) = 3.2768×10⁶⁴. Consequently its printed J_closed = 2.352724×10⁻⁶⁰ disagrees with J_numeric = 7.709407×10⁴, and U_closed = 1.075914×10⁻⁶³ disagrees with U_numeric = 35.25556. The answer and `solve.py` use K consistently; the scored numerical H and C do not depend on these faulty auxiliary expressions.
- `final.py:252` assigns the variable used for shear strain to −em instead of exy. Its equivalent-strain check therefore prints ratios 0.997466, 0.960337, and 0.956520 instead of unity. This is a checker error, not evidence that the computed field violates the constitutive law: auditing the same displacement functions with the actual shear strain gives a maximum relative tensor error of 1.63×10⁻⁸ at those points. The separate angular check in `check_constit.py` also reproduces the stated normalization ratio 50.
- Some preliminary checks remain incorrect: `verify.py` omits the odd-component sign reversal and changes normalization with the evaluation sample; its supposed eigenvalue update also ignores the trial s argument. `check_constit.py`'s final Jhat calculation treats polar stresses as Cartesian components and extrapolates the IVP outside its solved angular interval, returning −0.02232011 instead of approximately 4.81837909. These checks do not support a positive contour flux. The final contour implementation fixes the reflection and normalization problems.

**Unsupported claims**

- The successful closed-form-versus-numerical dimensional comparisons advertised by `final.py` are contradicted by its own output. Agreement is obtained with the correct prefactor, but its printed comparisons do not establish it.
- The claim of numerical amplitude invariance for C at A ∈ {0.5, 1, 2} is not demonstrated by the supplied loops: they vary amplitude only for H. C's amplitude cancellation is supported analytically, and its radius invariance is demonstrated numerically.

**Presentation problems**

- The displacement ansatz in the answer stops at `u_i = (alpha eps0/… )`; dimensional versus angular quantities and the raw versus normalized f(0) are not always clearly distinguished. The code resolves these ambiguities.
- Several exploratory scripts and debugging comparisons are retained without clearly identifying which checks remain valid. This obscures an otherwise correct final calculation.

Verification evidence is retained in `scratch/solve_output.txt`, `scratch/final_output.txt`, `scratch/elastic_output.txt`, `scratch/constitutive_output.txt`, and `scratch/audit_output.txt`. The `final.py` rerun changes only the radial sample count in U_numeric from 400 to 7 to limit memory; W r is independent of r, so this preserves the radial quadrature. The copied reference rerun reproduces its recorded results; no reference error was found, so no alternative corrected-reference grade is needed. The rubric provides no additional deduction for failed supplementary checks when the required successful checks and independently correct calculations are demonstrated.

**Tolerance verdict:** all requested numbers are within tolerance: p error 3.3×10⁻¹¹ < 10⁻⁵; H error 1.3×10⁻¹¹ < 0.005; C error 7.4×10⁻¹¹ < 0.002.
