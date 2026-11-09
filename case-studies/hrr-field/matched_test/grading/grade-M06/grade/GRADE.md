| rubric item | points awarded / max | reason |
|---|---:|---|
| Model: compatibility, equilibrium, and nonlinear shear law | 8 / 8 | Uses engineering shear γ = ∇w, div τ = 0, and the correct inverse power law, giving the p-Laplacian with exponent 9/7. The submitted ODE implements this formulation. |
| Model: admissible radial behavior | 8 / 8 | Derives 1 + (s − 1)(n + 1)/n = 0 from finite, nonzero, radius-independent J, obtaining s = 2/9 and stresses proportional to r^(−2/9). |
| Model: symmetry, free faces, and positive slip | 6 / 6 | Implements odd f and even f′, the free-face condition f′(π) = 0, and positive amplitude/slip. Rerun gives f′(π) = 1.0842 × 10^−16. All three conventions earn 2 points each. |
| Model: magnitude fixed by J | 8 / 8 | Uses stored work W = n/(n + 1) τ·γ and both terms of the supplied flux. Rerun gives I_J = 2.019595277307725 and A = 0.0009256090617515118, matching the key. |
| Δw at 0.0200 mm | 8 / 8 | Answer: +0.00131813 mm; key: +0.001318126778 mm. Relative error 0.00024446%, within 2%. |
| Δw at 0.0800 mm | 8 / 8 | Answer: +0.00179369 mm; key: +0.001793693738 mm. Relative error 0.00020840%, within 2%. |
| Signed τ_rz at 0.0200 mm, 60° | 6 / 6 | Answer: +38.587 MPa; key: +38.58659686 MPa. Relative error 0.00104476%, within 2%, with the correct sign. |
| Signed τ_θz at 0.0200 mm, 60° | 6 / 6 | Answer: +138.011 MPa; key: +138.0106166 MPa. Relative error 0.00027778%, within 2%, with the correct sign. |
| Engineering-shear contour radius at 0° | 4 / 4 | Answer: 0.0468944 mm; key: 0.04689442305 mm. Relative error 0.00004916%, within 2%. |
| Engineering-shear contour radius at 90° | 4 / 4 | Answer: 0.0250661 mm; key: 0.02506612349 mm. Relative error 0.00009371%, within 2%. |
| Engineering-shear contour radius at 180° | 4 / 4 | Answer: 0.0133984 mm; key: 0.01339840659 mm. Relative error 0.00004916%, within 2%. |
| Completion: seven identifiable numerical outputs | 5 / 5 | All seven requested values are present and reproduced by the submitted code. |
| Conventions: units, engineering shear, and face limit | 5 / 5 | Correct mm/MPa units (2), γ = ∇w rather than tensor shear (2), and explicit θ → π⁻ face-limit convention (1). |
| Verification: genuine requirement | 5 / 5 | Code checks the dimensional J-integral using Cartesian stresses and a finite-difference displacement derivative; it also checks Cartesian equilibrium and face traction. These address genuine field requirements. |
| Verification: numerical result | 5 / 5 | Rerun gives J = 0.120000000000820, 0.120000000000842, and 0.120000000000880 N/mm on the three circles; max divergence residual is 1.4211 × 10^−8 MPa/mm and face shear is 1.0787 × 10^−11 MPa. |
| Verification: interpretation | 5 / 5 | Correctly interprets agreement with J = 0.120 N/mm as amplitude/path-independence verification and the small face traction as satisfying the free-face condition. These supported checks suffice for all 15 points. |
| Recognition of missing fracture physics | 5 / 5 | States that initiation requires a fracture criterion, identifying critical J or a process-zone stress/separation criterion. The stationary field alone does not determine initiation. |

**Total: 100 / 100.** The reference rerun passes its independent-route and linear-limit assertions and reproduces the recorded key; no reference correction or alternative score is warranted.

**(a) Scientific errors**

- `submission/ANSWER.md:49` incorrectly reports f′(π/2) = 0.7636948. That is f′(π/3); the submitted code gives f′(π/2) = 0.54180841765 and g(π/2) = 0.61435299915. The stated derivative is also incompatible with the accompanying g because g ≥ |f′|. The code and requested 90° radius are correct. This extra intermediate-value error does not remove credit from any fulfilled rubric item.

**(b) Unsupported claims**

- The claim that the submitted solver numerically returns the n = 1 limiting solution is not demonstrated: `submission/solve.py` fixes n = 3.5 and contains no linear-limit run. The analytic limiting formulas are valid with their respective normalizations, but a claimed numerical limiting test is not credited.
- The claimed max-step/tolerance refinement study is not shown. The submitted code contains only one integration setting and no refinement results. Verification credit rests on the reproducible flux and residual checks above.

**(c) Presentation problems**

- The linear-limit paragraph switches between f′(0) = 1 and f′(0) = 1/2 without clearly separating the corresponding values. Under the first normalization, f(π) = 2 and I_J = π; under the second, f(π) = 1 and I_J = π/4.
- The incorrect angular label in the intermediate-value line makes the written angular data internally inconsistent, even though the calculation uses the correct values.

**Tolerance:** All seven requested signed/numerical outputs are within the stated 2% relative tolerance; the largest relative error is 0.00104476%.
