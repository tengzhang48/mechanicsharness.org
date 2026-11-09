| rubric item | points awarded / max | reason |
|---|---:|---|
| Elastic strip bending coupled to closing face forces | 7 / 7 | Implements Euler–Bernoulli bending with δ'''' = −bR t(δ), rather than treating the retained faces as traction-free. |
| Individual stiffnesses and opening-compliance coupling | 10 / 10 | EI₁ = 0.0495 N·m² and EI₂ = 0.2885625 N·m² are correct (5 points); R = 1/EI₁ + 1/EI₂ = 23.6674739599 (N·m²)⁻¹ is correct (5 points). |
| Local softening law and triangular connector work | 7 / 7 | Uses the full relative gap in t = σ_b(1 − δ/δ_c) and correctly obtains connector work 20.0 J/m². |
| Front displacement/slope matching and intrinsic work | 8 / 8 | Under the coordinate change s = L − x, the implemented δ(s=L) = δ'(s=L) = 0 supplies matching (4 points). The selected branch has front curvature 1.287137591 m⁻¹ = √(2CΓ₀), retaining Γ₀ = 1.40 J/m² (4 points), through the equivalent energy formulation. |
| Zero rear shear and correct moment/width convention | 8 / 8 | At the starter-cut root, the implemented δ'''(s=0) = 0 and δ''(s=0) = MR correctly match pure couples. R uses full-width stiffnesses; M is the total moment on each arm. |
| M, total on each strip | 15 / 15 | Reported 0.2126 N·m versus key 0.21262578737672128 N·m: relative error 0.0121281%. The more precise reported 0.21262579 N·m and rerun result also pass. |
| L | 15 / 15 | Reported 10.80 mm versus key 10.800102370897305 mm: relative error 0.000947870%. The more precise reported 10.8001 mm and rerun result also pass. |
| Full relative opening δ(L/2) | 10 / 10 | Reported 38.22 µm versus key 38.218897926093284 µm: relative error 0.00288358%. The more precise reported 38.2189 µm and rerun result also pass. |
| Verification | 6 / 12 | Executed boundary residual checks qualify for 6. The claimed energy check sets M from the same imposed balance and prints Γ₀ + σ_bδ_c/2 as “supply”; it does not execute a supply-versus-work check for the computed solution. The compliance check only confirms the algebraic coefficient. No independent solution, limit, or convergence check is demonstrated; the separate alternative is invalid. |
| Prescribed elastic regime and first connector-failure event | 4 / 4 | Uses elastic beam bending without stretching or shear. The selected profile stays between 0 and δ_c and first reaches the cutoff at the starter-cut root; no post-failure solution is used. |
| Dependence of zone/profile on both arms and geometry | 4 / 4 | Uses both stiffnesses in the local beam equation and solves for the length/profile separately from the connector-work area. |

**Total: 94 / 100** (mechanics 40/40; numerical outputs 40/40; verification 6/12; scope 8/8).

The endpoint labels in the answer contradict the requested coordinate, but the main implemented solution is mechanically equivalent after s = L − x. The grader's audit finds a maximum mirrored-profile discrepancy of 4.30×10⁻⁹ µm and implied intrinsic front work 1.40000000022 J/m². These grader-added checks establish modeling equivalence; they do not earn submission verification points. The unchanged reference program reproduces the recorded key; no reference error was found.

**(a) Scientific errors**

- The prose puts δ = δ_c at the advancing front and δ = 0 at the starter-cut root. Physically these endpoints are reversed: the faces meet at the advancing front, and the oldest retained connector fails at the starter-cut root. The correct main calculation is accepted in its reversed coordinate.
- The unused alternative in solve.py attempts zero curvature at the bonded front, inconsistent with nonzero Γ₀. The separate solve2.py encodes the cutoff condition incorrectly and returns negative gaps: midpoint −198.7083 µm and rear −472.2804 µm. It provides no valid alternative mechanics check.

**(b) Unsupported claims**

- The stated energy verification restates the equation used to determine M. In verify.py, the printed “G supply” is the prescribed resistance, rather than energy release evaluated from the computed M.
- A sign-changing bracket does not establish uniqueness or rule out other roots. The scalar equation also has roots at approximately 13.2446, 48.7057, and 52.0199 mm. Those branches violate the permitted gap range, but the answer does not demonstrate that exclusion.
- Calling the length-independent growth criterion “stable” is unsupported by a stability analysis.

**(c) Presentation problems**

- The coordinate reversal is never declared, making the stated front/root boundary conditions misleading.
- solve.py prints the correct main results, then exits with a TypeError because its appended fsolve calculation returns four residuals for three unknowns. This late failure does not erase the completed main calculation.
- verify.py is described as checking the ODE but contains no ODE residual evaluation. check.py's finite-difference ODE residual is 4.151×10⁴ m⁻³, not zero; the answer does not explain the numerical cancellation in that check. The analytic main profile nevertheless satisfies the ODE.

**Tolerance:** All three requested numbers, including the rounded results table, are within ±0.5% of the unrounded reference.
