# Grade

| rubric item | points awarded / max | reason |
|---|---:|---|
| Elastic strip bending coupled to closing face forces | 7 / 7 | Uses both Euler–Bernoulli arms and the local closing traction. The implemented equation for u = 1 − δ/δc, u'''' = λ⁴u, is equivalent to δ'''' = −Ct(δ). |
| Individual bending stiffnesses and combined opening compliance | 10 / 10 | Correct individual values d₁ = 1.98 and d₂ = 11.5425 N·m earn 5/5. The harmonic stiffness d̄ = (1/d₁ + 1/d₂)⁻¹ = 1.69008319 N·m correctly adds opening compliances, earning 5/5. |
| Local softening curve and triangular connector work | 7 / 7 | Applies t = σb(1 − δ/δc) locally to the full relative gap. Correctly obtains connector work σbδc/2 = 20 J/m². |
| Front displacement/slope matching and intrinsic work | 8 / 8 | Enforces δ(0) = δ′(0) = 0 (4/4) and the nonzero intrinsic-work curvature δ″(0) = √(2Γ₀/d̄) (4/4). Both are implemented in the exact solution coefficients. |
| Zero rear shear and moment/width convention | 8 / 8 | Enforces u'''(L) = 0 for pure couples, computes m = d̄δ″(L), and converts to total moment on each strip with M = bm. Correctly excludes the starter-cut length from the retained-zone calculation. |
| M on each strip | 15 / 15 | Reported 0.2126 N·m has relative error 0.0121281% against the unrounded key, within 0.5%. The parenthetical 0.21263 N·m also qualifies. |
| L | 15 / 15 | Reported 10.80 mm has relative error 0.000947870%, within 0.5%. |
| Full opening δ(L/2) | 10 / 10 | Reported 38.2 µm has relative error 0.0494465%, within 0.5%. The parenthetical 38.22 µm also qualifies. |
| Verification | 12 / 12 | Executed code obtains the moment from rear curvature and independently checks m²/(2d̄) = Γ₀ + σbδc/2 = 21.4 J/m², with relative discrepancy about 1.5 × 10⁻¹⁵. The traction-moment quadrature also agrees within 8.0 × 10⁻¹⁵ N. This satisfies the full-credit energy-balance criterion. |
| Scope and event recognition | 8 / 8 | Uses the prescribed elastic, small-rotation regime and identifies first failure at the original root, imposing δ(L) = δc and checking the spatial opening is monotone (4/4). The length and profile are solved using both arm stiffnesses and beam boundary conditions, rather than the curve area alone (4/4). |

**Total: 100 / 100.**

The submitted code was copied to `scratch/solve.py` and executed successfully. Its unrounded results are M = 0.2126257873767215 N·m, L = 10.800102370897184 mm, and δ(L/2) = 38.21889792609383 µm. Each agrees with `reference_results.json` to better than 1.5 × 10⁻¹²% relative error. Evidence is retained in `scratch/solve_output.txt` and `scratch/audit_results.json`. No reference error was found; alternative grading is unnecessary. The ancillary errors below do not invalidate the scored formulation or the executed energy balance.

## (a) Scientific errors

- The claimed Γ₀ → 0 limit s → π/2 is wrong. At fixed connector parameters, k → 0 and the characteristic equation becomes 1 + cosh(s)cos(s) = 0, whose first admissible root is s = 1.875104068712. Substitution of π/2 leaves a residual of 1, not zero.
- Calling k → 0 “no bridging” while σbδc → ∞ reverses its physical meaning: k² = 2Γ₀/(σbδc) measures intrinsic work relative to connector work. Also, k → 0 alone does not establish L → 0; L = s(d̄δc/σb)¹/⁴ depends on the parameter path.

## (b) Unsupported claims

- The “sensible limits” bullet is presented as a check, but the submitted script executes no limit calculations; its stated limits are incorrect as described above.
- The finite-difference ODE diagnostic does not establish a small numerical residual. Its h = 10⁻⁶ m stencil returns about 5.21 × 10⁸ m⁻⁴, approximately 109% of the ODE right-hand side. Cancellation makes this diagnostic unreliable: changing h to 2 × 10⁻⁵ m reduces the relative residual to about 4.57 × 10⁻⁷. The analytical solution satisfies the ODE, and the separate energy check remains valid.

## (c) Presentation problems

- The displayed energy check writes m²/2d̄ without grouping the denominator; write m²/(2d̄). The code uses the correct expression.
- The starter-cut interval is written L < x < a without defining a. With the stated front-based coordinate, the free ends are at x = L + 20 mm.
- The code initially labels the dimensionless u(L/2) value as δ(L/2), before displaying the correctly converted opening after an arrow.

**Tolerance:** All three requested numbers, including the rounded final quoted values, are within the required ±0.5% of the unrounded reference.
