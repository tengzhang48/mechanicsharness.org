# Grade

| rubric item | points awarded / max | reason |
|---|---:|---|
| Elastic bending coupled to closing face forces | 7 / 7 | Uses two Euler–Bernoulli arms and the local closing traction; solves the bridged opening equation rather than assigning the connector work as an isolated-crack toughness. |
| Individual stiffnesses and opening-compliance coupling | 10 / 10 | Correct total-width stiffnesses: 0.04950 and 0.2885625 N·m² (5 points). Correct sum of inverse stiffnesses, C_M = 23.6674739599 (N·m²)⁻¹ (5 points). |
| Local softening law and triangular work | 7 / 7 | Implements δ'''' = b C_M σ_b(δ/δ_c − 1) for the event profile and integrates W(δ) = σ_b δ − σ_b δ²/(2δ_c), giving W(δ_c) = 20.0 J/m². |
| Front displacement/slope matching and intrinsic work | 8 / 8 | Enforces δ(0) = δ'(0) = 0 (4 points). The implemented growth relation retains Γ_0 = 1.40 J/m² and is equivalent to the front-energy condition (4 points): the event solution gives δ''(0) = 1.2871375906 m⁻¹ = √(2b C_M Γ_0). |
| Rear shear, curvature, and moment convention | 8 / 8 | Enforces δ'''(L) = 0 and δ''(L) = C_M M with total-width stiffnesses. The energy formula M = √[2b(Γ_0 + W)/C_M] returns the total couple on each strip. |
| Numerical M on each arm | 15 / 15 | Reported 0.212626 N·m; relative error 0.0001000% against the unrounded key. Rerun: 0.2126257873767213 N·m. |
| Numerical L | 15 / 15 | Reported 10.8001 mm; relative error 0.000021953% against the unrounded key. Rerun: 10.8001023708972 mm. |
| Numerical full opening δ(L/2) | 10 / 10 | Reported 38.219 µm; relative error 0.000267077% against the unrounded key. Rerun: 38.2188979260938 µm. |
| Verification | 12 / 12 | Executed potential-energy differentiation at the event gives −∂Π/∂L = bΓ_0 = 0.035 J/m, with relative discrepancy 2.69×10⁻⁶. Independent solve_ivp shooting reproduces the midpoint and rear openings. These valid checks meet the full-credit criterion; flawed supplementary checks listed below earn no additional credit. |
| Scope: prescribed elastic regime and first failure | 4 / 4 | Uses bending-only, small-rotation mechanics and identifies the first inactive connector at the original cut root. The event profile increases from zero to 200 µm without an earlier interior cutoff. Continuation on the admissible branch reaches this event from initiation. |
| Scope: both arms and geometry determine the profile | 4 / 4 | Uses both arm compliances and a spatial boundary-value solution to determine L and the midpoint opening; uses the curve area to determine the event moment. |

**Total: 100 / 100.**

The reference rerun reproduces the recorded key and its independent checks; no reference correction is warranted. The required event calculation satisfies every rubric row. The supplementary faults below do not invalidate that calculation or its successful verification, and the rubric provides no separate general deduction for them. Rerun evidence is in [scratch/solve.log](scratch/solve.log), [scratch/verify.log](scratch/verify.log), and [scratch/audit.json](scratch/audit.json).

**(a) Scientific errors**

- ANSWER.md says the negative foundation stiffness produces “exponential rather than oscillatory modes.” The characteristic roots are ±μ and ±iμ, so both hyperbolic and trigonometric modes occur. The implemented basis is correct.
- The supplementary `solve.py` function `M_of_L` selects the smallest positive algebraic moment instead of the admissible continuation branch. At the correct event length it prints M = 0.18480559 N·m and rear opening 97.69023 µm. That solution has negative front curvature and a minimum gap of −1.24456 µm, so it cannot describe the prescribed opening path. The main event solve and `verify.py` select the correct branch.
- Supplementary diagnostics use the linear traction law outside its prescribed opening range. `verify.py` continues to L = 12 mm and δ_t = 251.949 µm without making failed connectors inactive; its declining moment beyond the event is therefore not a physical continuation. The extra energy sample at L = 14.0401 mm also has negative gaps (minimum −8.08559 µm), outside the stated traction law's domain.

**(b) Unsupported claims**

- The printed “σ_b→0 check” does not decrease σ_b or recompute μ. It compares the original finite-traction result, 20.864886 µm, with the traction-free result, 29.584342 µm. The analytical zero-traction limit is correct, but this code does not execute that limit check.
- Positivity of the quoted derivative G_M alone does not establish stability under moment control. For F(M,L) = J_far − Γ_0 − W(δ_t), the argument also needs the sign of F_L. The pre-event increasing branch is consistent with stability, but that necessary detail is absent from the stated criterion.

**(c) Presentation problems**

- The 50 µm bending-only midpoint comparison does not state what is held fixed. It is correct for a parabola with the same rear opening, δ(L) = 200 µm. At the same M and L as the event, the traction-free midpoint is 73.37255 µm, as the submitted code itself prints. Specifying the comparison would remove the ambiguity.

**Requested numbers within tolerance:** Yes—all three meet ±0.5%, including the coarser bold values; their relative errors are 0.0121281% for M, 0.000947870% for L, and 0.0494465% for δ(L/2).
