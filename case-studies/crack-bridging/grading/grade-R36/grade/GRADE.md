| rubric item | points awarded / max | reason |
|---|---:|---|
| Elastic strip bending coupled to closing face forces | 7 / 7 | Uses Euler–Bernoulli bending with `D̃ δ'''' = −b t(δ)` in the retained zone and a traction-free starter cut. |
| Individual arm bending stiffnesses | 5 / 5 | Correct total-width stiffnesses: `D₁ = 0.0495 N·m²`, `D₂ = 0.2885625 N·m²`. These equal the reference per-width stiffnesses multiplied by `b`. |
| Coupling of opening compliances | 5 / 5 | Correctly uses `D̃ = D₁D₂/(D₁+D₂) = 0.04225207986688851 N·m²`, equivalent to adding the two opening compliances. |
| Local softening law and triangular connector work | 7 / 7 | Applies the specified linear softening law to the full relative gap. `bridging_work` integrates it correctly; complete connector work is `σb δc/2 = 20 J/m²`. The computed first-event profile lies within the law's active range. |
| Front displacement and slope matching | 4 / 4 | Explicitly enforces `δ(0) = δ'(0) = 0` in the beam solution. |
| Nonzero intrinsic front work | 4 / 4 | Enforces `G_far − bridging_work(δ(L)) = Γ₀ = 1.4 J/m²`. Together with the beam first integral and front matching, this is equivalent to the reference curvature condition. The executed solution has `δ''(0) = 1.287137590621529 m⁻¹`, giving front work `1.399999999999984 J/m²`; it does not impose zero tip moment. |
| Pure-couple rear conditions and width convention | 8 / 8 | Enforces zero shear and `D̃ δ'' = M` at the free ends, with continuity through the unloaded starter cut. This yields the same conditions at the rear of the retained zone. `M` is the total couple on each strip, and the energy formula correctly includes `b`. |
| Numerical output: M | 15 / 15 | Reports `0.2126258 N·m`; relative error against the unrounded key is `0.00000593685%`. Main-code rerun gives `0.21262578737672133 N·m`. |
| Numerical output: L | 15 / 15 | Reports `10.80010 mm`; relative error is `0.00002195254%`. Main-code rerun gives `10.800102370897035 mm`. |
| Numerical output: δ(L/2) | 10 / 10 | Reports the full relative gap as `38.21890 µm`; relative error is `0.00000542639%`. Main-code rerun gives `38.218897926092424 µm`. |
| Verification | 12 / 12 | Executed `check2.py` independently constructs the reduced complex-exponential solution and reproduces all three outputs. `verify.py` also reproduces the event profile using numerical shooting and prints the correct `21.4 − 20 = 1.4 J/m²` balance. Its separate erroneous length search and repeated identical runs are not counted as successful checks. The valid checks satisfy the rubric's full-credit criterion. |
| Prescribed elastic regime and first complete failure | 4 / 4 | Uses the prescribed elastic beam model and identifies first failure as `δ(L) = δc` at the rear. The event profile has nonnegative slope, maximum opening `200 µm`, and rear relative rotation `0.04338 rad`. |
| Dependence on both arms and geometry beyond curve area | 4 / 4 | Uses both arm stiffnesses and solves the spatial beam problem for length and midpoint opening. Explicitly distinguishes the energy determination of M from the beam determination of L and the profile. |

**Total: 100 / 100.**

The reference program was rerun and reproduced all three stored targets exactly; no reference correction or alternate score is warranted. The issues below concern additional checks and continuation beyond the requested event. The rubric's required model, outputs, and valid verification are fully demonstrated, so these issues do not create an additional deduction.

**(a) Scientific errors**

- The separate shooting length search at the end of `verify.py` returns **48.705653175900125 mm**, contradicting the correct first-event length of **10.8001023708973 mm**. Its bracket spans nonphysical branches, and it does not validate the first admissible root or check `fsolve` convergence. The run emits a convergence warning. This is a failed length verification, although its earlier shooting evaluation at the correct event reproduces the requested profile. See `scratch/verify.out` and `scratch/verify.err`.
- The extended R-curve calculations continue to apply `σb(1 − δ/δc)` throughout `0 < x < L` after gaps exceed `δc`. This produces negative traction where connectors should be permanently inactive. The continuation is therefore not a valid post-failure mechanics calculation. Its output also contradicts the claimed globally monotonic continuation: `trace_rcurve.py` gives rear gaps of `725.2359 µm` at `11.9612 mm` and `216.4137 µm` at `13.6782 mm`; `final_state.py` switches to a lower-moment branch after `11.10 mm`. The actual first-event solution is unaffected.

**(b) Unsupported claims**

- `verify.py` advertises a comparison between clamped and free-moment front conditions in its module description, but implements no such comparison.
- No submitted sweep records the claimed starter-cut sensitivity check over `5–100 mm`. The grader independently confirmed the stated invariance at `5`, `20`, and `100 mm`; this is missing evidence of the claimed executed check, not an incorrect numerical or mechanical claim. No verification credit depends on that sweep.

**(c) Presentation problems**

- `verify.py` labels its numerical method “finite difference”/“FD,” but it actually uses shooting with `solve_ivp`. Its `n` argument is unused, and the loop over `400`, `800`, and `1600` performs the same calculation three times; this supplies no convergence evidence.
- The answer's description of the checks does not disclose the contradictory length printed by the separate shooting root search.

**Tolerance:** All three requested numbers are within ±0.5% of the unrounded reference; even the headline rounded values pass (errors: M `0.0121281%`, L `0.000947870%`, midpoint gap `0.00288358%`).
