# Scientific review profile for the analytical film crack task

This rubric applies to the stated problem and the independently established
reference. Review a neutral submission without model labels or desired
verdicts. Preserve original mathematics; identify substantive omissions or
incorrect assumptions even when numerical values agree.

For each item return **right**, **partly**, **wrong** or **missing**, with a
short reason, evidence location and the deciding check. A correct partial
calculation need not make its entire item right. Distinguish unsupported
claims, algebra errors, implementation errors, missing delivery and limits
of the physical model. Do not convert command success into scientific
correctness.

1. **Wake formulation and fields:** governing reduction, decay length,
   displacement/stress fields, total crack-face traction, and states ahead
   of and behind each geometry.
2. **Energy release:** complete stated energy, correct two geometries,
   translating-tip argument, period width and moving-tip/thickness
   normalization, and both energy-release expressions.
3. **Fracture modes and intensities:** justified symmetry and use of the
   supplied plane-stress relation, with the required normalized intensities.
4. **Propagation thresholds and limits:** critical reference stress,
   small/large-spacing limits with fixed quantities stated, and a mechanical
   interpretation consistent with the expressions.
5. **Requested numerical values:** sixteen samples with the requested
   normalization and locations. Report agreement count separately from this
   scientific assessment; state any corrected mapping or unresolved output.
   (26 Oct: the agreement count is `checks.py`, 0.5 % relative; the statement
   gives no accuracy and leaves the method open.)
6. **Verification and reproducibility:** actual executed checks, their
   outputs, tolerances, failure sensitivity and shared assumptions; agreement
   between code and the claimed answer; reproducibility of the declared
   command. Perform independent deciding-quantity checks rather than merely
   rerunning the same formula.
7. **Physical scope and uncertainty:** applicability of the foundation
   reduction, what calibration establishes, and which physical conclusions
   require more information or a richer model. Avoid treating every formal
   small-spacing limit as automatically validated or invalidated for a real
   substrate.

Provide GRADE.md and grade.json with the seven judgments, independent-check
details and actual output files, unresolved checks, errors, uncertainty and
source recognition found in the submission. The reviewer writes the
scientific judgments; the coordinator does not supply expected grades.

Administrative checks may record file presence, JSON parsing and raw numeric
agreement. Those do not replace this profile and do not constitute physical
validation. Preserve uncompleted review obligations explicitly if the review
budget ends.
