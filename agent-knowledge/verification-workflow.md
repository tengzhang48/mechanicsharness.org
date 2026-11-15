# MechanicsHarness — Optional Verification Workflow

**Version:** 0.1.0 (candidate)  
**Audience:** AI solvers in the *physics-plus-verification* test condition  
**Input:** a mechanics task and, when the experimental design calls for it, the separate `fracture-101.md` physics reference.

This file is a process intervention, **not** additional mechanical formulas or an answer key. Use it only when the benchmark organizer explicitly includes it in the workspace. All task-specific instructions and output contracts take precedence.

## 1. Establish the model

1. Write down geometry, reference coordinates, signs, dimensions, loads, boundary/initial conditions and requested quantities.
2. Identify constitutive assumptions and the range of validity of the proposed simplifications.
3. Distinguish mathematical limits (e.g., near tip versus large scale, early versus late time) before using asymptotic results.

## 2. Produce a primary solution

1. Derive the governing relations that are actually needed; identify every result taken as a standard identity instead of derived.
2. Implement and run any numerical computation needed. Save code and parameters when the task requires reproducibility.
3. Record units, normalization, selected branch and approximations. Do not fit an unknown quantity by secretly using a reference answer.

## 3. Check independently

Choose checks that do **not** merely reuse the same discretized equations or fitted constants:

- Dimensional consistency and physical signs.
- A known exactly solvable limit of the same governing problem, when one exists.
- Energy, momentum, mass, or work balance, when applicable.
- Reciprocal or variational identities for conservative elastic systems, when the governing assumptions permit them.
- Alternative observable, formulation, or extraction route for quantities inferred from a numerical field.
- Grid, domain, quadrature, fitting-window or timestep convergence as appropriate.

When a check fails, locate the disagreement; do not relabel it an acceptable numerical error without evidence.

## 4. Report what was actually checked

For every claimed check, record:

- The tested statement and its domain of validity.
- The script or derivation used and its inputs.
- The observed value or residual and the acceptance rule.
- Whether the check is independent of the primary route.
- Any unresolved discrepancy, approximation, or time/tool limitation.

Distinguish **computed**, **derived**, **recalled**, **assumed**, and **not verified**. Never claim to have executed an unavailable tool or checked a result against data you were not given.

## 5. Submit to the task contract

Respect the original task's result format and scope. Keep the verification notes separate from required output data unless instructed otherwise. State uncertainty when an answer depends on an unresolved model choice.

**Boundary:** Do not search for or access the benchmark key, grader, other agent runs or hidden case materials. Do not use web or outside files when a task forbids them.
