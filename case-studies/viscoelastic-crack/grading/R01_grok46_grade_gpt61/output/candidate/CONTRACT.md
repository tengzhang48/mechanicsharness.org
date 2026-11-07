# Submission and accuracy requirements

Submit `ANSWER.md`, `solve.py`, `results.json` and any supporting code or verification data. The entry point is

```text
python solve.py case.json results.json
```

Use Python with NumPy and SciPy. Symbolic tools already present in the supplied environment are allowed. The organizer will identify the actual environment before a run. Do not require a commercial solver or package installation.

## Numerical output

The JSON object must have `critical_speed` (the supremum in units of V; use `null` if unbounded) and `cases`, with one entry per input speed. Each case must include:

- `V`, `admissible` (boolean), and, for an admissible state, `epsilon`, `a`, `R`;
- `bulk_dissipation_ratio`, computed by the spatial integration requested in the problem;
- `fields`, with one row for every pair of `sample_X` and `sample_y_over_h`. Each row has `X`, `y_over_h`, `T1`, `T2`, `E1`, `E2`.

Use finite JSON numbers. Missing or unresolved values must be `null` and explained in `ANSWER.md`; do not replace them with guessed numbers. The explanation may include additional formulas and limits that do not fit the JSON format.

## Accuracy and evidence

The accuracy requirements are known before the run:

- Each positive scalar `epsilon`, `a`, `R` and any finite `critical_speed` must be within 0.5% of the reference.
- Each normalized stress and strain component must have absolute error at most 0.005.
- `bulk_dissipation_ratio` must have absolute error at most \(0.01\max(1,\Gamma_D/\Gamma_0)\). Also require \(|\Gamma_D/\Gamma_0+1-R|/\max(1,R)\leq0.01\).
- The derivation and checks must support the reported values. Include the actual residuals or discrepancies, their normalization, and convergence or truncation evidence. A matching number alone does not establish a derivation or an independent energy calculation.

The reported outcome lists numerical accuracy, verification evidence and delivery separately, alongside a 100-point review score. Accuracy is a required gate; a high total score cannot compensate for failing it. Verification is reviewed for what it checks, without requiring a particular numerical method.

## Closed book conditions

Only the task directory and the supplied computing environment may be used. No internet, external models, source papers, other projects, organizer files or previous solutions are available. Do not search outside the workspace. Work as one agent. Record tools and resources used, elapsed work, known limitations and any recognition of a published source. Recognizing a source from memory is not itself a violation; retrieving it is outside this test.

The proposed pilot budget is six hours. Budget changes, environment details and any clarification answers will be fixed identically for all compared runs before launch. Save a reproducible partial answer if the budget expires.
