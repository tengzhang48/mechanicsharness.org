# Provided knowledge

- `TOPIC.md`: crack-tip and notch-tip asymptotic fields: governing equations and where each
  assumption enters, how the stress exponent is fixed, J and its conventions, and what changes
  outside the scope.
- `checks/crack_tip_checks.py`: checks you can run on **your own** solution, whatever method you
  used (`check_field(...)`; see the docstring). Each reports when it applies and its tolerance.
  Passing shows your field solves the stated model; whether the model fits your problem is a
  separate judgement. `checks/mutant_detection.txt` lists which checks caught which planted
  errors when the checks were validated (the validation tests are not included here).
- `entries/hutchinson_1968/`: a knowledge entry for one paper (its problem, results, conventions,
  misprints, scope) and `reference_data.json` with that paper's tabulated values.
