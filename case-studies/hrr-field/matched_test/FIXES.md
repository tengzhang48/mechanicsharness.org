# Fixes after the audit of 7 November 2026

Audit: `audit_2026-11-07/AUDIT.md` (verdicts V1 READY, V2 READY AFTER FIXES; both keys right). The statements in
`task/V1` and `task/V2`, the keys and the rubrics are kept exactly as the agents and graders had them (the staged copies
of 25 October, unchanged); a numeric checker was added. A1 of the same round is not in this folder (see `README.md`).

## Statements: kept, flaws named

- V2, lines 23 and 30: invented provenance. J is "supplied" by "an external calculation" and C_zz is "reported" by "a
  computation with this constitutive law". Neither value is in the paper; they are the data of an inverse variant, and
  the mathematics is unaffected (an admissible T exists and is unique, below).
- V2, line 44: "evaluated using the leading field" can also be read as putting the leading stresses into the complete
  law, with its σ_y term. That reading gives a lower endpoint of 0.00927 mm instead of 0.0100 mm, outside the 2% band.
  No graded answer used it; M22 considered it and set it aside (its ANSWER.md, lines 78–80). The key and rubric use the
  strains of the leading field.
- V2, line 37: the stipulation that "matching to the surrounding body" makes corrections negligible for r ≤ 0.200 mm is
  invented and no answer needs it; the complete law's own σ_y term changes ε_zz by 3.9 % to 17 % over the screened radii
  (`audit_2026-11-07/v2_sigma_y_correction.out`). Line 44 already says the screens are not error estimates.
- Cues (minor): V1 line 32, "whether these results alone determine", leans toward "no"; V2 line 21 ("no unspecified
  hardening exponent") and lines 44 and 46 point at the two decisive steps (the effective Poisson ratio and ε_zz as the
  largest principal strain).
- Relation to the paper: V1 is a variant (the antiplane crack is cited on p. 19, not solved; the key is our
  computation). V2 is adapted from Sec. 6 (Eq. 35, p. 24; Eqs. 39–42 and the β of p. 25). Its key keeps the exact
  effective Poisson ratio and so departs on purpose from the incompressible approximation of Eq. 48 (p. 27); inverting
  that approximation gives T 5.7 % high, outside the 1 % band.
- Neither statement asks for machine-readable numbers; the answers are free text.

A corrected V2 (the audit's fixes 1–3) is drafted in `audit_2026-11-07/V2_TASK_PROBLEM.md.proposed`. Under the rule for
this repository it would be a new version with new runs, never a change to this one.

## Keys

Unchanged. The audit recomputed every value with its own code, by routes that differ from the builder's at the step at
risk (V1: a stress function with compatibility governing and the exponent found by shooting, J on square and circular
contours; V2: the leading compliance as a numerical limit of the complete law, J on square contours without assuming
J = K²/T). Agreement with the key: 2.4 × 10⁻¹² or better (V1), 3.6 × 10⁻¹¹ or better (V2). V1 has one singular eigenvalue
in (0, 1) in the scanned range (`scan_V1_branches.out`); V2's inverse has exactly one root on 0 < T < E, since |C_zz|(T)
falls strictly from infinity to 8.09 × 10⁻⁴ mm^½ at T = E.

The builder's two routes share a step (V1: the amplitude function; V2: K = √(JT)), so they were only partly independent.
A check of V1 by Opus on 7 October is recorded as one line of a process log, without a script or output; it is not
counted.

## Rubrics

Unchanged. The consequential half-credit rule (V1/RUBRIC.md, lines 19 and 33; the same in V2) gave M20 76 of 100 with
none of its 7 numbers within tolerance, and M08 lost 2 points for an endpoint 134 % off. As for T1 and T2, report whether
all numbers are within tolerance beside the total.

## Checker (`key/check_V1_V2.py`, new)

The audit's checker, with two changes for new answers: a value that is not a number, or a block that is not an object,
fails instead of stopping the check; and an input without a V1 or V2 block exits with status 1 ("nothing was
checked") instead of exiting silently with status 0 (Codex's review, 8 Nov). It applies the stated tolerances (V1: all seven numbers 2 %; V2: T 1 %, the other
six 2 %; relative, signed) to `{"V1": {slip_0.02, slip_0.08, tau_rz, tau_tz, r_0, r_90, r_180}, "V2": {T, C_rr, C_tt,
C_g, C_y, r_min, r_max}}`. The key passes 14 of 14. Every one of the audit's 27 mutants fails the rule "all seven within
tolerance", but some single items pass for wrong models (ν_t = ½: C_γrθ; the plane-strain J–K relation: C_rr, C_γrθ and
C_y), so that is the pass rule. The justification, verification and scope rows of the rubrics have no numeric check.

## Saved runs checked

`audit_2026-11-07/check_V1_V2_audit.out` (numbers transcribed from each graded ANSWER.md, with their lines, in
`answers_extracted.json`): the checker agrees with every grader's tolerance verdict, 24 of 24. 22 answers have all seven
numbers within tolerance; M20 (V1, condition C, DeepSeek) has none, from a parity error in its J integral
(`runs/P015-V1-C-r2-deepseek/workspace/solve.py`, lines 44–46, reproduced by the audit's mutant); M08 (V2, B0, GLM) has
six, its lower screening endpoint 0.02336 mm instead of 0.0100 mm (the trace term left out of the in-plane strains). The
three complete Kimi K3 V1 answers, which were not graded, have seven of seven. No grade was changed.

Not examined by the audit: the rubrics' judgement rows, the knowledge files, the session logs for rule compliance (only
file times and byte identity of the graded copies were checked), and A1.
