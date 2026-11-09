# Boundary-check coding

- **B1 existence — present_right; after.** The solver reports finite-energy and crack-field admissibility checks. The stated r^(-1/2) field has locally integrable strain energy, so this check is physically appropriate. The numerical solve precedes `check_field` in the step log (lines 3–6).
  - submission/ANSWER.md:145: “finite energy, crack exponent fixed by J, crack faces open.”
  - submission/PROCESS.md:6: “ran check_field (all PASS) and direct J (0.8 in 3 form)”
- **B2 limiting cases — absent; unclear.** No independent known-limit benchmark, such as recovering the original elastic material at T = E, is shown. The uniaxial tangent identity and small-r substitution into the original law verify the constitutive reduction. Evidence: none qualifying.
- **B3 applicability — present_right; after.** The answer states the high-stress, small-strain, and outer-matching screens, distinguishes them from finite-r error estimates, and identifies missing fracture/process-zone physics. This is consistent with the reference. The discussion is documented in the final answer, written after the computations (submission/PROCESS.md:11).
  - submission/ANSWER.md:157: “still small-strain (`max|principal strain| ≤ 0.05`)”
  - submission/ANSWER.md:164: “a process-zone/damage model — explicitly excluded”
- **B4 degenerate cases — absent; unclear.** There is no examination of coinciding exponents, a singularity threshold, or another parameter-driven change of structure. Selecting the high-stress branch and computing the screening interval do not examine such a degeneracy. Evidence: none qualifying.

Notes: The original-law small-r convergence check and uniaxial tangent identity verify the constitutive reduction, rather than benchmark the crack-field method against an independent known limiting case such as T = E. Timing follows the step log: the field checks follow the initial numerical solve, and the applicability discussion is documented in the final answer after the numerical work.
