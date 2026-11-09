# Boundary-check coding

- **B1 existence: present_right; timing unclear.** The finite-energy and finite, nonzero-J conditions correctly select the square-root singularity for the derived linear leading compliance. Evidence: “Finite energy + finite non-zero J for a crack ⇒ p = 1/2” (`submission/ANSWER.md:35`).
- **B2 limiting cases: absent; timing unclear.** “Williams; HRR n→1 limit” (`submission/ANSWER.md:35`) names a known limit but does not show the solver checking its result through the same method in that limit. No such test appears in either submission file.
- **B3 applicability: present_right; timing unclear.** The answer states the loading scope and identifies omitted finite-strain/blunting physics; it also explains the leading-field screens and the missing initiation criterion. Evidence: “Loading is quasistatic, monotonic and proportional at the material points” (`submission/ANSWER.md:18–19`); “large-strain/blunting physics near the tip” (`submission/ANSWER.md:91`).
- **B4 degenerate cases: absent; timing unclear.** Neither file examines coinciding exponents or material-parameter thresholds where the singular field changes structure. There is no supporting quote for this absent check.

The step log records an analytic reduction followed by a bundled final write-up (`submission/PROCESS.md:2–3`), without locating the existence or applicability checks relative to numerical work. Naming the HRR n→1 limit does not demonstrate a limiting-case test, and the radius screens do not analyze parameter degeneracies.
