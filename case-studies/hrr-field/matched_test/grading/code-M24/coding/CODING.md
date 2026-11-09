# Boundary-check coding

- **B1 existence — present_right; timing unclear.** The claimed unique inverse and admissible T agree with the reference; the step log does not date this check.
  Evidence: “Solving (brentq; LHS monotone in T, unique root) gives **T = 4816.47 MPa**” (`submission/ANSWER.md:32`); “t = T/E = 0.068807” (`submission/ANSWER.md:7`).
- **B2 limiting cases — absent; timing unclear.** No known-limit verification, such as recovery of the elastic problem at T = E, is documented. Evidence: none.
- **B3 applicability — present_right; timing after.** The strain/matching restrictions and missing initiation criterion agree with the reference. The scope review follows the numerical results in the step log.
  Evidence: “infinitesimal strain still defensible (principal strains ≤ 5%), and leading-field validity (matching corrections negligible for r ≤ 0.2 mm)” (`submission/ANSWER.md:54`); “no failure criterion (critical J or K, CTOD, critical strain/stress, damage model) is supplied” (`submission/ANSWER.md:56`).
- **B4 degenerate cases — absent; timing unclear.** No check of proximity to a structural degeneracy is documented. Evidence: none.

The monotonicity/unique-root claim and reported 0 < T/E < 1 support B1, but its timing is not recorded; the scope review follows the main numerical work (submission/PROCESS.md:3, :5, :8). The J, kinematic, and finite-radius compliance checks do not test a known limiting case, and the principal-strain dominance comparison and nonempty annulus do not examine proximity to a structural degeneracy.
