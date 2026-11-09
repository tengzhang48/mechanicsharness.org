- **B1 existence — present_right; after.** Independent shooting confirms that the exponent required by finite nonzero J satisfies the traction-free angular boundary condition. Evidence: “independently shooting m to enforce g(π) = 0 recovers m = 0.777777777778 = n/(n+1) to 12 digits.” (submission/ANSWER.md:36); “First run: eigen-checks pass” (submission/PROCESS.md:2).

- **B2 limiting cases — present_right; after.** The same solver recovers the n = 1 linear angular solution and the correct linear-elastic J–K relation. Evidence: “the amplitude relation reduces to J = K²/(2μ), μ = τ₀/γ₀ = 40 000 MPa” (submission/ANSWER.md:38); “n=1 gives -cos(th/2)” (submission/PROCESS.md:2).

- **B3 applicability — present_right; after.** The answer identifies the outer limit of the leading field and the inner failure of small-strain physics. Evidence: “Higher-order (non-singular) terms from the outer geometry are neglected per the problem” (submission/ANSWER.md:43); “Small strain fails at the very tip (γ → ∞ as r → 0)” (submission/ANSWER.md:44).

- **B4 degenerate cases — absent; unclear.** Neither file explicitly examines coincident exponents, a singularity threshold, or another parameter-driven change of structure. Evidence: none.

B1 and B2 are first recorded as checks in the first numerical run (submission/PROCESS.md:2) and repeated in final verification; B3 is documented in the final answer written at submission/PROCESS.md:6, so these are coded after. The n = 1 comparison checks a known limit but does not explicitly examine a structural degeneracy or a singularity threshold.
