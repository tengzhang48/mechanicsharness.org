# Boundary-check coding

**B1 existence — present_right; after.** The finite, radius-independent J condition selects s = 1/(n+1) = 2/9, and the numerical eigenvalue and face residual confirm an admissible field. This agrees with the reference: W scales as 1/r, giving finite local stored work. The explicit verification follows the main results in the step log.

- submission/ANSWER.md:30: “shooting eigenvalue reproduces s* = 0.2222222222 with face residual Φ(π) = 1.5×10⁻¹⁶,”
- submission/PROCESS.md:5: “Verified: Phi(pi)=1.5e-16, J_direct=0.120000000 at r=0.02/0.08/0.20”

**B2 limiting cases — present_right; after.** The same pipeline is checked at n = 1 against the linear mode-III solution. Its modulus, slip formula, stresses, and strain-contour expression are consistent with J = K_III²/(2μ). The step log reports the comparison after the main results.

- submission/ANSWER.md:50: “setting n = 1 (μ = τ₀/γ₀ = 40 000 MPa) the whole pipeline reproduces the exact LEFM K-field”
- submission/PROCESS.md:5: “elastic limit matches LEFM to ~1e-16”

**B3 applicability — present_right; after.** The answer identifies the leading-field approximation and its outer-geometry limitation, and explains that stationary constitutive fields do not supply a fracture-growth criterion. These are valid scope limitations. The step log places the answer with caveats after the results and numerical verifications (submission/PROCESS.md:6).

- submission/ANSWER.md:56: “truncation to the leading (separable) field, assumed negligible for r ≪ 0.30 mm”
- submission/ANSWER.md:58: “without an additional growth criterion (e.g., a critical J, CTOD, or cohesive strength)”

**B4 degenerate cases — absent; unclear.** Neither file examines coinciding exponents or a parameter threshold where the singularity appears, disappears, or otherwise changes structure. There is no supporting quote and no check to time.

**Notes:** Timing follows the step log: results appear at line 4, explicit numerical verifications at line 5, and the answer with caveats at line 6; the earlier solver-writing entry does not establish that a check had been performed. The n = 1 comparison is a regular limiting-case check, with no examination of coinciding exponents or a parameter threshold where the field changes structure; no overall grade is assigned.
