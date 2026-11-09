**B1: present_right; timing: unclear.** Shooting confirms a real, admissible eigenvalue for the traction-free angular field, consistent with the reference exponent. Evidence: “s* = 0.22222222222220 in (0.03, 0.6), matching the known exponent 1/(n+1)” (submission/ANSWER.md:38); “Ran solve.py: s*=1/(n+1) to 2e-14, all J paths = 0.120, n=1 limit exact; got 7 numbers.” (submission/PROCESS.md:3).

**B2: present_right; timing: unclear.** The same numerical method recovers the correct linear-elastic exponent, angular normalization, and flux integral. Evidence: “the same pipeline returns s = 0.5, I_ang = π, f(π) = 2 exactly” (submission/ANSWER.md:66); “Ran solve.py: s*=1/(n+1) to 2e-14, all J paths = 0.120, n=1 limit exact; got 7 numbers.” (submission/PROCESS.md:3).

**B3: present_right; timing: after.** The stated proportional-loading condition and the missing fracture criterion correctly delimit the stationary constitutive model. Evidence: “proportional loading ⇒ deformation theory is exact here” (submission/ANSWER.md:54); “initiation requires an additional fracture criterion” (submission/ANSWER.md:88).

**B4: absent; timing: unclear.** Neither file checks coinciding exponents, a singularity threshold, or another parameter-dependent degeneracy. Evidence: none.

The eigenvalue and n = 1 checks share the main numerical-run entry, so their timing is unclear; applicability is documented in the later final answer (submission/PROCESS.md:5). The n = 1 check does not examine coalescing exponents or a singularity threshold, so it does not establish B4.
