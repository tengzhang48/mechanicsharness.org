Boundary checks only; no overall grade is assigned.

- **B1 — existence: present_right; before.** The solver selects the real crack exponent p = 1/2 using J before the numerical work and explicitly confirms finite-energy admissibility in the answer. Evidence: “energy admissibility p < 1 satisfied” (submission/ANSWER.md:43-44); “p = 1/2 (J fixes p)” (submission/PROCESS.md:3).

- **B2 — limiting cases: present_right; before.** Reducing the supplied tensor law to uniaxial tension recovers the known bilinear law with post-yield slope T. This is a correct closed-form special-case check, already recorded in the initial plan. Evidence: “check: uniaxially ε = σ/T + σ_y(1/E − 1/T), a bilinear curve with tangent modulus T” (submission/ANSWER.md:38); “derived: uniaxial eps = sigma/T + sigma_y(1/E-1/T)” (submission/PROCESS.md:3).

- **B3 — applicability: present_right; after.** The answer states the proportional-loading condition for its constitutive treatment and identifies physical effects excluded by the membrane idealization. This scope discussion appears in the answer drafted after the numerical confirmation. Evidence: “deformation theory is used as given (proportional loading is assumed), not flow theory.” (submission/ANSWER.md:60-61); “real 3D fronts, thickness, and possible face contact are outside the 2D idealization” (submission/ANSWER.md:106-107).

- **B4 — degenerate cases: absent; unclear.** No check assesses whether the calibrated parameters approach a change in the structure or admissibility of the leading singular field. Evidence: none.

Notes: B2 counts the closed-form uniaxial check; the initial exponent and uniaxial arguments occur in submission/PROCESS.md:3 before the numerical confirmation at line 4, while the applicability discussion appears in the answer drafted at line 5. The crack-face yield transition and ν* → 1/2 observation do not assess a structural degeneration of the leading singular field or proximity to one, so B4 is absent.
