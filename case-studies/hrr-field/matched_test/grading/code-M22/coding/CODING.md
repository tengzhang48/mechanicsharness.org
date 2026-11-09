Timing is relative to the first numerical solve at `submission/PROCESS.md:2` (03:35:47).

**B1 existence — present_right; after.** The solver checks compatibility of the assumed mode-I field and explicitly finds a nonempty radius range satisfying the leading-field screens. Both checks follow the first solve in the step log (lines 3 and 5), and are consistent with the reference.

Evidence: “checked that the mode-I displacement field reproduces the stress-law strains” (`submission/PROCESS.md:3`); “interval is non-empty: **[0.01000, 0.06133] mm**.” (`submission/ANSWER.md:70`).

**B2 limiting cases — present_right; unclear.** The reduced tensor compliance is compared with the closed-form uniaxial bilinear law: its leading compliance is 1/T, matching the tensile tangent of that special case. This is an analytic special-case check, but the step log does not establish its timing relative to the numerical solve.

Evidence: “so the uniaxial compliance is A+B = 1/T” (`submission/ANSWER.md:23`); “This reproduces the bilinear uniaxial law ε = σ_y/E + (σ−σ_y)/T, whose tangent is T” (`submission/ANSWER.md:28`).

**B3 applicability — present_right; after.** The answer states the high-stress, small-strain, and outer-matching screens, distinguishes them from finite-r error bounds, and identifies omitted crack-growth, damage, and inertia physics. It correctly explains that initiation requires additional fracture data. The log records the screens and scope discussion after the first solve (lines 5–6).

Evidence: “(ii) all principal strains ≤ 5% (infinitesimal kinematics are defensible)” (`submission/ANSWER.md:104`); “(crack growth, damage and inertia are explicitly excluded).” (`submission/ANSWER.md:111`).

**B4 degenerate cases — absent; unclear.** Neither document investigates proximity to a parameter threshold, coinciding exponents, or the appearance/disappearance of a singularity. Evidence: none.

Notes: B1 credits the compatibility check and explicit nonempty leading-field screening range; B2 credits the stated comparison with the closed-form uniaxial tangent law, whose timing is not logged. No T = E recovery test or parameter-driven change of field structure is recorded, and numerical consistency checks alone are not treated as such tests.
