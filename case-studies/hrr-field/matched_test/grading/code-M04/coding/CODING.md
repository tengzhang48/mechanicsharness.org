- **B1 existence: present_right; before.** The shooting check confirms the face boundary condition at the admissible exponent q = 1/(n+1), consistent with the reference field. Evidence: “shooting F(0)=0,F'(0)=1 gives F'(pi)=0 exactly at q=1/(n+1)=0.2222222” (submission/PROCESS.md:2); “F′(π) = 0 (residual −6×10⁻¹⁷) **exactly at q = 1/(n+1)**” (submission/ANSWER.md:53).

- **B2 limiting cases: present_right; after.** The same code is checked against the known linear mode-III solution at n = 1; the stated face slip agrees with that solution. Evidence: “n=1 limit equals K^2/(2mu)” (submission/PROCESS.md:4); “the closed-form Δw(0.02) = 5.5279×10⁻⁴ mm matches the code.” (submission/ANSWER.md:58).

- **B3 applicability: present_right; after.** The answer states the local-field and proportional-loading conditions and correctly identifies omitted physics. Evidence: “the **leading local field only**, which the task authorizes for r ≲ 0.30 mm, and used **deformation theory under proportional loading**” (submission/ANSWER.md:49); “Finite-strain blunting, the fracture process zone, and any elastic term are omitted” (submission/ANSWER.md:67).

- **B4 degenerate cases: absent; unclear.** Neither file examines coinciding exponents, a singularity threshold, or proximity to a structural change. No evidence quote.

Timing follows step-log order: the angular solve and amplitude calculation at submission/PROCESS.md:3 begin the main numerical work, and applicability is documented only in the final answer. The n = 1 comparison and n sweep do not analyze a degeneracy or singularity threshold.
