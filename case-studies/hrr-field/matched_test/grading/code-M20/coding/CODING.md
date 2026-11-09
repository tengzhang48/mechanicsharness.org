# Boundary-check coding

- **B1 existence — present_right; before.** The finite, radius-independent energy-flux condition selects the admissible exponent beta = 1/(n+1) = 2/9. This agrees with the reference. The step log records the derivation at line 2, before numerical integration at line 3 and load normalization and requested outputs at lines 4–6.
  - `submission/ANSWER.md:40-41`: “Path-independence/finiteness of $J$ for a single power law forces the HRR exponent $\beta=1/(n+1)$”
  - `submission/PROCESS.md:2`: “derived HRR mode-III power-law model: w=A r^beta f(theta), beta=1/(n+1)=0.2222”
- **B2 limiting cases — present_right; after.** The answer states the correct linear n = 1 angular form for the same model, and the log compares the computed face angular value with its exact expression. The logged comparison follows the requested numerical outputs at line 6; the n = 1 statement appears only in final verification.
  - `submission/PROCESS.md:7`: “cross-checked f(pi) against closed form (n+1)n^{-n/(n+1)}=1.69843; agree”
  - `submission/ANSWER.md:57`: “$n=1$ the model reduces to the classical linear mode-III field ($f\propto\sin(\theta/2)$).”
- **B3 applicability — present_right; after.** The final answer identifies the assumed local-field/outer-geometry scope and the missing fracture criterion. Those limitations are physically correct for the prescribed stationary calculation. These scope statements are documented only in the final answer, written after the numerical work.
  - `submission/ANSWER.md:62`: “geometry outside $r=0.30$ mm is negligible; no far-field/loaded-boundary condition was imposed.”
  - `submission/ANSWER.md:68-69`: “include no fracture/toughness criterion (e.g. a critical $J_c$, crack-tip slip, or critical strain)”
- **B4 degenerate cases — absent; unclear.** Neither file examines coinciding exponents, a singularity threshold, or another parameter-driven change of solution structure. There is no supporting quotation and no check to time.

Notes: These codes assess only the boundary checks, not the requested numerical answers; timing follows the step log, with final-answer-only checks counted as after. The correct n=1 limiting form is stated without a separately logged numerical n=1 run, and no parameter threshold or change of solution structure is examined.
