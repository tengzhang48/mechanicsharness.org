# Grading process

- 03:35:16 UTC — Read TASK.md, the problem, rubric, and reference; began reviewing the submission. All work stays in this folder, with no web access.
- 03:35:30 UTC — Copied the submission scripts and reference solver into scratch/ for reruns.
- 03:35:47 UTC — Both reruns finished successfully: submission PID 812911 and reference PID 812912 (recorded in scratch/*.pid). All seven solver outputs agree with the key; the reference's assertions pass. Submission circle fluxes are 0.120 N/mm and the independent ODE face derivative is -2.835e-17. No process was stopped or killed.
- 03:36:23 UTC — Completed the per-number tolerance audit (scratch/audit.json; audit PID 814385, exit 0). Largest relative error in the main displayed answers is 0.008819%, below 2%. Recomputed dimensional circle fluxes to full precision. Confirmed the auxiliary debug formula omits A in w_,x and yields 16.4130064 N/mm; solve.py is correct. Assessed every rubric subitem and prepared a 100/100 grade with scientific, unsupported-claim, and presentation notes.
- 03:37:02 UTC — Wrote GRADE.md and checked its required sections, all 19 rubric rows, and the 100/100 total. All seven tolerance checks pass. Completed within the 25-minute limit; used no web search or fetch and killed no processes.
