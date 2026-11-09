# Grading process

- 17:07:38 UTC — Started grading; read the task, problem, rubric, reference results, and submission answer and primary code. All subsequent work is confined to this folder; no network access is used.
- 17:07:50 UTC — Finished inspecting all submitted scripts and the reference implementation; prepared to copy submitted Python files to scratch/ for execution.
- 17:08:05 UTC — Reference assertions and the submitted final_check.py completed successfully. The main solution matches the unrounded targets; the independent dislocation-density integration reproduces the root opening. Recorded the spectral comparison and identified the incorrect claim that the tip intensity vanishes at initiation.
- 17:08:41 UTC — Audited the spectral solution's actual residuals: relative tip residual −1.03×10⁻¹⁴ and relative opening residual 8.33×10⁻⁹. Its stress/length disagreements are 0.0227408% / 0.167729%. Confirmed the reference rerun exactly reproduces its recorded targets and the event lies on the increasing-load branch.
- 17:09:41 UTC — Wrote and reviewed GRADE.md: all seven rubric criteria are satisfied, totaling 100/100. Included the scientific errors, unsupported convergence claim, presentation issues, and explicit tolerance decisions. Preserved rerun evidence in scratch/; no submitted or reference files were modified, and no processes were killed.
