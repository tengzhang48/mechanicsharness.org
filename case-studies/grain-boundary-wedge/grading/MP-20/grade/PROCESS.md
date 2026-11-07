# Grading process

- 19:23:37 UTC — Started review; read TASK.md, TASK_PROBLEM.md, KEY.md, and AUTO_CHECKS.json. Will work only in this folder and use no web access.
- 19:24:37 UTC — Copied submission into scratch/submission; ran its exact reproduction command (PID recorded in scratch/solve.pid). Exit 0; all five scalar claims and results.json reproduced exactly. Output saved in scratch/solve.stdout.
- 19:26:01 UTC — Inspected all delivered source files and ran scratch/audit.py. Found extra 2/pi elasticity normalization, omitted interface curvature BC, nonreal/unstable discarded spectrum, and equilibrium residual about 0.546. Auxiliary methods and claimed checks require separate assessment from the final reproduction.
- 19:27:54 UTC — Reran t_expm.py on the copy: exit 0, but zero fitted decay slope on most grids; no confirmation of the reported rate. Outputs and the process PID are saved in scratch/.
- 19:29:08 UTC — Reran t_time.py on the copy: exit 0, enormous extrapolated traction values and fitted lambda1 = -0.5461. Confirmed auxiliary symbolic/diagonal errors with local calculations.
- 19:35:33 UTC — Wrote GRADE.md. A diagnostic wrapper confirmed that t_time.py's integrator failed at tau = 0.0012612703451757136 (step size below machine spacing), while the script evaluated dense output through tau = 5. Recorded status in scratch/t_time_status.stdout.
- 19:37:34 UTC — Completed review and validation: all nine grading sections present; final JSON parses and is the last block; error/check counts match the report. Answer profile: (a) partly, (b) partly, (c) wrong; supporting state wrong; claims partly supported. All started processes completed. Worked only in this folder, retained original submission files, and used no web access.
