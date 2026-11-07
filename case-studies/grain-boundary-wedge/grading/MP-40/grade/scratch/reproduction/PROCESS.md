21:02:47  START: read TASK.md; problem = GB stress relaxation by diffusion, thin film on substrate.
21:58:52  Diagnosed: earlier lambda1=11.8 came from buggy null_space on ill-scaled constraint rows; built well-conditioned basis {phi0}U{x(1-x)^3 P_j}; convergence study running (explore31).
22:10:06  KEY: same verified operator (Nval==quadrature) gives lambda1=24.8 (x(1-x)^j basis) vs ~11.8 (Legendre basis) -> polynomial collocation unreliable; eigenfunctions likely contain log boundary layers.
22:28:51 Ran explore38/out38d.txt: printed 8 eigenvalues per N; N>=24 spectra collapse (4.06, 0.0034) => Vn ill-conditioned, not a precision issue.
22:28:51 explore39.py retried N=16..36 at dps=130: identical collapse at N>=24 => basis conditioning is the limit.
22:28:51 explore40.py: N=20 eigenvector mode shape f(0.5)/f(1)=0.800103, lambda1=12.147477.
22:28:51 solve.py written (standalone, self-contained); output: lambda1=11.871 (N=16/20 Richardson), lambda2=257.493, mode1_shape_half=0.8001, K_inf_norm=1.1215, relaxation_time_t0=0.0842.
22:28:51 Wrote results.json and ANSWER.md (model, method, checks, uncertainties).
