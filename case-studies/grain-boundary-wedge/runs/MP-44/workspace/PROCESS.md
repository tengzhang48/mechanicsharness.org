# Process log (MP-44)

21:03:57 Read TASK.md. Problem: stress relaxation by GB diffusion in a thin film on a substrate; plane strain half-space, GB from z=0 (free surface) to z=h (film/substrate interface). Need lambda1, lambda2, mode1_shape_half=f(0.5), K_inf_norm, relaxation_time_t0, method.
21:25:59 Ran cheb5.py (weighted-Chebyshev, family A = (Z(1-Z))^-1/2 T_k on f'''): R[f]/f constant at 83-88 across Z -> genuine eigenmode; family B gives garbage. Discretisation works.
21:26:42 Ran cheb6.py: family A eigenvalues vs N: lambda1 -> 88.33, residual dropping to ~1e-4. Slow convergence traced to quadrature near Z->0.
21:27:53 Ran cheb7.py (N=60, mixed A+B): lambda1=88.33175, resid 3.6e-7. Adding family B does not lower lambda1 (pure A is cleanest).
21:28:57 Ran solve.py (first version): lambda1=88.311; drifted with quadrature size -> corner at Z=0 needs treatment.
21:30:10 Ran sens.py: lambda1 drifts 88.376->88.258 as regular-part Gauss points 200->3200; drift O(1/n). Confirms corner singularity; need paneled quadrature.
21:31:39 Ran solve2.py (paneled) -> Rf wrong (phi range bug); abandoned.
21:32:27 Ran brute.py (mpmath PV, no identity): true R[f_1](0.5)=4.833361. cheb5's Rval=2.420171 = exactly half -> cheb5 had a uniform factor-2 in BOTH R and M, so its eigenvalues are still valid.
21:33:46 Ran solve3.py (correct factors: Cauchy= -2pi U, theta_int x2): Rf matches brute to <0.1%; lambda1=88.24907 at N=60, well converged.
21:34:55 Ran validate2.py: independent brute-force R[f]/f = 88.20-88.21 (vs 88.249) -> confirms lambda1.
21:35:53 Ran solve3b.py: lambda1 = 88.24905, 88.24904, 88.24904, 88.24905 for N=64,70,80,90 -> lambda1=88.249 (5 digits); lambda2 -> 1361.85.
21:37:37 Ran solve.py (final): lambda1=88.2491, lambda2=1361.871, mode1_shape_half=0.9705, K_inf_norm=1.1215, relaxation_time_t0=0.01133. Eigenmode residual 1.07e-07; mode coefficients decay geometrically.
21:38:37 solve.py with added independent paneled PV check (no Hilbert identity): agrees to 4.37e-11 -> operator verified. Wrote results.json, ANSWER.md.
21:40:29 Reran reproduce command python3 solve.py -> lambda1=88.2491 lambda2=1361.871 mode1_shape_half=0.9705 K_inf_norm=1.1215 relaxation_time_t0=0.01133, matches results.json. Fixed ANSWER.md paragraph order; deliverables complete.
