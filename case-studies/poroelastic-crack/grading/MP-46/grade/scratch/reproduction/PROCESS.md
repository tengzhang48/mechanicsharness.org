# PROCESS

21:41:34 read TASK.md; fixed the conventions: plane-strain Biot given by the task, crack -L<=x<=L, two loadings 1 (sigma_n=-H, p=0) and 2 (sigma_n=0, p=H).
21:41:34 wrote solve.py with the Rice-Cleary crack-line reduction: elastic crack K=P*sqrt(pi L), V=pi(1-nu)P L^2/G; loading 1 -> P=+1 (time-independent K, V from undrained pi(1-nu_u) to drained pi(1-nu)), small-time coefficient 2*sqrt(pi)(nu_u-nu), p_faces_0=2(nu_u-nu)/(alpha(1-2nu)).
21:41:34 loading 2 -> skeleton effective face traction t'=alpha p n = (0,-alpha) on the upper face pulls the crack closed; equivalent opening traction P=-alpha, giving K2=-alpha sqrt(pi L), V2_inf=-pi(1-nu)alpha.
21:41:34 wrote an independent transient plane-strain FEM (fem_poro.py): linear triangles, double-noded slit, backward-Euler storage, single LU factorisation, vectorised solves.
21:41:34 ran the FEM (background, PID recorded) for an elastic calibration, loading 1 and loading 2 -> fem_out.txt.
21:41:34 check A (elastic calibration, nu=0.2): FEM V=1.65178 vs exact 2.51327, COD-K=1.429 vs 1.0 -> FEM too coarse/confined, usable only for signs and trends.
21:41:34 check B (loading 1): V monotone 0.936 -> 1.6498 plateau; V1_0/V1_inf=0.75 = (1-nu_u)/(1-nu); drained-modulus COD-K plateaus at the calibration value -> K time-independent, consistent with K1_0=K1_inf.
21:41:34 check C (loading 2): V<0 for all t (closes), w(0)<0, nearly time-independent K -> confirms the closing sign and K2_0=K2_inf.
21:41:34 wrote solve.py results to results.json (written by the reproduce command) and verified by running it.