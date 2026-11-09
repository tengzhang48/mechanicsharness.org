# HRR plane-strain Mode I crack tip, n = 6.5 - answer (CORRECTED)

## Requested numbers

| quantity | value | tolerance |
|---|---|---|
| **p** | **0.1333333** | +/-1e-5 |
| **H** | **1.4706**    | +/-0.005 |
| **C** | **1.6978**    | +/-0.002 |

Full precision: p = 0.1333333333, H = 1.470562, C = 1.697757.

## Method

**Scaling / exponent p.** With the total-deformation law
eps_ij = (3/2) alpha eps0 (sigma_e/sigma0)^{n-1} s_ij/sigma0 and
sigma_e ~ r^{-p}, the strain scales as eps ~ r^{-p n} = r^{-(1-p)} (because
n p = n/(n+1) = 1-p), so u ~ r^{1-(1-p)} = r^{p} and
W = n/(n+1) sigma:eps ~ r^{-1}. Then every term of J (W cos(theta) and
t_i u_{i,x}) scales as r^{0}, so J is constant in r and
**p = 1/(n+1) = 2/15 = 0.13333333** exactly.

**Angular ODE (Airy).** Phi = r^{k} phi(theta), k = 2-p, so
  sigtilde_rr = A = k phi + phi'',  sigtilde_tt = B = k(k-1) phi,
  sigtilde_rt = C = -(k-1) phi',  se = sqrt((3/4)(A-B)^2 + 3 C^2).
With u_i = r^{p} U_i(theta) and q = se^{n-1}, the constitutive law + kinematics:
  (I)   p U_r            = (3/4) q (A-B)
  (II)  U_r + U_t'       = (3/4) q (B-A)
  (III) U_r' + (p-1) U_t = 3 q C
(III) uses 2 eps_rtheta = u_{r,theta}/r + u_{theta,r} - u_theta/r, i.e.
r^{p-1}(U_r' + (p-1) U_t); cf. TASK.md. Eliminating U_r, U_t and using
(A-B)'' = k p phi'' + phi''''  (since 2-k = p) gives a 4th-order ODE for phi
(a 4-dim first-order system), integrated with Radau (hrr_clean.py).

**Boundary conditions.** Mode I symmetry at theta = 0 gives only
sigma_rtheta(0) = 0, i.e. C(0)=0, i.e. phi'(0)=0. The two free-surface
conditions at theta = pi, **A(pi)=0 and C(pi)=0**, determine the two free
initial values phi''(0), phi'''(0) (p is already fixed). Shooting from
phi(0)=1, phi'(0)=0:
  **phi''(0) = -0.9142244740,  phi'''(0) = +0.0328463401**,
  giving A(pi) = -6.5e-12, C(pi) = +4.3e-12.
Two independent solvers (shooting fsolve + solve_bvp) agree.

**The three numbers (RAW normalization, phi(0)=1, NO re-scaling of se).**
  se(0)=0.576198, se_max=1.067351 (at theta=1.653), se(pi)=0.783865,
  B(0)=1.617778, A(0)=0.952442,
  I_n = int_{-pi}^{pi} (Wtilde cos(theta) - ttilde_i utilde_{i,x}) dtheta = 2.045340,
  Wbar = int_{-pi}^{pi} Wtilde dtheta = 3.472491.
With J = alpha eps0 sigma0 K^{n+1} I_n and
sigma_tt(r,0) = sigma0 K r^{-p} B(0) = H sigma0 [J/(alpha eps0 sigma0 r)]^{p}
(note [J/(alpha eps0 sigma0 r)]^{p} = [K^{n+1} r^{-1}]^{p} = K r^{-p}):
  **H = B(0) / I_n^{p} = 1.470562.**
With U(R) = int_0^R int W r dtheta dr = alpha eps0 sigma0 K^{n+1} Wbar R and
J R = alpha eps0 sigma0 K^{n+1} I_n R:
  **C = U(R)/(J R) = Wbar / I_n = 1.697757.**

## Verification (independent checks)

1. **Scale-invariance of H and C.** Rescaling phi -> lambda phi (tested
   lambda = 1, 2, 5) leaves H and C bit-identical (1.470562, 1.697757) while
   B(0) and I_n scale. These are the only scale-free combinations.
2. **ODE validated against the n=1 analytic linear-elastic (Williams) limit.**
   Substituting the exact biharmonic eigenfunctions phi = cos(theta/2),
   cos(3 theta/2) into the RHS reproduces phi'''' to ~1e-16 (validate3.py).
3. **Independent Cartesian contour J.** Building physical sigma_ij, traction
   t_i = sigma_ij N_j, du_i/dx (with polar-basis rotation:
   u_{r,x} = cos(theta) u_{r,r} - (sin(theta)/r)(u_{r,theta}-u_theta)) gives
   J = 2.045340, **identical for R = 0.02 ... 0.7** (R-independent), matching
   the angular integral (jcheck.py).
4. **Constitutive/energy consistency.** p U_r - (3/4) q(A-B) = 2.8e-17,
   U_r + U_t' - (3/4) q(B-A) = 1.1e-16, U_r' + (p-1)U_t - 3 q C = 3.3e-16.
   Wtilde = (n/(n+1)) sigma:eps = se^{n+1} identically; W >= 0.
5. **Free-surface BCs.** A(pi) = -6.5e-12, C(pi) = +4.3e-12.

## Notes on the earlier (superseded) values

Earlier drafts reported H = 2.15092, C = 0.29178 (and, in one draft,
I_n = 7.7687). Those are **invalid**: (i) an intermediate harness had a sign
error in the eps_rtheta kinematic relation, and (ii) I_n is **not invariant**
under re-scaling se (the energy term ~ se^{n+1} and the traction term ~ se^{n}
scale differently), so computing I_n from se_max-normalized fields corrupts it.
Both were fixed: the ODE is validated against the n=1 Williams limit, and
I_n, B(0), Wbar are all taken from the same raw fields. The corrected invariant
results are **p = 0.1333333, H = 1.4706, C = 1.6978.**

## Cross-check data (for reconciliation)

- phi''(0) = -0.9142244740, phi'''(0) = +0.0328463401.
- RAW: se(0)=0.576198, se_max=1.067351, se(pi)=0.783865, A(0)=0.952442, B(0)=1.617778.
- I_n = 2.045340 (angular integral == Cartesian contour J, R-independent).
- Wbar = 3.472491.  H = B(0)/I_n^p = 1.470562.  C = Wbar/I_n = 1.697757.
- Key discriminators: ODE reproduces n=1 Williams exactly; epsilon_rtheta
  sign; se-rescaling destroys I_n (must use raw fields).
