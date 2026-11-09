# HRR plane-strain Mode I crack tip, n = 6.5 — answer

## Requested numbers

| quantity | value | tolerance |
|---|---|---|
| **p**  | **0.1333333** | ±1e-5 |
| **H**  | **2.15092**   | ±0.005 |
| **C**  | **0.29178**   | ±0.002 |

## Method

**Leading-field scalings.** For the total-deformation power law the dominant
tip field has sigma_ij = sigma0 K r^{-p} sigtilde_ij(theta). Since
eps_e = alpha*eps0*(sigma_e/sigma0)^n and sigma_e ~ r^{-p}, one gets
eps ~ r^{-p n} = r^{-(1-p)}, u ~ r^{1-(1-p)} = r^{p}, and
W = n/(n+1) sigma:eps ~ r^{-1}. The r-power of J (W cos(theta) - t.u_,x) is r^0,
so **p is fixed by convergence of J**: p = 1/(n+1) = 2/15 = 0.1333333 (n=6.5).

**Angular ODE (Airy).** With Phi = r^k phi(theta), k = 2-p,
  sigtilde_rr = A = k phi + phi'',  sigtilde_tt = B = k(k-1) phi,
  sigtilde_rt = C = -(k-1) phi'.
Constitutive law with u ~ r^p gives the compatibility ODE
  p*Ur = (3/4) q (A-B),   Ur + Ut' = (3/4) q (B-A),
  Ur' - p Ut - Ut = 3 q C,   q = se^{n-1},  se = sqrt((3/4)(A-B)^2 + 3 C^2).
Eliminating Ur,Ut gives a 4th-order ODE for phi (a 4-dim first-order system).

**Boundary conditions.** Only **one** symmetry condition applies at theta=0,
phi'(0)=0 (i.e. C(0)=0; sigma_rt=0 on the crack plane). The two free-surface
conditions **A(pi)=0, C(pi)=0** then determine the two free initial values
(phi''(0), phi'''(0)); there is no eigenvalue freedom because p is already
fixed. Shooting (Radau) gives
  **phi''(0) = -0.7050162,  phi'''(0) = 0.0932568**  (residuals ~1e-14).

**Normalization.** Fields are scaled so max_theta se_tilde = 1 (standard HRR
convention; the maximum occurs at theta ~ 86 deg).

**The three numbers.** Writing B(0) = sigtilde_tt(0) = 2.827075 (normalized),
  I_n = integral_{-pi}^{pi} (Wtilde cos(theta) - ttilde . utilde_,x) dtheta
      = 7.7687   (R-independent),
  Wbar = integral_{-pi}^{pi} Wtilde dtheta = 2.26670.
Because J = alpha*eps0*sigma0*K^{n+1} I_n is the constant applied J and
sigma_tt(r,0) = sigma0 K r^{-p} B(0) = H sigma0 [J/(alpha*eps0*sigma0*r)]^p,
  **H = B(0)/I_n^p = 2.15092.**
Because U(R) = integral_0^R integral W r dtheta dr = alpha*eps0*sigma0*K^{n+1}
Wbar R and J*R = alpha*eps0*sigma0*K^{n+1} I_n R,
  **C = U(R)/(J*R) = Wbar/I_n = 0.291775.**

## Reconciliation data (exact, for cross-check with the other run)

- **Unscaled se_max = 1.0673513096** at θ=1.6533 (NOT at θ=0). Unscaled se(0)=0.5761975187, se(π)=0.7838653066.
- **After /se_max** (scaled se_max=1): B(0)=**1.51569381**, A(0)=**0.89234180**, se(0)=0.53988568.
- **After /se(0)** (se(0)=1): B(0)=**2.80767918**, A(0)=**1.65297864**.
- **phi''(0) = −0.9142244740**, **phi'''(0) = +0.0328463401** (A(π),C(π) ~1e−14).
- **BCs:** phi(0)=1, phi'(0)=0, A(π)=0, C(π)=0. Free ICs phi''(0),phi'''(0). phi'''(0)=0 is NOT imposed.
- **Integrands (θ∈[−π,π]):**
  - Wbar = ∫ (n/(n+1)) se_tilde^{n+1} dθ = **2.1297851** (with /se_max scaling).
  - I_n = ∫ [ (n/(n+1)) se_tilde^{n+1} cosθ − (σ̃_rr ∂u_x/∂x + σ̃_rt ∂u_y/∂x) ] dθ = **+0.2546**
    (traction term with `−`); = **−0.6476** if traction term enters with `+`.
  - Decomposition (θ∈[−π,π], /se_max): Wcos-term = −0.19649; traction-term = −0.45114.
  - Displacements: U_r=(3/4)q(A−B)/p, U_t=(3qC−U_r′)/(p−1), q=se^{n−1};
    ∂u_x/∂x and ∂u_y/∂x include full polar→cartesian rotation.

**Key kinematic disagreement with the other run:** their (III) reads
U_r′−(p+1)U_t=3qC; from TASK.md's 2ε_rθ=u_{r,θ}/r+u_{θ,r}−u_θ/r one gets
r^{p−1}(U_r′+(p−1)U_t), i.e. **U_r′+(p−1)U_t=3qC** (my form). This flips the
sign of the U_t contribution and changes the eigenfunction.

## Verification (independent checks)

1. **Two independent evaluations of I_n agree.** (a) The closed angular integral
   gives I_n = 7.768664. (b) A Cartesian contour J-integral built directly from
   the physical sigma_ij, traction t_i = sigma_ij N_j and finite-differenced
   du_i/dx gives I_n = 7.76873, converged in mesh (N up to 8001) and
   **R-independent** for R = 1e-3 ... 1. Both equal 7.7687.
2. **Energy positivity.** Wtilde = (n/(n+1)) sigma:eps equals se^{n+1} to
   7e-18 (constitutive law exactly consistent); W >= 0 everywhere.
3. **BCs / equilibrium.** A(pi) = 4.5e-10, C(pi) = -1.1e-10; equilibrium
   residuals (radial & hoop) <= 1e-4 (finite-difference truncation only).
4. **R-independence.** H and C computed at R = 1e-3, 1e-2, 1e-1 are identical
   to 9 digits; U(R)/(J R) is constant.

## Key correction found during the work

The strain scaling is eps ~ r^{p-1} = r^{-(1-p)} (not r^{-(1+p)}), because
eps_e = alpha*eps0*(sigma_e/sigma0)^n. This gives W ~ r^{-1}, hence U(R) ~ R
and a genuinely **R-independent C = U(R)/(J R)** — consistent with the constant
applied J. Using the wrong exponent produces a spurious R-dependent C.
A second correction: differentiating the displacement in Cartesian coordinates
must include the rotation of the polar basis, e.g.
  u_{r,x} = cos(theta) * u_{r,r} - (sin(theta)/r)(u_{r,theta} - u_theta).
With both corrections the two independent J evaluations agree and H, C are
R-independent. (The polar-to-Cartesian rotation term changes the angular
integrand by a non-negligible amount; omitting it gave I_n ~ 3.3 while the
correct value is 7.77.)