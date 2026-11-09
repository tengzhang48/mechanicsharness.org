# HRR plane-strain Mode I angular solution, n = 6.5 — session B findings

(Because ANSWER.md is being written concurrently by another run, my independent
findings are recorded here. See also PROCESS.md "Session B".)

**Exact exponents (J-integral dimensional analysis):**
p = 1/(n+1) = 1/7.5 = 0.13333333…, k = 2 − p = 1.86666666…

**Airy formulation** Φ(r,θ)=r^k φ(θ):
- A = σ̃_rr = kφ + φ″ ; B = σ̃_tt = k(k−1)φ ; C = σ̃_rt = −(k−1)φ′
- s_rr=(A−B)/2, s_tt=(B−A)/2, s_rt=C, se² = (3/4)(A−B)² + 3C²

**Constitutive (u∼r^p, ε∼r^{p−1}), from TASK kinematics:**
- (I)  pU_r = (3/4)q(A−B)
- (II) U_r + U_t′ = (3/4)q(B−A)
- (III) U_r′ + (p−1)U_t = 3qC ,  q=se^{n−1}

(Note: 2ε_rθ = u_{r,θ}/r + u_{θ,r} − u_θ/r = r^{p−1}(U_r′+(p−1)U_t).)

## 1. Boundary-value analysis and robust solution

✅ The ODE is 4th order. The only axis condition at θ=0 is C(0)=0 ⟺ φ′(0)=0.
Imposing φ‴(0)=0 too over-determines it (no even solution exists; verified by a
determinant test). 

**Robust method:** shoot from θ=0 with φ(0)=1, φ′(0)=0 and **two free** ICs
φ″(0), φ‴(0); impose A(π)=0, C(π)=0 (2 eqns / 2 unknowns). Converges by
`fsolve` AND independently by `solve_bvp`. 

**φ″(0) = −0.9142244740 , φ‴(0) = +0.0328463401 , A(π)=1.9e−14 , C(π)=−9.2e−15.**

## 2. Angular functions, normalized se(θ=0)=1

| quantity | value (se(0)=1) | value (se_max=1) |
|---|---|---|
| A(0)=σ̃_rr(0) | **1.65298** | 0.89234 |
| B(0)=σ̃_tt(0) | **2.80768** | 1.51569 |
| C(0)=σ̃_rt(0) | 0.00000 | 0.00000 |
| se(0) | 1.00000 | 0.53989 |

Unscaled: se(0)=0.5761975187, se_max=1.0673513096 (at θ=1.653), se(π)=0.7838653.

Angular shapes (se(0)=1):
```
 th      A        B        C       se
0.000  +1.65298 +2.80768  0.00000  1.00000
0.262  +1.70379 +2.72126 +0.34795  1.06754
0.524  +1.78920 +2.47770 +0.63484  1.25083
0.785  +1.74562 +2.11555 +0.83304  1.47800
1.047  +1.49712 +1.67514 +0.95961  1.66923
1.309  +1.13382 +1.18589 +1.03507  1.79336
1.571  +0.72456 +0.67050 +1.06718  1.84900
1.833  +0.31489 +0.14973 +1.05692  1.83621
2.094  -0.04610 -0.35502 +1.00046  1.75338
2.356  -0.27952 -0.81874 +0.88597  1.60402
2.618  -0.29479 -1.20739 +0.68811  1.43008
2.880  -0.15311 -1.47406 +0.38612  1.32512
3.142   0.00000 -1.57087  0.00000  1.36041
```

## 3. Displacement angular functions
U_r = (3/4)q(A−B)/p , U_t = (3qC − U_r′)/(p−1) , q=se^{n−1}.
(II)≡(III) hold to ~1e−16 with ODE-consistent U_r″ (`diag3.py`).

## 4. Integral I_n
I_n = ∫_{−π}^{π}[ (n/(n+1)) se^{n+1} cosθ − (σ̃_rr ∂u_x/∂x + σ̃_rt ∂u_y/∂x) ]dθ
with FULL polar→cartesian rotation of ∂u_i/∂x.

- /se_max scaling: Wbar=∫(n/(n+1))se^{n+1}dθ = **2.12979**;
  Wcos-term = −0.19649; traction-term = −0.45114; **I_n = +0.2546**
  (traction with `−`); −0.6476 if `+`.
- The magnitude of I_n is HIGHLY sensitive to the displacement normalization
  convention (the other run reports I_n=7.77 — see reconciliation).

## Verification (independent checks)
1. Corrected ODE reproduces the exact n=1 Williams biharmonic eigenfunctions
   φ=cos(θ/2) and cos(3θ/2) to 1e−16 (`validate3.py`); the old hrr2.py gave O(1)
   residuals (bug: wrong displacement exponent r^{−p} instead of r^{p}).
2. Two independent solvers (shooting `fsolve`; `solve_bvp`) give identical
   φ″(0), φ‴(0), A(π)≈0, C(π)≈0.
3. ODE self-consistency: dφ‴/dθ from finite differences = RHS to ~1e−11.
4. (II)≡(III) constitutive consistency to ~1e−16.

## Open discrepancy with concurrent run — RESOLVED by the n=1 test
They use (III) U_r′−(p+1)U_t=3qC (U_t sign flipped vs TASK kinematics) and
obtain φ″(0)=−0.705, φ‴(0)=0.093 — which is NOT a root of my ODE (residual 0.44).

**Discriminating test (decisive):** at n=1 the ODE must be satisfied by the exact
biharmonic (Williams) eigenfunctions φ=cos(θ/2), cos(3θ/2).
- My (III) U_r′+(p−1)U_t=3qC: residuals ~1e−16 (exact). ✅
- Their (III) U_r′−(p+1)U_t=3qC: residuals 0.66 and 1.6. ❌
So my (III) is the correct one (it is the unique choice reducing to linear
elasticity). Kinematic origin from TASK 2ε_rθ=u_{r,θ}/r+u_{θ,r}−u_θ/r:
r^{p−1}(U_r′+pU_t−U_t)=r^{p−1}(U_r′+(p−1)U_t).

## Evenness
For every n>1, the even ansatz (φ‴(0)=0) with the single free DOF φ″(0) cannot
satisfy both A(π)=0 and C(π)=0 (C(π)≠0 at each A-root; `even_scan.py`). Hence
the HRR Mode I Airy eigenfunction genuinely has φ‴(0)≠0 — not an artifact.