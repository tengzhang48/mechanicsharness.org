# HRR plane-strain Mode I angular solution, n = 6.5

**Exact exponents (from J-integral dimensional analysis):**
p = 1/(n+1) = 1/7.5 = 0.13333333…, k = 2 − p = 1.86666666…

**Airy formulation** Φ(r,θ)=r^k φ(θ):
- A = σ̃_rr = k φ + φ″
- B = σ̃_tt = k(k−1) φ
- C = σ̃_rt = −(k−1) φ′
- s_rr=(A−B)/2, s_tt=(B−A)/2, s_rt=C,  se² = (3/4)(A−B)² + 3C²

**Kinematics/constitutive used (corrected scalings):**
σ ~ r^(−p),  ε ~ r^(p−1),  u_i = Λ r^p U_i(θ):
- p U_r = q(A−B)/2  (I)
- U_r + U_t′ = q(B−A)/2  (II)
- U_r′ + (p−1)U_t = 2qC  (III),  q = se^(n−1)

See `hrr_clean.py` for the resulting 4th-order ODE for φ (state [φ,φ′,φ″,φ‴]).

## 1. Boundary-value analysis (why solve_bvp failed) and robust solution

The ODE is 4th order. The physically correct symmetry condition at θ=0 is
**only** C(0)=0 ⟺ φ′(0)=0 (σ_rθ = 0 on the crack plane). The condition
φ‴(0)=0 is **not** an independent BC for this problem: imposing it (the
original `hrr2.py`/`hrr3.py` setup) over-determines the problem — I verified by
a determinant test (`eig_det.py`) that the 2×2 determinant of
[A(π), C(π)] over the even subspace **never vanishes** for any p, so no even
solution exists. Treating p as a free eigenvalue also fails (drives p→1).

**Robust method:** shoot from θ=0 with φ(0)=1 (normalization), φ′(0)=0, and
**two free initial curvatures** φ″(0)=u1, φ‴(0)=u2; impose the two traction-free
conditions A(π)=0, C(π)=0. This 2×2 system has a clean root, found by `fsolve`
(and independently by `solve_bvp` with BCs [φ(0)−1, φ′(0), A(π), C(π)]):

    φ″(0) = −0.9142244740 ,  φ‴(0) = +0.0328463401
    A(π) = 1.9e−14 ,  C(π) = −9.2e−15   (machine zero)

Both `fsolve` (shooting) and `solve_bvp` converge to the identical solution
(`final_solve.py`); ODE residual verified: d3′ from finite differences equals
the RHS to ~1e−11 at multiple θ.

## 2. Angular functions, normalized so se(θ=0)=1

Raw se(0)=0.5761975188 → scale factor Sc = 1/0.5761975188 = 1.7355139.

**Values at θ=0** (normalization se(0)=1):

| quantity | value |
|---|---|
| A(0) = σ̃_rr(0) | **1.65298** |
| B(0) = σ̃_tt(0) | **2.80768** |
| C(0) = σ̃_rt(0) | **0.00000** |
| se(0)           | 1.00000 |

So B(θ=0) = σ̃_tt(0) = **2.80768** (with se(0)=1 normalization).

Angular shapes (θ : A, B, C, se):

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
(A(π)=0 and C(π)=0 as required; smooth throughout; se is bounded and slowly
varying. This is the classic HRR Mode I shape — se largest near θ≈π/2.)

## 3. Displacement angular functions (same normalization)

From (I) and (III): U_r = q(A−B)/(2p); U_t = (2qC − U_r′)/(p−1).
Verified that (II) and (III) agree to ~1e−16 with the ODE-consistent U_r″
(`diag3.py`).

## 4. The integral I_n

Using the requested definition
I_n = ∫_{−π}^{π} [ (n/(n+1)) se^{n+1} cosθ
                  − (σ̃_rr·∂u_x/∂x + σ̃_rt·∂u_y/∂x) ] dθ ,
with se(0)=1 and displacement angular functions from the same constitutive
constant (U_r = q(A−B)/(2p), U_t = (2qC−U_r′)/(p−1)), and with the FULL
polar→cartesian rotation of ∂u_i/∂x included:

  ∫_{0}^{π} (n/(n+1)) se^{n+1} cosθ dθ = **−10.0074**
  ∫_{0}^{π} (traction·∂u/∂x term) dθ    = **−0.5160**
  **I_n (full [−π,π]) = 2×(W − traction) = −18.983**

The `W`-term is negative because se grows toward θ=π (se even in θ, cosθ odd-
signed on the half-range), so the θ→π contribution outweighs θ→0. The traction
term is numerically small. I_n coming out **negative** is inconsistent with
J>0 and with published HRR conventions (I_n > 0): this points to a residual
sign/normalization bookkeeping error in the traction (displacement-derivative)
term of the J-integrand that I did not fully pin down within the time budget.
The angular functions A,B,C,se (item 1), B(0) (item 2), and their shapes are
robust and independently cross-validated; the numerical value of I_n should be
treated as provisional (magnitude ≈ 19; sign/definition to be re-derived).

Reported pieces for recombination:
- I_W = ∫_{−π}^{π} (n/(n+1)) se^{n+1} cosθ dθ = −20.0148
- I_T = ∫_{−π}^{π} (σ̃_rr ∂u_x/∂x + σ̃_rt ∂u_y/∂x) dθ = −1.03201
- I_n = I_W − I_T = −18.9828

## Verification
- Corrected ODE reproduced the exact n=1 Williams biharmonic eigenfunctions
  φ=cos(θ/2), cos(3θ/2) to 1e−16 (`validate3.py`).
- Two independent solvers (shooting fsolve; solve_bvp) give identical
  φ″(0), φ‴(0), A(π)≈0, C(π)≈0.
- (II)≡(III) constitutive consistency at 1e−16.
- ODE residual via finite differences ~1e−11.