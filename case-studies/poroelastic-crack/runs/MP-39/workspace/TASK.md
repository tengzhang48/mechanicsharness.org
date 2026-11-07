# Task

Solve the problem below. Budget: 50 minutes. Write `ANSWER.md` early and keep it updated, and
make sure it is complete before the budget ends. You may write and run code in this folder
(`python3` with numpy and scipy).

**Deliver three things:**

1. **`ANSWER.md`.** Your answers, the model and method you chose and why, your derivation in brief,
   and what is missing or uncertain. If anything in the problem is unclear or missing, list it under
   "Questions" with the assumption you used, and continue.
2. **`results.json`.** It holds every quantity listed under "Report" in the problem, with the names
   given there (numbers as numbers, formulas as strings), plus `reproduce_command`: the command of
   item 3.
3. **One command that reproduces your numbers,** named in ANSWER.md (for example `python3 solve.py`).
   It must print the values in `results.json`.

**For every check you make:** report what you compared, the output (pasted, not retyped), and how the
check could have failed.

Keep `PROCESS.md`: one line per step, each starting with a time pasted from `date -u +%H:%M:%S`
(never type a time).

Rules: use only this folder and your own knowledge. Do not open files outside this folder, and do
not use web search or web fetch. The session log will be checked. Never stop or kill a process you
did not start. Kill only by the PID you recorded when you started it; never use pkill, killall, or
pgrep-based loops (this machine is shared with other runs).

---

# Problem: a suddenly pressurized crack in a fluid-saturated porous solid

An infinite, homogeneous, isotropic, fluid-saturated porous solid deforms in plane strain
(x–y plane). A crack occupies −L ≤ x ≤ L, y = 0. The solid obeys Biot's linear poroelasticity, with
shear modulus G, drained Poisson's ratio ν, undrained Poisson's ratio ν_u (ν < ν_u ≤ ½), Biot
coefficient α (0 < α ≤ 1) and permeability coefficient κ.

In terms of the total stress σ_ij (tension positive), the pore pressure p, the strain e_ij and the
variation of fluid content per unit volume ζ:

- 2G e_ij = σ_ij − ν σ_kk δ_ij + α(1 − 2ν) δ_ij p
- 2G ζ = α(1 − 2ν) σ_kk + α²(1 − 2ν)² p / (ν_u − ν)
- equilibrium: σ_ij,j = 0
- Darcy: q_i = −κ p,i
- mass balance: ∂ζ/∂t + q_i,i = 0

(indices 1, 2; summation implied). The diffusivity is c = 2κG(1 − ν)(ν_u − ν) / [α²(1 − 2ν)²(1 − ν_u)].

Stresses and pore pressure are zero everywhere for t < 0 and vanish at infinity. Two loadings act on
both crack faces (−L ≤ x ≤ L) for t > 0. Here σ_n = σ_yy is the normal stress on the faces and H is
the unit step.

- **Loading 1:** σ_n = −H(t), p = 0.
- **Loading 2:** σ_n = 0, p = H(t).

Let K be the mode I stress intensity factor at the crack tips, and V the crack volume per unit
thickness (the integral of the opening over the crack; positive for an open crack).

(a) For loading 1, find K and V just after loading (t = 0⁺) and at long times (t → ∞). Find also the
leading-order growth of V at small times, ct/L² ≪ 1.

(b) For loading 2, find K at t = 0⁺ and as t → ∞, and V as t → ∞.

(c) Evaluate your results for ν = 0.2, ν_u = 0.4 and α = 0.89.

**Report in `results.json`** (for (c); K in units of √(πL), V in units of L²/G, per unit loading
magnitude):
- `K1_0`, `K1_inf`, `V1_0`, `V1_inf`: loading 1 at t = 0⁺ and t → ∞;
- `V1_short_coeff`: the coefficient A in V ≈ V1_0 + A √(ct/L²) for loading 1 at small times;
- `p_faces_0`: in loading 1, the pore pressure induced in the solid immediately next to the crack
  faces at t = 0⁺, before it diffuses;
- `K2_0`, `K2_inf`, `V2_inf`: loading 2;
- `formulas`: your general expressions for all of the above, as strings.
