# Task

Solve the problem below. Budget: 170 minutes. Write `ANSWER.md` early and keep it updated, and
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

# Problem: stress relaxation by diffusion into a grain boundary of a thin film

A polycrystalline film of thickness h lies on a substrate with the same isotropic elastic constants
(Young's modulus E, Poisson's ratio ν, E* = E/(1 − ν²)), so the two together form an elastic
half-space in plane strain. The film's free surface is z = 0. A grain boundary runs perpendicular to
the surface from z = 0 to the film/substrate interface at z = h. Initially the normal traction on the
grain boundary is uniform, σ₀ > 0 (tension).

Atoms diffuse from the free surface into the grain boundary and plate onto it. This inserts a wedge
of material with opening 2u(z, t), equivalent to climb-edge dislocations on the boundary. Neither
diffusion nor sliding occurs at the film/substrate interface, and the free surface stays flat
(surface diffusion is fast; grooving is neglected).

- **Elasticity.** The traction on the boundary due to the wedge is
  σ_gb(z,t) = σ₀ − (E*/2π) ∫₀^h S(z,ζ) ∂u(ζ,t)/∂ζ dζ, with
  S(z,ζ) = 1/(z − ζ) − 1/(z + ζ) − 2ζ(z − ζ)/(z + ζ)³ (principal value at ζ = z).
- **Diffusion.** The chemical potential is μ = −σ_gb Ω. The atomic flux along the boundary is
  j = −(δ_gb D_gb/kTΩ) ∂μ/∂z. Mass conservation gives 2 ∂u/∂t = −Ω ∂j/∂z. Here Ω is the atomic
  volume and δ_gb D_gb the boundary thickness times its diffusivity.
- **Boundary and initial conditions:** σ_gb(0,t) = 0; ∂σ_gb/∂z = ∂²σ_gb/∂z² = 0 at z = h; and
  σ_gb(z,0) = σ₀.

Use t₀ = 4πkT h³/(E* δ_gb D_gb Ω) as the unit of time.

(a) Describe how σ_gb(z,t) evolves. Give the rate at which the traction decays at late times, as a
dimensionless number λ₁ (the decay ∝ e^(−λ₁t/t₀)), and the next rate, λ₂.

(b) What does the grain boundary become as t → ∞? Give the stress intensity factor at z = h in that
limit, in units of σ₀√(πh).

(c) On what time scale, in units of t₀, does the traction relax?

**Report in `results.json`:**
- `lambda1`, `lambda2`;
- `mode1_shape_half`: the normalized traction profile of the slowest-decaying part, f(z/h) with
  f(1) = 1, evaluated at z/h = 0.5;
- `K_inf_norm`: K(t → ∞)/(σ₀√(πh));
- `relaxation_time_t0`: your relaxation time scale in units of t₀;
- `method`: the method you used, as a string.
