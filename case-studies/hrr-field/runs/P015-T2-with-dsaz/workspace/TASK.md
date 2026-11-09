# Task

Solve the problem below. Write your complete answer to `ANSWER.md` in this folder: the requested
numbers or decisions, your method, and your checks. You may write and run code in this folder
(`python3` with numpy and scipy). Budget: 45 minutes; write `ANSWER.md` early, update it, and make
sure it is complete by 40 minutes. Keep `PROCESS.md`: one line per step, each starting with a time
pasted from `date -u +%H:%M:%S` (never type a time).

Rules: use only this folder (including `./knowledge/`) and your own knowledge. Do not open files
outside this folder, and do not use web search or web fetch. The session log will be checked.
Never stop or kill a process you did not start. Kill only by the PID you recorded when you started it; never use pkill, killall, or pgrep-based loops (this machine is shared with other runs).

Provided material: `./knowledge/` contains a knowledge unit about a mechanics paper (an entry with a
"Using this unit" section, an extraction, a report, Python code and its tests). Use it as you see fit.

---

# A tensile notch

A homogeneous body contains a sharp, traction-free notch. In polar coordinates centered at its tip, the **material occupies −2π/3 < θ < 2π/3**; the missing wedge is centered on the negative x-axis and has a 120° opening. θ = 0 points into the material along the symmetry line. The body is in plane strain. Remote tension acts in the y direction, increases proportionally, and excites the leading symmetric opening field with nonzero amplitude.

The material is isotropic and pressure-insensitive. Under monotonic proportional loading, its small-strain response is

\[
\epsilon_{ij}=\frac{1+\nu}{E}\sigma_{ij}
 -\frac{\nu}{E}\sigma_{kk}\delta_{ij}
 +\frac32\alpha\frac{\sigma_0}{E}
 \left(\frac{\sigma_e}{\sigma_0}\right)^{n-1}\frac{s_{ij}}{\sigma_0},
\quad s_{ij}=\sigma_{ij}-\tfrac13\sigma_{kk}\delta_{ij},
\quad \sigma_e=\sqrt{\tfrac32s_{ij}s_{ij}}.
\]

Use E = 200 GPa, ν = 0.30, σ₀ = 300 MPa, α = 0.02, and n = **4**. Indices run over three spatial directions, repeated indices are summed, and shear strain is tensor shear. Plane strain imposes εzz = εrz = εθz = 0. Treat the sharp tip within this small-strain continuum idealization, with finite energy in any finite neighborhood.

For the **leading local asymptote as r → 0**, find:

1. p in σθθ(r,0) ∝ r⁻ᵖ, with absolute tolerance **0.0003**.
2. R = σrr(r,0)/σθθ(r,0), with absolute tolerance **0.002**. Both stresses are at the same radius on θ = 0, and the opening amplitude is chosen positive.

Give the two numbers, a short physical justification of your model and boundary conditions, and one quantitative verification. The unspecified remote-load magnitude cancels from these two quantities. Allow about 45 minutes with NumPy/SciPy.

General computational identities, if useful: an Airy function Φ gives σrr = Φ,r/r + Φ,θθ/r², σθθ = Φ,rr, σrθ = −∂r(Φ,θ/r). The strain-displacement relations are εrr = ur,r, εθθ = (ur + uθ,θ)/r, and 2εrθ = ur,θ/r + uθ,r − uθ/r.
