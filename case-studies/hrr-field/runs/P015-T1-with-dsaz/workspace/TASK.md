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

# A noninteger hardening exponent

A stationary, sharp crack occupies the negative x-axis in a homogeneous body in **plane strain**, with traction-free faces at θ = ±π. Polar coordinates are centered at its tip; θ = 0 points ahead of it. Loading is symmetric opening, increases proportionally, and has positive energy flux J. Use small-strain kinematics and the following idealized, fully plastic total-deformation law:

\[
\epsilon_{ij}=\frac32\alpha\epsilon_0
 \left(\frac{\sigma_e}{\sigma_0}\right)^{n-1}\frac{s_{ij}}{\sigma_0},
\qquad s_{ij}=\sigma_{ij}-\frac13\sigma_{kk}\delta_{ij},
\qquad \sigma_e=\sqrt{\frac32s_{ij}s_{ij}}.
\]

Here n = **6.5**, α = 0.02, σ₀ = 400 MPa, E = 200 GPa, and ε₀ = σ₀/E. E only sets the reference strain; do not add an elastic strain term to this limiting law. Indices run over all three spatial directions, repeated indices are summed, and εᵣθ is tensor shear strain. Plane strain means εzz = εrz = εθz = 0.

Let W = n σᵢⱼεᵢⱼ/(n+1) be strain energy per volume. On a counterclockwise circular contour, with outward radial normal N and traction tᵢ = σᵢⱼNⱼ, define

\[
J=\int_{-\pi}^{\pi}\left(W\cos\theta-t_i\frac{\partial u_i}{\partial x}\right)r\,d\theta.
\]

Crack-face contributions vanish. J has units J/m². Let U(R) = ∫₀ᴿ∫₋π^π W r dθ dr, the energy inside radius R **per unit out-of-plane thickness**, in J/m. Work with the leading tip field and take R within its idealized range.

Compute these three dimensionless numbers:

1. p in σθθ(r,0) ∝ r⁻ᵖ, to absolute tolerance **10⁻⁵**.
2. H defined by σθθ(r,0) = H σ₀ [J/(α ε₀ σ₀ r)]ᵖ, to absolute tolerance **0.005**.
3. C = U(R)/(J R), to absolute tolerance **0.002**.

Submit the numbers, a brief method, and one verification. No particular code format is required; allow about 45 minutes using NumPy/SciPy.

Computational identities, if useful: an Airy function Φ gives σrr = Φ,r/r + Φ,θθ/r², σθθ = Φ,rr, σrθ = −∂r(Φ,θ/r). In polar components, εrr = ur,r, εθθ = (ur + uθ,θ)/r, and 2εrθ = ur,θ/r + uθ,r − uθ/r. Cartesian differentiation of displacement must include rotation of the polar basis.
