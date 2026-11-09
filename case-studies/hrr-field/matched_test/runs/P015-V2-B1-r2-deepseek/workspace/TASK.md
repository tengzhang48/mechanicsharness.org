# Task

Solve the problem below. Write your complete answer to `ANSWER.md` in this folder, with these
parts: the requested numbers or decisions; the model you chose and why; any assumption you
changed and why; how you verified the result; and what is missing or uncertain. You may write
and run code in this folder (`python3` with numpy and scipy). Budget: 45 minutes; write
`ANSWER.md` early, update it, and make sure it is complete by 40 minutes. Keep `PROCESS.md`:
one line per step, each starting with a time pasted from `date -u +%H:%M:%S` (never type a time).

Rules: use only this folder and your own knowledge. Do not open files outside this folder, and do
not use web search or web fetch. The session log will be checked.
Never stop or kill a process you did not start. Kill only by the PID you recorded when you started it; never use pkill, killall, or pgrep-based loops (this machine is shared with other runs).

Provided material: `./knowledge/` contains a mechanics topic (crack-tip and notch-tip asymptotic
fields), checks you can run on your own solution, and a knowledge entry for one paper. See
`./knowledge/README.md`. Use it as you see fit; you are free to choose your own method.

---

# Infer a material stiffness from a transverse strain field

Use mm, N, and MPa (= N/mm²). Strains are dimensionless tensor strains; an engineering shear is twice the corresponding tensor shear. A stationary, straight slit lies along x < 0, y = 0 in a homogeneous, isotropic membrane. Its tip is at the origin. Loading is quasistatic, monotonically increasing, symmetric tension across the slit, and proportional at the material points under consideration. The faces are traction free. There is no bending. The two-dimensional constraint is σ_zz = σ_xz = σ_yz = 0, while ε_zz is free to change. Use infinitesimal kinematics and exclude inertia and damage.

For this problem the **complete total-strain law** is

\[
\epsilon_{ij}=
 \frac{1+\nu}{E}\sigma_{ij}
 -\frac{\nu}{E}\sigma_{kk}\delta_{ij}
 +\frac32\left(\frac1T-\frac1E\right)
       \left\langle 1-\frac{\sigma_y}{\sigma_e}\right\rangle s_{ij},
\]

\[
s_{ij}=\sigma_{ij}-\tfrac13\sigma_{kk}\delta_{ij},\qquad
\sigma_e=\sqrt{\tfrac32 s_{ij}s_{ij}},\qquad
\langle q\rangle=\max(q,0).
\]

Repeated indices sum over x, y, z, including both off-diagonal tensor entries. The last term is zero at σ_e = 0. The constants are E = 70,000 MPa, ν = 0.30, and σ_y = 10.0 MPa. The one unknown material parameter is T, in MPa, with 0 < T < E. There is no unspecified hardening exponent or additional material parameter.

An external calculation supplies J = 0.800 N/mm. Its convention is

\[
J=\int_\Gamma (W N_x-t_i u_{i,x})\,ds,\qquad
t_i=\sigma_{ij}N_j,\qquad W=\int \sigma_{ij}\,d\epsilon_{ij},
\]

where Γ is a circular arc from the lower to the upper face, N is the radial outward normal, ds is positive arc length, and work on the faces vanishes. Along the ray x > 0, a computation with this constitutive law reports the leading coefficient

\[
C_{zz}(0)=\lim_{r\to0}\sqrt r\,\epsilon_{zz}(r,0)
       =-0.005000\ {\rm mm}^{1/2}.
\]

The limit defines a coefficient of the formal local field in this small-strain model; it is not a claim that arbitrarily small radii are physically accessible. Matching to the surrounding body is known to make corrections to this leading field negligible for r ≤ 0.200 mm.

Determine:

1. T, in MPa.
2. At θ = 60°, the three coefficients C_rr = lim √r ε_rr, C_θθ = lim √r ε_θθ, and C_γrθ = lim √r (2ε_rθ), all in mm^1/2. Here x = r cos θ, y = r sin θ, e_r = (cos θ,sin θ), and e_θ = (−sin θ,cos θ).
3. C_y(90°) = lim u_y(r,90°)/√r, in mm^1/2, for the same leading field. Displacements have no additive rigid translation.
4. The interval of radii **on θ = 0** that passes all three of the following screening conditions, evaluated using the leading field: σ_e ≥ 10σ_y; the largest absolute value of the three principal strains is ≤ 0.0500; and r ≤ 0.200 mm. Report both endpoints in mm, or demonstrate that the interval is empty. These screens are defined checks, not error estimates for the full finite-r solution.

Justify the field and the constitutive reduction you use, including the treatment of ε_zz. Include one quantitative verification that is not just substitution into the measured C_zz. Explain briefly what these screens establish and whether the supplied constitutive and loading data also determine crack initiation.

T has a **1% relative tolerance**. Each of the other six requested numbers has a **2% relative tolerance**. Any correct formulation earns full credit. This is intended to be solvable in 45 minutes using Python with numpy and scipy.
